#!/usr/bin/env python3
"""Build VCF from per-arm variant calls. Genotype = arm A vs arm B per sample."""
import os
import pickle
from collections import Counter

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
EXCLUDE_SAMPLES = {"HG03456"}
WORK_DIR = "work"
MIN_SAMPLE_FRACTION = 0.80  # site must be genotypeable in >= this fraction of samples
MULTIALLELIC_LOG = True


def build_vcf(palindrome: str) -> None:
    pkl = os.path.join(WORK_DIR, palindrome, "variants", "variants.pkl")
    with open(pkl, "rb") as f:
        data: dict = pickle.load(f)  # {sample: {'A': {pos: (ref,alt)}, 'B': {pos: (ref,alt)}}}

    # Keep only samples with both arms, excluding blacklisted samples
    samples = sorted(s for s, arms in data.items() if "A" in arms and "B" in arms and s not in EXCLUDE_SAMPLES)
    min_count = int(len(samples) * MIN_SAMPLE_FRACTION)

    # Collect all variant positions and alleles
    # At each position: gather all alleles from all arms across all samples
    pos_alleles: dict[int, Counter] = {}
    pos_ref: dict[int, str] = {}  # reference allele (from CS tag, same across all records)

    for sample in samples:
        for arm in ("A", "B"):
            for pos, (ref_allele, query_allele) in data[sample][arm].items():
                if pos not in pos_alleles:
                    pos_alleles[pos] = Counter()
                    pos_ref[pos] = ref_allele
                pos_alleles[pos][query_allele] += 1

    # For each position, the reference allele itself is carried by arms with no variant
    # Add implicit reference allele counts
    for pos in pos_alleles:
        ref_allele = pos_ref[pos]
        n_arms_with_variant = sum(pos_alleles[pos].values())
        n_arms_total = len(samples) * 2
        pos_alleles[pos][ref_allele] += n_arms_total - n_arms_with_variant

    multiallelic_sites = []
    vcf_records = []

    for pos in sorted(pos_alleles):
        allele_counts = pos_alleles[pos]
        unique_alleles = [a for a, _ in allele_counts.most_common()]

        if len(unique_alleles) > 2:
            multiallelic_sites.append((pos, dict(allele_counts)))
            continue

        if len(unique_alleles) < 2:
            continue  # monomorphic

        # REF = majority allele, ALT = minority
        ref_allele, alt_allele = unique_alleles[0], unique_alleles[1]

        # Build genotypes
        genotypes = []
        n_genotyped = 0
        for sample in samples:
            arm_a = data[sample]["A"]
            arm_b = data[sample]["B"]
            a_allele = arm_a[pos][1] if pos in arm_a else pos_ref[pos]
            b_allele = arm_b[pos][1] if pos in arm_b else pos_ref[pos]

            def to_gt(allele: str) -> str:
                if allele == ref_allele:
                    return "0"
                elif allele == alt_allele:
                    return "1"
                return "."

            gt_a, gt_b = to_gt(a_allele), to_gt(b_allele)
            if "." not in (gt_a, gt_b):
                n_genotyped += 1
            genotypes.append("/".join(sorted([gt_a, gt_b])))

        if n_genotyped < min_count:
            continue

        vcf_records.append((pos + 1, ref_allele, alt_allele, genotypes))  # 1-based POS

    # Write VCF
    out = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    with open(out, "w") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write("\t".join(samples) + "\n")
        for pos, ref, alt, genotypes in vcf_records:
            f.write(f"{palindrome}\t{pos}\t.\t{ref}\t{alt}\t999\tPASS\t.\tGT\t")
            f.write("\t".join(genotypes) + "\n")

    # Log multiallelic sites
    if MULTIALLELIC_LOG and multiallelic_sites:
        log = os.path.join(WORK_DIR, palindrome, "vcf", "multiallelic.txt")
        with open(log, "w") as f:
            f.write("pos\tallele_counts\n")
            for pos, counts in multiallelic_sites:
                f.write(f"{pos + 1}\t{counts}\n")

    print(
        f"{palindrome}: {len(vcf_records)} sites, {len(samples)} samples"
        + (f", {len(multiallelic_sites)} multiallelic skipped" if multiallelic_sites else "")
    )


def main() -> None:
    for palindrome in PALINDROMES:
        build_vcf(palindrome)


if __name__ == "__main__":
    main()
