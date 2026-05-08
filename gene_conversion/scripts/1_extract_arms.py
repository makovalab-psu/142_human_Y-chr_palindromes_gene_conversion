#!/usr/bin/env python3
"""Extract palindrome arm sequences from masked FASTA files."""
import os
import subprocess
from multiprocessing import Pool
from collections import defaultdict

BED = "data/output_palindromes20260505_145745.bed"
FASTA_DIR = "data/masked/masked_sequences"
# PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
# PALINDROMES = ["P1Long","P2Long"]
PALINDROMES = ["RED"]
BED_NAME_MAP = {"Q6": "P9"}   # BED file uses old names; map to canonical names
EXCLUDE_SAMPLES = {"HG03456"}
WORK_DIR = "work"
NCPUS = 10


def parse_bed(bed_file: str) -> tuple[dict, dict[str, set]]:
    """
    Return ({palindrome: {sample: {'A': ..., 'B': ...}}}, {palindrome: {duplicate_samples}}).
    Samples with more than one arm-A or arm-B entry for the same palindrome are flagged as
    duplicates (they have two copies of that palindrome) and must be excluded.
    """
    raw: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    with open(bed_file) as f:
        for line in f:
            chrom, start, end, name = line.strip().split("\t")
            parts = name.split("|")
            palindrome = BED_NAME_MAP.get(parts[-1], parts[-1])
            if palindrome not in PALINDROMES:
                continue
            arm = "A" if parts[0].endswith("A") else "B"
            sample = chrom.split("_chrY")[0]
            raw[palindrome][sample][arm].append((chrom, int(start), int(end)))

    data: dict = defaultdict(lambda: defaultdict(dict))
    duplicates: dict[str, set] = defaultdict(set)
    for palindrome, samples in raw.items():
        for sample, arms in samples.items():
            if any(len(entries) > 1 for entries in arms.values()):
                duplicates[palindrome].add(sample)
                continue
            for arm, entries in arms.items():
                data[palindrome][sample][arm] = entries[0]
    return data, duplicates


def extract_arm(args: tuple) -> None:
    sample, chrom, start, end, out_fa, fasta_dir = args
    if os.path.exists(out_fa):
        return
    fa = os.path.join(fasta_dir, f"{sample}_masked.fa")
    region = f"{chrom}:{start + 1}-{end}"  # samtools is 1-based
    result = subprocess.run(
        ["samtools", "faidx", fa, region],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.splitlines(keepends=True)
    lines[0] = f">{sample}\n"
    with open(out_fa, "w") as f:
        f.writelines(lines)


def main() -> None:
    all_samples = {
        os.path.basename(f).replace("_masked.fa", "")
        for f in os.listdir(FASTA_DIR)
        if f.endswith(".fa") and not f.endswith(".fai")
    }

    data, duplicates = parse_bed(BED)
    tasks: list[tuple] = []

    for palindrome in PALINDROMES:
        arms_dir = os.path.join(WORK_DIR, palindrome, "arms")
        os.makedirs(arms_dir, exist_ok=True)
        os.makedirs(os.path.join(WORK_DIR, palindrome, "variants"), exist_ok=True)
        os.makedirs(os.path.join(WORK_DIR, palindrome, "vcf"), exist_ok=True)
        os.makedirs(os.path.join(WORK_DIR, palindrome, "tree"), exist_ok=True)
        os.makedirs(os.path.join(WORK_DIR, palindrome, "events"), exist_ok=True)
        os.makedirs(os.path.join(WORK_DIR, palindrome, "plots"), exist_ok=True)

        exclude = EXCLUDE_SAMPLES | duplicates.get(palindrome, set())
        if duplicates.get(palindrome):
            print(f"{palindrome}: excluding {sorted(duplicates[palindrome])} (multiple palindrome copies)")

        present = set(data[palindrome].keys())
        missing = all_samples - present

        with open(os.path.join(WORK_DIR, palindrome, "missing.txt"), "w") as f:
            for s in sorted(missing):
                f.write(s + "\n")

        incomplete = []
        for sample, arms in data[palindrome].items():
            if sample in exclude:
                continue
            if "A" not in arms or "B" not in arms:
                incomplete.append(sample)
                continue
            for arm in ("A", "B"):
                chrom, start, end = arms[arm]
                out_fa = os.path.join(arms_dir, f"{sample}_{arm}.fa")
                tasks.append((sample, chrom, start, end, out_fa, FASTA_DIR))

        if incomplete:
            print(f"{palindrome}: {len(incomplete)} samples with incomplete arm pairs — skipped")

    print(f"Extracting {len(tasks)} arm sequences with {NCPUS} workers...")
    with Pool(NCPUS) as pool:
        pool.map(extract_arm, tasks)
    print("Done.")


if __name__ == "__main__":
    main()
