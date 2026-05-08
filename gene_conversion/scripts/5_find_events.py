#!/usr/bin/env python3
"""Parsimony-based classification of mutations vs gene conversion events per palindrome."""
import os
from collections import defaultdict
from multiprocessing import Pool

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR = "work"
NCPUS = 6  # one per palindrome

TRANS_MATRIX = {
    "A": {"A": "None", "C": "Transversion", "G": "Transition", "T": "Transversion"},
    "C": {"A": "Transversion", "C": "None", "G": "Transversion", "T": "Transition"},
    "G": {"A": "Transition", "C": "Transversion", "G": "None", "T": "Transversion"},
    "T": {"A": "Transversion", "C": "Transition", "G": "Transversion", "T": "None"},
}


def gene_conversion_or_mutation(ancestral: str, new: str) -> str:
    if ancestral in ("0/0", "1/1") and new == "0/1":
        return "mutation"
    if ancestral == "0/1" and new in ("0/0", "1/1"):
        return "geneconversion"
    return "unknown"


def traverse_tree(node: list, ancestral: str, results: list) -> list:
    if len(node[4]) != 1 and ancestral in node[4]:
        new_node = [ancestral]
    elif len(node[4]) != 1 and ancestral not in node[4]:
        new_node = ["0/1"]
    else:
        new_node = node[4]

    if ancestral != new_node[0]:
        results.append([ancestral, new_node[0], gene_conversion_or_mutation(ancestral, new_node[0]), node[0], node[3]])

    if node[2]:
        for child in node[1]:
            traverse_tree(child, new_node[0], results)

    return results


def make_node(node1: list, node2: list, name: str) -> list:
    genotype_counter: dict = defaultdict(int)
    genotype: list = []

    if len(node1[4]) == 1 and len(node2[4]) == 1:
        if node1[4] == node2[4]:
            genotype = node1[4]
        else:
            combined = node1[4] + node2[4]
            if "0/0" in combined and "1/1" in combined:
                genotype = ["0/1"]
            else:
                genotype = combined
    elif (len(node1[4]) != 1) != (len(node2[4]) != 1):  # XOR: one ambiguous
        for g in node1[4] + node2[4]:
            genotype_counter[g] += 1
        if len(genotype_counter) == 3:
            genotype = ["0/1"]
        else:
            genotype = [max(genotype_counter, key=genotype_counter.get)]
    else:
        genotype = list(set(node1[4] + node2[4]))

    return [name, [node1, node2], True, node1[3] + node2[3], genotype]


def initialize_nodes(names: list, genotypes: list, topology_file: str) -> tuple[str, dict]:
    nodes: dict = {}
    for name, geno in zip(names, genotypes):
        nodes[name] = [name, None, False, [name], [geno]]

    has_parents: dict = {}
    with open(topology_file) as f:
        for line in f:
            nodename, child1, child2 = line.strip().split("\t")
            nodes[nodename] = make_node(nodes[child1], nodes[child2], nodename)
            has_parents[child1] = nodename
            has_parents[child2] = nodename

    root = next(k for k in nodes if k not in has_parents)
    return root, nodes


def conversion_details(
    a: str, b: str, ref: str, alt: str, ancestral: str
) -> tuple[str, str, str, str]:
    nuc = [ref, alt]

    ancestral_base = "het"
    if ancestral in ("0/0", "1/1"):
        ancestral_base = nuc[int(ancestral[0])]

    if a == "0/1":  # gene conversion
        if b == "0/0":
            converting_from, converting_to = alt, ref
            if ancestral == "0/0":
                anc_der = "derived->ancestral"
            elif ancestral == "1/1":
                anc_der = "ancestral->derived"
            else:
                anc_der = "ancestral_heterozygote"
        else:  # b == "1/1"
            converting_from, converting_to = ref, alt
            if ancestral == "0/0":
                anc_der = "ancestral->derived"
            elif ancestral == "1/1":
                anc_der = "derived->ancestral"
            else:
                anc_der = "ancestral_heterozygote"
    else:  # mutation
        if a == "0/0":
            converting_from, converting_to = ref, alt
        else:
            converting_from, converting_to = alt, ref
        anc_der = "ancestral->derived"

    nuc_change = f"{converting_from}->{converting_to}"
    ts_tv = TRANS_MATRIX.get(converting_from, {}).get(converting_to, "Unknown")
    return nuc_change, anc_der, ts_tv, ancestral_base


def process_palindrome(palindrome: str) -> None:
    vcf_file = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    tree_file = os.path.join(WORK_DIR, palindrome, "tree", f"{palindrome}.tree.txt")
    out_file = os.path.join(WORK_DIR, palindrome, "events", f"{palindrome}.events.txt")

    stats: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    ref_dict: dict = defaultdict(dict)
    alt_dict: dict = defaultdict(dict)
    individuals: list = []

    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                individuals = line.strip().split()[9:]
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, _, _, *genotypes = line.strip().split()
            ref_dict[chrom][pos] = ref
            alt_dict[chrom][pos] = alt
            for gt, ind in zip(genotypes, individuals):
                gt = gt.split(":")[0].replace("|", "/")
                if gt == "1/0":
                    gt = "0/1"
                stats[chrom][pos][ind] = gt

    with open(out_file, "w") as out:
        out.write(
            "palindrome\tpos\tfrom\tto\tbase_change\tancestral\ttype\t"
            "node\tindividuals\tancestral_state\tTransition_transversion\tfilter\tancestral_base\n"
        )

        for palindrome_id in stats:
            for pos in sorted(stats[palindrome_id].keys(), key=int):
                genotypes = [stats[palindrome_id][pos].get(x, "./.") for x in individuals]
                ref = ref_dict[palindrome_id][pos]
                alt = alt_dict[palindrome_id][pos]

                root, nodes = initialize_nodes(individuals, genotypes, tree_file)

                best_results: list = []
                lowest_mutations = 999
                best_genotype = ""

                for genotype in nodes[root][4]:
                    results: list = []
                    mutations = 0
                    results = traverse_tree(nodes[root], genotype, results)
                    for event in results:
                        if event[2] == "mutation":
                            mutations += 1
                    if mutations < lowest_mutations:
                        lowest_mutations = mutations
                        best_results = results
                        best_genotype = genotype

                if lowest_mutations >= 2:
                    continue

                for event in best_results:
                    nuc_change, anc_der, ts_tv, anc_base = conversion_details(
                        event[0], event[1], ref, alt, best_genotype
                    )
                    out.write(
                        "\t".join([
                            palindrome_id, pos,
                            event[0], event[1],
                            nuc_change, anc_der, event[2],
                            event[3], ",".join(event[4]),
                            best_genotype, ts_tv, "pass", anc_base,
                        ]) + "\n"
                    )

    n_events = sum(1 for line in open(out_file) if not line.startswith("palindrome"))
    print(f"{palindrome}: {n_events} events → {out_file}")


def main() -> None:
    with Pool(NCPUS) as pool:
        pool.map(process_palindrome, PALINDROMES)


if __name__ == "__main__":
    main()
