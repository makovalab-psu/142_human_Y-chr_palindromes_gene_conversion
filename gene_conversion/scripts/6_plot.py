#!/usr/bin/env python3
"""Generate per-site phylogeny plots (Graphviz PDF) for each palindrome."""
import os
import subprocess
from collections import defaultdict
from multiprocessing import Pool

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR = "work"
NCPUS = 6

COLORS = {"0/0": "#F5E9CF", "0/1": "#7DB9B6", "1/1": "#E96479"}

DOTFILE_TEMPLATE = "work/{palindrome}/plots/temp_{pos}.dot"


def gene_conversion_or_mutation(ancestral: str, new: str) -> str:
    if ancestral in ("0/0", "1/1") and new == "0/1":
        return "MUT"
    if ancestral == "0/1" and new in ("0/0", "1/1"):
        return "GC"
    return "?"


def traverse_tree_graph(node: list, ancestral: str, ancestral_name: str, out) -> None:
    if len(node[4]) != 1 and ancestral in node[4]:
        new_node = [ancestral]
    elif len(node[4]) != 1 and ancestral not in node[4]:
        new_node = ["0/1"]
    else:
        new_node = node[4]

    pretty = ",".join(new_node)
    note_info = f'"{node[0]}\\n{pretty}"'
    anc_info = f'"{ancestral_name}\\n{ancestral}"'

    if ancestral != new_node[0]:
        label = gene_conversion_or_mutation(ancestral, new_node[0])
        if anc_info != note_info:
            out.write(f"{anc_info} -> {note_info} [label=<<B>{label}</B>>]\n")
    else:
        if anc_info != note_info:
            out.write(f"{anc_info} -> {note_info}\n")

    if node[2]:
        for child in node[1]:
            traverse_tree_graph(child, new_node[0], node[0], out)
    else:
        color = COLORS.get(pretty, "white")
        out.write(f'{note_info} [style=filled fillcolor="{color}" shape="box"]\n')


def traverse_tree(node: list, ancestral: str, ancestral_name: str, results: list) -> list:
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
            traverse_tree(child, new_node[0], node[0], results)

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
    elif (len(node1[4]) != 1) != (len(node2[4]) != 1):
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


def plot_palindrome(palindrome: str) -> None:
    vcf_file = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    tree_file = os.path.join(WORK_DIR, palindrome, "tree", f"{palindrome}.tree.txt")
    plots_dir = os.path.join(WORK_DIR, palindrome, "plots")

    stats: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    individuals: list = []

    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                individuals = line.strip().split()[9:]
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, _, _, *genotypes = line.strip().split()
            for gt, ind in zip(genotypes, individuals):
                gt = gt.split(":")[0].replace("|", "/")
                if gt == "1/0":
                    gt = "0/1"
                stats[chrom][pos][ind] = gt

    n_plotted = 0
    for palindrome_id in stats:
        for pos in sorted(stats[palindrome_id].keys(), key=int):
            genotypes = [stats[palindrome_id][pos].get(x, "./.") for x in individuals]
            root, nodes = initialize_nodes(individuals, genotypes, tree_file)

            lowest_mutations = 999
            best_genotype = ""
            for genotype in nodes[root][4]:
                results: list = []
                results = traverse_tree(nodes[root], genotype, nodes[root][0], results)
                mutations = sum(1 for e in results if e[2] == "MUT")
                if mutations < lowest_mutations:
                    lowest_mutations = mutations
                    best_genotype = genotype

            if lowest_mutations >= 2 or not best_genotype:
                continue

            dot_file = os.path.join(plots_dir, f"temp_{pos}.dot")
            pdf_file = os.path.join(plots_dir, f"{palindrome_id}_{pos}.pdf")

            with open(dot_file, "w") as out:
                out.write('digraph {\nfontsize="10"\n')
                traverse_tree_graph(nodes[root], best_genotype, nodes[root][0], out)
                out.write("}\n")

            subprocess.run(
                ["dot", "-T", "pdf", dot_file, "-Nshape=rect", "-Gsize=10,10", "-o", pdf_file],
                check=False,
            )
            os.remove(dot_file)
            n_plotted += 1

    print(f"{palindrome}: {n_plotted} plots → {plots_dir}")


def main() -> None:
    with Pool(NCPUS) as pool:
        pool.map(plot_palindrome, PALINDROMES)


if __name__ == "__main__":
    main()
