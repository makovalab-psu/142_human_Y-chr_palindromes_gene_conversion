#!/usr/bin/env python3
"""Parse Nexus tree, prune to per-palindrome samples, output flat topology table."""
import os
import re

NEXUS = "data/143males_HGSVC_HPRC_CEPH_110925_150M.nex"
PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR = "work"


# ── Minimal Newick parser ────────────────────────────────────────────────────

def parse_newick(s: str) -> dict:
    """
    Parse a Newick string into a nested dict:
        {'name': str, 'children': [child_dict, ...]}
    Leaf names are their sample names; internal names are None until assigned.
    """
    s = s.strip().rstrip(";")

    def parse(pos: int) -> tuple[dict, int]:
        children = []
        if s[pos] == "(":
            pos += 1  # consume '('
            while True:
                child, pos = parse(pos)
                children.append(child)
                if s[pos] == ",":
                    pos += 1
                elif s[pos] == ")":
                    pos += 1
                    break
        # Read label (may include branch length like 'name:0.1')
        end = pos
        while end < len(s) and s[end] not in (",", ")", ";"):
            end += 1
        label = s[pos:end]
        name = label.split(":")[0].strip() if label else None
        pos = end
        node = {"name": name or None, "children": children}
        return node, pos

    root, _ = parse(0)
    return root


def get_leaves(node: dict) -> set[str]:
    if not node["children"]:
        return {node["name"]}
    leaves = set()
    for child in node["children"]:
        leaves |= get_leaves(child)
    return leaves


def prune(node: dict, keep: set[str]) -> dict | None:
    """
    Return a pruned copy of the tree containing only tips in keep.
    Returns None if this subtree has no kept tips.
    """
    if not node["children"]:
        return {"name": node["name"], "children": []} if node["name"] in keep else None

    kept_children = [pruned for c in node["children"] if (pruned := prune(c, keep)) is not None]

    if not kept_children:
        return None
    if len(kept_children) == 1:
        # Collapse single-child node — just return the child
        return kept_children[0]

    return {"name": node["name"], "children": kept_children}


def name_internal_nodes(node: dict, counter: list) -> None:
    """Assign node1, node2, ... to internal nodes in preorder (level-like)."""
    if node["children"]:
        counter[0] += 1
        node["name"] = f"node{counter[0]}"
        for child in node["children"]:
            name_internal_nodes(child, counter)


def write_topology(node: dict, out) -> None:
    """Write parent \\t child1 \\t child2 in postorder (children before parent)."""
    if node["children"]:
        for child in node["children"]:
            write_topology(child, out)
        if len(node["children"]) == 2:
            out.write(f"{node['name']}\t{node['children'][0]['name']}\t{node['children'][1]['name']}\n")


def count_tips_and_internals(node: dict) -> tuple[int, int]:
    if not node["children"]:
        return 1, 0
    tips, internals = 0, 1
    for child in node["children"]:
        t, i = count_tips_and_internals(child)
        tips += t
        internals += i
    return tips, internals


# ── Main ─────────────────────────────────────────────────────────────────────

def load_newick(nexus_path: str) -> str:
    with open(nexus_path) as f:
        content = f.read()
    return re.search(r"=\s*\[&R\]\s*(.+);", content).group(1).strip()


def get_vcf_samples(palindrome: str) -> list[str]:
    vcf = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    with open(vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.strip().split("\t")[9:]
    return []


def main() -> None:
    nwk = load_newick(NEXUS)
    full_tree = parse_newick(nwk)
    all_tips = get_leaves(full_tree)

    for palindrome in PALINDROMES:
        samples = get_vcf_samples(palindrome)
        keep = set(samples) & all_tips

        missing_from_tree = set(samples) - all_tips
        if missing_from_tree:
            print(f"{palindrome}: {len(missing_from_tree)} VCF samples not in tree: {missing_from_tree}")

        pruned = prune(full_tree, keep)
        name_internal_nodes(pruned, [0])

        out_path = os.path.join(WORK_DIR, palindrome, "tree", f"{palindrome}.tree.txt")
        with open(out_path, "w") as f:
            write_topology(pruned, f)

        n_tips, n_internal = count_tips_and_internals(pruned)
        print(f"{palindrome}: {n_tips} tips, {n_internal} internal nodes → {out_path}")


if __name__ == "__main__":
    main()
