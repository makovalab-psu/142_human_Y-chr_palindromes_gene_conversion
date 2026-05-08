#!/usr/bin/env python3
"""Plot P3 node3 8.5kb conversion tract: one tree, branch annotated with all 35 events."""
import os
from collections import defaultdict
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Tract definition ──────────────────────────────────────────────────────────
PALINDROME  = "P3"
TRACT_NODE  = "node3"
TRACT_POS   = [
    133711,133951,134589,134806,134893,135063,135084,135632,135851,136439,
    136804,136901,137192,137443,137781,138164,138332,138348,138715,138774,
    138822,138839,138999,139213,139363,139620,139698,139781,139849,140342,
    140398,140812,141351,141523,142198,
]
REP_POS     = str(TRACT_POS[len(TRACT_POS) // 2])   # middle site for genotypes

WORK_DIR    = "work"
OUTPUT_DIR  = "output"
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")

C_REF_INT  = "#00AAB5"
C_HET_INT  = "#FFD100"
C_ALT_INT  = "#F41C54"
C_REF_LEAF = "#80D5DA"
C_HET_LEAF = "#FFE880"
C_ALT_LEAF = "#F98EAA"
C_MISS     = "#D0D0D0"
LEAF_W, LEAF_H = 0.40, 0.28
INT_W,  INT_H  = 0.22, 0.15


# ── Tree helpers (same as 9_plot_trees.py) ────────────────────────────────────

def parse_tree(tree_file: str) -> tuple[str, dict, dict]:
    children: dict = {}
    parent: dict = {}
    with open(tree_file) as f:
        for line in f:
            node, c1, c2 = line.strip().split("\t")
            children[node] = [c1, c2]
            parent[c1] = parent[c2] = node
    root = next(n for n in children if n not in parent)
    return root, children, parent


def layout(root: str, children: dict) -> tuple[dict, list, int]:
    def dfs(n: str) -> list:
        return [n] if n not in children else dfs(children[n][0]) + dfs(children[n][1])
    leaves = dfs(root)
    depths: dict = {}
    def _d(n: str, d: int) -> None:
        depths[n] = d
        if n in children:
            for c in children[n]: _d(c, d + 1)
    _d(root, 0)
    max_d = max(depths.values())
    pos: dict = {lf: (i, max_d) for i, lf in enumerate(leaves)}
    def _assign_x(n: str) -> None:
        if n not in children: return
        l, r = children[n]
        _assign_x(l); _assign_x(r)
        pos[n] = ((pos[l][0] + pos[r][0]) / 2, depths[n])
    _assign_x(root)
    return pos, leaves, max_d


def _bottom_up(root: str, children: dict, leaf_genos: dict) -> dict:
    memo: dict = {}
    def _up(node: str) -> None:
        if node in memo: return
        if node not in children:
            memo[node] = [leaf_genos.get(node, "./.")]
            return
        l, r = children[node]
        _up(l); _up(r)
        g1, g2 = memo[l], memo[r]
        if len(g1) == 1 and len(g2) == 1:
            if g1 == g2: memo[node] = g1[:]
            else:
                combined = g1 + g2
                memo[node] = ["0/1"] if ("0/0" in combined and "1/1" in combined) else combined
        elif (len(g1) != 1) != (len(g2) != 1):
            ctr: dict = defaultdict(int)
            for g in g1 + g2: ctr[g] += 1
            memo[node] = ["0/1"] if len(ctr) == 3 else [max(ctr, key=ctr.get)]
        else:
            memo[node] = list(set(g1 + g2))
    _up(root)
    return memo


def _resolve(node: str, anc: str, children: dict, bu: dict, out: dict) -> None:
    g = bu[node]
    if len(g) != 1 and anc in g: resolved = anc
    elif len(g) != 1 and anc not in g: resolved = "0/1"
    else: resolved = g[0]
    out[node] = resolved
    if node in children:
        for c in children[node]: _resolve(c, resolved, children, bu, out)


def _count_muts(node: str, anc: str, children: dict, ng: dict) -> int:
    cur = ng[node]
    muts = 1 if (anc in ("0/0", "1/1") and cur == "0/1") else 0
    if node in children:
        for c in children[node]: muts += _count_muts(c, cur, children, ng)
    return muts


def reconstruct(root: str, children: dict, leaf_genos: dict) -> dict:
    bu = _bottom_up(root, children, leaf_genos)
    best: dict = {}
    best_m = 9999
    for root_state in bu[root]:
        ng: dict = {}
        _resolve(root, root_state, children, bu, ng)
        m = _count_muts(root, root_state, children, ng)
        if m < best_m: best_m, best = m, ng
    return best


def load_vcf_pos(vcf_file: str, pos: str) -> dict:
    genos: dict = {}
    with open(vcf_file) as f:
        samples = []
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.strip().split()[9:]
            if line.startswith("#"): continue
            parts = line.strip().split()
            if parts[1] != pos: continue
            for gt, samp in zip(parts[9:], samples):
                gt = gt.split(":")[0].replace("|", "/")
                if gt == "1/0": gt = "0/1"
                genos[samp] = gt
    return genos


def _node_color(geno: str, is_leaf: bool) -> str:
    if geno == "0/0": return C_REF_LEAF if is_leaf else C_REF_INT
    if geno == "1/1": return C_ALT_LEAF if is_leaf else C_ALT_INT
    if geno == "0/1": return C_HET_LEAF if is_leaf else C_HET_INT
    return C_MISS


# ── Main plot ─────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    vcf_file  = os.path.join(WORK_DIR, PALINDROME, "vcf",  f"{PALINDROME}.vcf")
    tree_file = os.path.join(WORK_DIR, PALINDROME, "tree", f"{PALINDROME}.tree.txt")

    root, children, parent = parse_tree(tree_file)
    pos_layout, leaves, _  = layout(root, children)
    n_leaves = len(leaves)
    max_y    = max(y for _, y in pos_layout.values())

    leaf_genos = load_vcf_pos(vcf_file, REP_POS)
    node_genos = reconstruct(root, children, leaf_genos)

    span_bp   = TRACT_POS[-1] - TRACT_POS[0]
    n_sites   = len(TRACT_POS)
    pos_range = f"{TRACT_POS[0]:,}–{TRACT_POS[-1]:,} bp"

    fig_w = max(14, n_leaves * 0.13)
    fig_h = 8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── Edges ─────────────────────────────────────────────────────────────────
    for node in children:
        nx, ny = pos_layout[node]
        l, r   = children[node]
        lx, ly = pos_layout[l]
        rx, ry = pos_layout[r]
        kw = dict(color="#BBBBBB", lw=0.4, zorder=1)
        ax.plot([lx, rx], [ny, ny], **kw)
        ax.plot([lx, lx], [ny, ly], **kw)
        ax.plot([rx, rx], [ny, ry], **kw)

    # ── Highlight node3 branch ────────────────────────────────────────────────
    if TRACT_NODE in pos_layout and TRACT_NODE in parent:
        nx, ny = pos_layout[TRACT_NODE]
        py     = pos_layout[parent[TRACT_NODE]][1]
        mid_y  = (ny + py) / 2
        # thicker highlighted branch segment
        ax.plot([nx, nx], [ny, py], color="#E96479", lw=2.0, zorder=2)
        ax.text(
            nx + 0.12, mid_y,
            f"GC tract\n{n_sites} sites\n{span_bp:,} bp\n({pos_range})",
            color="#E96479", fontsize=4.5, fontweight="bold",
            va="center", ha="left", zorder=5,
            bbox=dict(facecolor="white", edgecolor="#E96479",
                      linewidth=0.5, boxstyle="round,pad=0.3", alpha=0.9),
        )

    # ── Nodes ─────────────────────────────────────────────────────────────────
    for node, (nx, ny) in pos_layout.items():
        is_leaf = node not in children
        geno    = node_genos.get(node, "./.")
        color   = _node_color(geno, is_leaf)
        w, h    = (LEAF_W, LEAF_H) if is_leaf else (INT_W, INT_H)
        ax.add_patch(mpatches.Rectangle(
            (nx - w / 2, ny - h / 2), w, h,
            facecolor=color, edgecolor="white", linewidth=0.3, zorder=3,
        ))
        if is_leaf:
            ax.text(nx, ny + LEAF_H / 2 + 0.05, node,
                    rotation=90, ha="center", va="bottom",
                    fontsize=3, zorder=4,
                    bbox=dict(facecolor=color, edgecolor="none", pad=1.5, alpha=1.0))

    ax.set_xlim(-1, n_leaves)
    ax.set_ylim(-2.5, max_y + 2.5)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(
        f"P3 — node3 conversion tract   {n_sites} sites   {span_bp:,} bp   "
        f"(genotypes shown at pos {REP_POS})",
        fontsize=6,
    )
    ax.legend(
        handles=[
            mpatches.Patch(facecolor=C_REF_LEAF, label="0/0 (leaf)"),
            mpatches.Patch(facecolor=C_REF_INT,  label="0/0 (internal)"),
            mpatches.Patch(facecolor=C_HET_LEAF, label="0/1 (leaf)"),
            mpatches.Patch(facecolor=C_HET_INT,  label="0/1 (internal)"),
            mpatches.Patch(facecolor=C_ALT_LEAF, label="1/1 (leaf)"),
            mpatches.Patch(facecolor=C_ALT_INT,  label="1/1 (internal)"),
            mpatches.Patch(facecolor=C_MISS,     label="Missing"),
        ],
        loc="upper right", fontsize=4, framealpha=0.7, ncol=2, borderpad=0.5,
    )
    plt.tight_layout()

    out = os.path.join(OUTPUT_DIR, f"P3_node3_tract_{TIMESTAMP}.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
