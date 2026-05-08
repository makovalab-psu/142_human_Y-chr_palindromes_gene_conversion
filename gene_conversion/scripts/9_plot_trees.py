#!/usr/bin/env python3
"""Plot per-site parsimony trees with colored nodes and MUT/GC event labels."""
import os
from collections import defaultdict
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR    = "work"
OUTPUT_DIR  = "output"
TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")

# Colors  (internal = darker, leaf = lighter)
C_REF_INT  = "#00AAB5"   # 0/0 internal – RGB(0,170,181)
C_HET_INT  = "#FFD100"   # 0/1 internal – RGB(255,209,0)
C_ALT_INT  = "#F41C54"   # 1/1 internal – RGB(244,28,84)
C_REF_LEAF = "#80D5DA"   # 0/0 leaf (lightened)
C_HET_LEAF = "#FFE880"   # 0/1 leaf (lightened)
C_ALT_LEAF = "#F98EAA"   # 1/1 leaf (lightened)
C_MISS     = "#D0D0D0"   # grey – missing / unknown
C_MUT      = "#222222"   # MUT label
C_GC       = "#222222"   # GC label

# Node rectangle sizes (in leaf-spacing data units)
LEAF_W, LEAF_H = 0.40, 0.28
INT_W,  INT_H  = 0.22, 0.15


# ── Parsimony ─────────────────────────────────────────────────────────────────

def _bottom_up(root: str, children: dict, leaf_genos: dict) -> dict:
    """Return {node: [possible genotypes]} via Fitch parsimony bottom-up pass."""
    memo: dict = {}

    def _up(node: str) -> None:
        if node in memo:
            return
        if node not in children:
            memo[node] = [leaf_genos.get(node, "./.")]
            return
        l, r = children[node]
        _up(l); _up(r)
        g1, g2 = memo[l], memo[r]
        if len(g1) == 1 and len(g2) == 1:
            if g1 == g2:
                memo[node] = g1[:]
            else:
                combined = g1 + g2
                memo[node] = ["0/1"] if ("0/0" in combined and "1/1" in combined) else combined
        elif (len(g1) != 1) != (len(g2) != 1):  # exactly one ambiguous
            ctr: dict = defaultdict(int)
            for g in g1 + g2:
                ctr[g] += 1
            memo[node] = ["0/1"] if len(ctr) == 3 else [max(ctr, key=ctr.get)]
        else:
            memo[node] = list(set(g1 + g2))

    _up(root)
    return memo


def _resolve(node: str, anc: str, children: dict, bu: dict, out: dict) -> None:
    """Top-down pass: resolve ambiguous nodes, record resolved state in out."""
    g = bu[node]
    if len(g) != 1 and anc in g:
        resolved = anc
    elif len(g) != 1 and anc not in g:
        resolved = "0/1"
    else:
        resolved = g[0]
    out[node] = resolved
    if node in children:
        for c in children[node]:
            _resolve(c, resolved, children, bu, out)


def _count_muts(node: str, anc: str, children: dict, ng: dict) -> int:
    cur = ng[node]
    muts = 1 if (anc in ("0/0", "1/1") and cur == "0/1") else 0
    if node in children:
        for c in children[node]:
            muts += _count_muts(c, cur, children, ng)
    return muts


def reconstruct(root: str, children: dict, leaf_genos: dict) -> dict:
    """Full parsimony reconstruction; returns {node: resolved_genotype_str}."""
    bu = _bottom_up(root, children, leaf_genos)
    best: dict = {}
    best_m = 9999
    for root_state in bu[root]:
        ng: dict = {}
        _resolve(root, root_state, children, bu, ng)
        m = _count_muts(root, root_state, children, ng)
        if m < best_m:
            best_m, best = m, ng
    return best


# ── Tree layout ────────────────────────────────────────────────────────────────

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
    """
    Cladogram layout: root at y=0, all leaves at y=max_depth.
    Returns {node: (x, y)}, ordered leaves, max_depth.
    """
    def dfs(n: str) -> list:
        return [n] if n not in children else dfs(children[n][0]) + dfs(children[n][1])

    leaves = dfs(root)

    depths: dict = {}
    def _d(n: str, d: int) -> None:
        depths[n] = d
        if n in children:
            for c in children[n]:
                _d(c, d + 1)
    _d(root, 0)
    max_d = max(depths.values())

    pos: dict = {lf: (i, max_d) for i, lf in enumerate(leaves)}

    def _assign_x(n: str) -> None:
        if n not in children:
            return
        l, r = children[n]
        _assign_x(l); _assign_x(r)
        pos[n] = ((pos[l][0] + pos[r][0]) / 2, depths[n])
    _assign_x(root)

    return pos, leaves, max_d


# ── Plotting ──────────────────────────────────────────────────────────────────

def _node_color(geno: str, is_leaf: bool) -> str:
    if geno == "0/0":
        return C_REF_LEAF if is_leaf else C_REF_INT
    if geno == "1/1":
        return C_ALT_LEAF if is_leaf else C_ALT_INT
    if geno == "0/1":
        return C_HET_LEAF if is_leaf else C_HET_INT
    return C_MISS


def plot_site_tree(
    palindrome: str,
    site_pos: str,
    events: list[dict],
    node_genos: dict,
    pos: dict,
    children: dict,
    parent: dict,
    leaves: list,
) -> plt.Figure:
    n_leaves = len(leaves)
    max_y = max(y for _, y in pos.values())

    fig_w = max(12, n_leaves * 0.13)
    fig_h = 8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── Edges ─────────────────────────────────────────────────────────────────
    for node in children:
        nx, ny = pos[node]
        l, r   = children[node]
        lx, ly = pos[l]
        rx, ry = pos[r]
        kw = dict(color="#BBBBBB", lw=0.4, zorder=1)
        ax.plot([lx, rx], [ny, ny], **kw)
        ax.plot([lx, lx], [ny, ly], **kw)
        ax.plot([rx, rx], [ny, ry], **kw)

    # ── Event labels on branches ───────────────────────────────────────────────
    event_map: dict = {e["node"]: e["type"] for e in events}
    for node_name, evt_type in event_map.items():
        if node_name not in pos:
            continue
        nx, ny = pos[node_name]
        py = pos[parent[node_name]][1] if node_name in parent else ny - 1
        mid_y  = (ny + py) / 2
        label  = "MUT" if evt_type == "mutation" else "GC"
        color  = C_MUT if evt_type == "mutation" else C_GC
        ax.text(nx + 0.08, mid_y, label,
                color=color, fontsize=4.5, fontweight="bold",
                va="center", ha="left", zorder=5)

    # ── Nodes ─────────────────────────────────────────────────────────────────
    for node, (nx, ny) in pos.items():
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

    # ── Axes / labels ──────────────────────────────────────────────────────────
    ax.set_xlim(-1, n_leaves)
    ax.set_ylim(-2.5, max_y + 2.5)
    ax.invert_yaxis()
    ax.axis("off")

    event_types = sorted({e["type"] for e in events})
    ax.set_title(
        f"{palindrome}   pos={site_pos}   [{', '.join(event_types)}]",
        fontsize=6,
    )

    # ── Legend ────────────────────────────────────────────────────────────────
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
    return fig


# ── Data loading ──────────────────────────────────────────────────────────────

def load_vcf(vcf_file: str) -> tuple[list, dict]:
    samples: list = []
    genos: dict   = defaultdict(dict)
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.strip().split()[9:]
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            pos   = parts[1]
            for gt, samp in zip(parts[9:], samples):
                gt = gt.split(":")[0].replace("|", "/")
                if gt == "1/0":
                    gt = "0/1"
                genos[pos][samp] = gt
    return samples, genos


def load_events(events_file: str) -> dict:
    by_pos: dict = defaultdict(list)
    with open(events_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            by_pos[fields[1]].append({"type": fields[6], "node": fields[7]})
    return by_pos


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for palindrome in PALINDROMES:
        vcf_file    = os.path.join(WORK_DIR, palindrome, "vcf",    f"{palindrome}.vcf")
        tree_file   = os.path.join(WORK_DIR, palindrome, "tree",   f"{palindrome}.tree.txt")
        events_file = os.path.join(WORK_DIR, palindrome, "events", f"{palindrome}.events.txt")
        out_pdf     = os.path.join(OUTPUT_DIR, f"{palindrome}_trees_{TIMESTAMP}.pdf")

        root, children, parent = parse_tree(tree_file)
        pos_layout, leaves, _  = layout(root, children)

        _, genos_by_pos  = load_vcf(vcf_file)
        events_by_pos    = load_events(events_file)
        positions        = sorted(events_by_pos, key=int)

        print(f"{palindrome}: {len(positions)} positions → {out_pdf}")

        with PdfPages(out_pdf) as pdf:
            for i, site_pos in enumerate(positions, 1):
                leaf_genos = genos_by_pos.get(site_pos, {})
                node_genos = reconstruct(root, children, leaf_genos)
                events     = events_by_pos[site_pos]

                fig = plot_site_tree(
                    palindrome, site_pos, events,
                    node_genos, pos_layout, children, parent, leaves,
                )
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

                if i % 50 == 0:
                    print(f"  {i}/{len(positions)}")

        print(f"  saved {out_pdf}")


if __name__ == "__main__":
    main()
