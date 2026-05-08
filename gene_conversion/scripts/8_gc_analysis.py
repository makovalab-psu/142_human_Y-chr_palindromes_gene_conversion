#!/usr/bin/env python3
"""
Two analyses:
1. GC event frequency per position along each palindrome arm (hotspot map)
2. Estimated conversion tract length distribution per palindrome
"""
import os
from collections import defaultdict
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR = "work"
OUTPUT_DIR = "output"
MAX_DIST = 1000  # bp — max gap to consider sites part of the same tract
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Data loading ─────────────────────────────────────────────────────────────

# ── Analysis 1: GC frequency per position ────────────────────────────────────

def gc_frequency_per_position(events: list[dict]) -> dict[int, int]:
    """Count number of independent GC events per position."""
    counts: dict[int, int] = defaultdict(int)
    for e in events:
        if e["type"] == "geneconversion":
            counts[e["pos"]] += 1
    return counts


# ── Analysis 2: Tract lengths ─────────────────────────────────────────────────

def cluster_positions(positions: list[int], max_dist: int) -> list[list[int]]:
    if not positions:
        return []
    clusters = [[positions[0]]]
    for pos in positions[1:]:
        if pos - clusters[-1][-1] <= max_dist:
            clusters[-1].append(pos)
        else:
            clusters.append([pos])
    return clusters


def tract_lengths(events: list[dict]) -> list[int]:
    """
    For each node, find clusters of GC events within max_dist.
    Return the span (bp) of every cluster with >= 2 sites.
    """
    node_positions: dict[str, list[int]] = defaultdict(list)
    for e in events:
        if e["type"] == "geneconversion":
            node_positions[e["node"]].append(e["pos"])

    spans = []
    for positions in node_positions.values():
        positions = sorted(positions)
        for cluster in cluster_positions(positions, MAX_DIST):
            if len(cluster) >= 2:
                spans.append(cluster[-1] - cluster[0])
    return sorted(spans)


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_gc_frequency(palindrome: str, gc_counts: dict[int, int], arm_length: int) -> None:
    """Bar plot of GC event frequency along the arm."""
    if not gc_counts:
        return

    positions = sorted(gc_counts)
    counts = [gc_counts[p] for p in positions]

    fig, ax = plt.subplots(figsize=(12, 3))
    ax.bar(positions, counts, width=max(arm_length / 500, 50), color="#7DB9B6", linewidth=0)
    ax.set_xlim(0, arm_length)
    ax.set_xlabel("Position in arm (bp)")
    ax.set_ylabel("GC events")
    ax.set_title(f"{palindrome} — GC event frequency along arm")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x/1000)}kb"))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    plt.tight_layout()

    out = os.path.join(OUTPUT_DIR, f"{palindrome}_gc_frequency_{TIMESTAMP}.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_tract_lengths(palindrome: str, spans: list[int]) -> None:
    """Histogram of estimated conversion tract lengths."""
    if not spans:
        print(f"  {palindrome}: no tracts to plot")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, max(spans) + 100, min(40, len(set(spans)) + 1))
    ax.hist(spans, bins=bins, color="#E96479", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Tract span (bp)")
    ax.set_ylabel("Count")
    ax.set_title(f"{palindrome} — Estimated conversion tract lengths (n={len(spans)})")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x/1000)}kb" if x >= 1000 else str(int(x))))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    median = int(np.median(spans))
    ax.axvline(median, color="black", linestyle="--", linewidth=1, label=f"median={median}bp; n={len(spans)}")
    ax.legend(fontsize=9)
    plt.tight_layout()

    out = os.path.join(OUTPUT_DIR, f"{palindrome}_tract_lengths_{TIMESTAMP}.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved: {out}")


def _style_violin(parts: dict, color: str) -> None:
    for pc in parts["bodies"]:
        pc.set_facecolor(color)
        pc.set_alpha(0.75)
        pc.set_linewidth(0)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.5)
    for key in ("cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_visible(False)


def _despine(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def sample_count_from_vcf(palindrome: str) -> int:
    vcf = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    with open(vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return len(line.strip().split()) - 9
    return 0


def plot_combined_annotated(
    all_gc_counts: dict,
    all_arm_lengths: dict,
    all_events: dict,
) -> None:
    """
    3-column layout, (n palindromes + 1 summary) rows.
      col 1 — GC frequency bar plots in kb, shared x-axis
      col 2 — horizontal violin plots of tract lengths, shared x-axis
      col 3 — statistics in two lines (row1: N/rev/fix, row2: GCbias/tracts/med)
    Summary row: single text line spanning full width, no plots.
    """
    from matplotlib.gridspec import GridSpec

    plt.rcParams.update({"font.family": "Arial", "font.size": 5})

    n       = len(PALINDROMES)
    max_arm = max(all_arm_lengths[p] for p in PALINDROMES)
    kb_fmt  = mticker.FuncFormatter(lambda x, _: f"{x/1000:.0f}kb" if x >= 1000 else f"{int(x)}")
    bp_fmt  = mticker.FuncFormatter(lambda x, _: f"{int(x/1000)}kb" if x >= 1000 else str(int(x)))

    # Compute ruler ticks once so all subplots share the same positions
    _ruler_loc  = mticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
    ruler_ticks = [t for t in _ruler_loc.tick_values(0, max_arm) if 0 <= t <= max_arm]

    MM      = 1 / 25.4
    fig_w   = 180 * MM
    fig_h   = fig_w * (0.9 * (n + 1) / 18)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = GridSpec(
        n + 1, 3, figure=fig,
        width_ratios=[3.5, 2, 2],
        height_ratios=[1] * n + [0.15],
        hspace=0.40,
        wspace=0.10,
    )

    bar_axes  = [fig.add_subplot(gs[i, 0]) for i in range(n)]
    viol_axes = [fig.add_subplot(gs[i, 1]) for i in range(n)]
    ax_table  = fig.add_subplot(gs[0:n, 2])
    ax_ruler  = fig.add_subplot(gs[n, 0])

    # Share y within col 1; share x within col 2
    # Bar axes do NOT share x so each can have independent tick positions
    for ax in bar_axes[1:]:
        ax.sharey(bar_axes[0])
    for ax in viol_axes[1:]:
        ax.sharex(viol_axes[0])

    # Remove x-axis labels from all bar subplots (ruler provides the scale)
    for ax in bar_axes:
        ax.tick_params(labelbottom=False)
    for ax in viol_axes[:-1]:
        ax.tick_params(labelbottom=False)

    cum_spans: list[int] = []
    cum_gc = cum_anc = cum_der = cum_at_gc = cum_gc_at = 0
    table_rows: list[list[str]] = []

    for i, palindrome in enumerate(PALINDROMES):
        ax_bar  = bar_axes[i]
        ax_viol = viol_axes[i]

        gc_counts  = all_gc_counts[palindrome]
        arm_length = all_arm_lengths[palindrome]
        events     = all_events[palindrome]

        # ── col 1: frequency bar plot ─────────────────────────────────────────
        if gc_counts:
            positions = sorted(gc_counts)
            counts    = [gc_counts[p] for p in positions]
            ax_bar.bar(positions, counts, width=max(max_arm / 500, 50),
                       color="#7DB9B6", linewidth=0, align="center")
        ax_bar.axvline(arm_length, color="#888888", linewidth=0.5, linestyle="--", zorder=2)
        ax_bar.set_xlim(0, max_arm)
        n_samples = sample_count_from_vcf(palindrome)
        ax_bar.set_ylabel(f"{palindrome}\n(N={n_samples})", rotation=0, labelpad=22, va="center", fontsize=5)
        ax_bar.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=3))
        ax_bar.tick_params(labelsize=5)
        ax_bar.spines["bottom"].set_bounds(0, arm_length)
        ax_bar.set_xticks([t for t in ruler_ticks if t <= arm_length])
        _despine(ax_bar)

        # ── col 2: horizontal violin ──────────────────────────────────────────
        spans = tract_lengths(events)
        cum_spans.extend(spans)
        ax_viol.set_yticks([])
        ax_viol.set_ylim(0.4, 1.6)
        if len(spans) >= 3:
            _style_violin(
                ax_viol.violinplot([spans], positions=[1], vert=False, showmedians=True),
                "#E96479",
            )
        elif spans:
            ax_viol.scatter(spans, [1] * len(spans), color="#E96479", s=4, alpha=0.6, zorder=3)
        if spans:
            med_v = int(np.median(spans))
            ax_viol.text(med_v, 1.25, f"{med_v}bp n={len(spans)}",
                         ha="left", va="bottom", fontsize=4.5, color="#333333")
        ax_viol.tick_params(axis="x", labelsize=5)
        _despine(ax_viol)

        # ── col 3: statistics — two lines ────────────────────────────────────
        n_gc               = sum(1 for e in events if e["type"] == "geneconversion")
        to_anc, to_der     = gc_resolution_fate(events)
        at_to_gc, gc_to_at = gc_bias(events)
        pct_bias           = 100 * at_to_gc / (at_to_gc + gc_to_at) if (at_to_gc + gc_to_at) else 0
        n_tracts           = len(spans)
        med_span           = int(np.median(spans)) if spans else 0

        cum_gc    += n_gc
        cum_anc   += to_anc
        cum_der   += to_der
        cum_at_gc += at_to_gc
        cum_gc_at += gc_to_at

        ev_per_kb = n_gc / (arm_length / 1000)
        table_rows.append([
            str(n_gc),
            f"{ev_per_kb:.2f}",
            str(to_anc),
            str(to_der),
            f"{pct_bias:.0f}",
            # str(n_tracts),
        ])

    # ── ruler axis below all bar plots ───────────────────────────────────────
    ax_ruler.set_xlim(0, max_arm)
    ax_ruler.set_ylim(0, 1)
    ax_ruler.set_yticks([])
    ax_ruler.set_xticks(ruler_ticks)
    ax_ruler.xaxis.set_major_formatter(kb_fmt)
    ax_ruler.tick_params(axis="x", labelsize=5, top=False)
    ax_ruler.spines["top"].set_visible(False)
    ax_ruler.spines["left"].set_visible(False)
    ax_ruler.spines["right"].set_visible(False)

    # viol_axes[-1].set_xlabel("Tract span")
    # viol_axes[-1].xaxis.set_major_formatter(bp_fmt)

    # ── table: totals row + render ────────────────────────────────────────────
    pct_bias_all  = 100 * cum_at_gc / (cum_at_gc + cum_gc_at) if (cum_at_gc + cum_gc_at) else 0
    total_arm_kb  = sum(all_arm_lengths[p] for p in PALINDROMES) / 1000
    ev_per_kb_all = cum_gc / total_arm_kb
    table_rows.append([
        str(cum_gc),
        f"{ev_per_kb_all:.2f}",
        str(cum_anc),
        str(cum_der),
        f"{pct_bias_all:.0f}",
        # str(len(cum_spans)),
    ])

    col_labels = ["events", "events/kb", "reversions", "fixations", "GC bias (%)"]
    row_labels = PALINDROMES + ["Sum"]
    n_rows = len(row_labels)
    n_cols = len(col_labels)

    ax_table.axis("off")
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(0, 1)

    import matplotlib.patches as mpatches

    row_lbl_w = 0.18
    col_w     = (1 - row_lbl_w) / n_cols
    hdr_h     = 0.22
    row_h     = (1 - hdr_h) / n_rows
    fs        = 5

    tr = ax_table.transAxes

    # top-left empty header cell
    ax_table.add_patch(mpatches.FancyBboxPatch(
        (0, 1 - hdr_h), row_lbl_w, hdr_h,
        boxstyle="square,pad=0", facecolor="#DDDDDD", edgecolor="#AAAAAA", lw=0.4,
        transform=tr, clip_on=False))

    # column header cells — rotated text anchored near cell bottom
    for j, lbl in enumerate(col_labels):
        x = row_lbl_w + j * col_w
        y = 1 - hdr_h
        ax_table.add_patch(mpatches.FancyBboxPatch(
            (x, y), col_w, hdr_h,
            boxstyle="square,pad=0", facecolor="#DDDDDD", edgecolor="#AAAAAA", lw=0.4,
            transform=tr, clip_on=False))
        ax_table.text(x + col_w / 2, y + 0.03, lbl,
                      rotation=90, va="bottom", ha="center",
                      fontsize=fs, fontweight="bold",
                      transform=tr, clip_on=False)

    # data rows
    for i, (row_lbl, row_data) in enumerate(zip(row_labels, table_rows)):
        y   = 1 - hdr_h - (i + 1) * row_h
        is_total = (i == n_rows - 1)
        bg  = "#F0F0F0" if is_total else ("white" if i % 2 == 0 else "#FAFAFA")
        fw  = "bold" if is_total else "normal"

        # row label cell
        ax_table.add_patch(mpatches.FancyBboxPatch(
            (0, y), row_lbl_w, row_h,
            boxstyle="square,pad=0", facecolor=bg, edgecolor="#AAAAAA", lw=0.4,
            transform=tr, clip_on=False))
        ax_table.text(row_lbl_w - 0.02, y + row_h / 2, row_lbl,
                      va="center", ha="right", fontsize=fs, fontweight="bold",
                      transform=tr, clip_on=False)

        # data cells
        for j, val in enumerate(row_data):
            x = row_lbl_w + j * col_w
            ax_table.add_patch(mpatches.FancyBboxPatch(
                (x, y), col_w, row_h,
                boxstyle="square,pad=0", facecolor=bg, edgecolor="#AAAAAA", lw=0.4,
                transform=tr, clip_on=False))
            ax_table.text(x + col_w - 0.02, y + row_h / 2, val,
                          va="center", ha="right", fontsize=fs, fontweight=fw,
                          transform=tr, clip_on=False)

    ax_table.set_title("Summary statistics", pad=3)

    # ── column headers ────────────────────────────────────────────────────────
    bar_axes[0].set_title("Gene conversion event frequency", pad=3)
    viol_axes[0].set_title("Tract length distribution", pad=3)

    fig.suptitle("Gene conversion events along palindrome arms")
    plt.tight_layout(rect=[0, 0.02, 1.0, 0.97], h_pad=0.3)

    out = os.path.join(OUTPUT_DIR, f"all_palindromes_gc_frequency_{TIMESTAMP}.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved combined: {out}")


# ── Analysis 3: GC resolution fate at mutated sites ──────────────────────────

def gc_resolution_fate(events: list[dict]) -> tuple[int, int]:
    """
    At sites where a mutation created heterozygosity (root = ancestral 0/0),
    count how GC events resolved the het:
      - back to ancestral (0/1 → 0/0): ancestral field == "derived->ancestral"
      - forward to derived  (0/1 → 1/1): ancestral field == "ancestral->derived"
    Excludes ancestral_heterozygote events (root was already 0/1).
    """
    to_ancestral = sum(1 for e in events if e["type"] == "geneconversion" and e["ancestral"] == "derived->ancestral")
    to_derived   = sum(1 for e in events if e["type"] == "geneconversion" and e["ancestral"] == "ancestral->derived")
    return to_ancestral, to_derived


# ── Analysis 4: GC bias in gene conversion events ────────────────────────────

def gc_bias(events: list[dict]) -> tuple[int, int]:
    """
    Among gene conversion events, count AT→GC vs GC→AT conversions.
    GC↔GC (G↔C) and AT↔AT (A↔T) transversions are excluded.
    Returns (at_to_gc, gc_to_at).
    """
    AT = {"A", "T"}
    GC = {"G", "C"}
    at_to_gc = 0
    gc_to_at = 0
    for e in events:
        if e["type"] != "geneconversion":
            continue
        from_allele, to_allele = e["base_change"].split("->")
        if from_allele in AT and to_allele in GC:
            at_to_gc += 1
        elif from_allele in GC and to_allele in AT:
            gc_to_at += 1
    return at_to_gc, gc_to_at


# ── Data loading (extended) ───────────────────────────────────────────────────

def load_events_full(palindrome: str) -> tuple[list[dict], int]:
    """Like load_events but includes all fields needed for analyses 3 and 4."""
    events_file = os.path.join(WORK_DIR, palindrome, "events", f"{palindrome}.events.txt")
    events = []
    with open(events_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            events.append({
                "pos": int(fields[1]),
                "from": fields[2],
                "type": fields[6],
                "node": fields[7],
                "base_change": fields[4],
                "ancestral": fields[5],
            })

    vcf_file = os.path.join(WORK_DIR, palindrome, "vcf", f"{palindrome}.vcf")
    max_pos = 0
    with open(vcf_file) as f:
        for line in f:
            if not line.startswith("#"):
                max_pos = max(max_pos, int(line.split("\t")[1]))

    return events, max_pos


# ── Analysis 5: Event and site counts by type ────────────────────────────────

def event_site_counts(events: list[dict]) -> dict:
    """
    Returns counts of events and sites for gene conversion and mutation types.
    'sites' = unique positions with at least one event of that type.
    """
    gc_sites  = {e["pos"] for e in events if e["type"] == "geneconversion"}
    mut_sites = {e["pos"] for e in events if e["type"] == "mutation"}
    return {
        "n_gc_events":  sum(1 for e in events if e["type"] == "geneconversion"),
        "n_mut_events": sum(1 for e in events if e["type"] == "mutation"),
        "n_gc_sites":   len(gc_sites),
        "n_mut_sites":  len(mut_sites),
        "n_both_sites": len(gc_sites & mut_sites),
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    all_gc_counts = {}
    all_arm_lengths = {}
    all_events = {}

    print("=== Analysis 1: GC frequency per position ===")
    for palindrome in PALINDROMES:
        events, arm_length = load_events_full(palindrome)
        all_events[palindrome] = events
        all_arm_lengths[palindrome] = arm_length

        gc_counts = gc_frequency_per_position(events)
        all_gc_counts[palindrome] = gc_counts

        n_sites = len(gc_counts)
        max_count = max(gc_counts.values()) if gc_counts else 0
        print(f"  {palindrome}: {n_sites} sites with GC events, max={max_count} independent events at one site")
        # plot_gc_frequency(palindrome, gc_counts, arm_length)

    plot_combined_annotated(all_gc_counts, all_arm_lengths, all_events)

    print("\n=== Analysis 2: Tract length distribution ===")
    for palindrome in PALINDROMES:
        events = all_events[palindrome]
        spans = tract_lengths(events)
        if spans:
            print(f"  {palindrome}: {len(spans)} tracts, median={int(np.median(spans))}bp, "
                  f"max={max(spans)}bp, mean={int(np.mean(spans))}bp")
        else:
            print(f"  {palindrome}: no multi-site tracts within {MAX_DIST}bp")
        # plot_tract_lengths(palindrome, spans)

    print("\n=== Analysis 3: GC resolution fate at mutated sites ===")
    print(f"  {'palindrome':<12} {'→ancestral':>12} {'→derived':>10} {'total':>8} {'%→ancestral':>13}")
    print("  " + "-" * 58)
    for palindrome in PALINDROMES:
        to_anc, to_der = gc_resolution_fate(all_events[palindrome])
        total = to_anc + to_der
        pct = 100 * to_anc / total if total else 0
        print(f"  {palindrome:<12} {to_anc:>12} {to_der:>10} {total:>8} {pct:>12.1f}%")

    print("\n=== Analysis 4: GC bias in gene conversion events ===")
    print(f"  {'palindrome':<12} {'AT→GC':>8} {'GC→AT':>8} {'total_GC':>10} {'%AT→GC':>8}")
    print("  " + "-" * 50)
    for palindrome in PALINDROMES:
        at_to_gc, gc_to_at = gc_bias(all_events[palindrome])
        total = at_to_gc + gc_to_at
        pct = 100 * at_to_gc / total if total else 0
        print(f"  {palindrome:<12} {at_to_gc:>8} {gc_to_at:>8} {total:>10} {pct:>7.1f}%")

    print("\n=== Analysis 5: Event and site counts by type ===")
    print(f"  {'palindrome':<12} {'GC events':>10} {'mut events':>11} {'GC sites':>9} {'mut sites':>10} {'both sites':>11}")
    print("  " + "-" * 65)
    totals = {"n_gc_events": 0, "n_mut_events": 0, "n_gc_sites": 0, "n_mut_sites": 0, "n_both_sites": 0}
    for palindrome in PALINDROMES:
        c = event_site_counts(all_events[palindrome])
        for k in totals:
            totals[k] += c[k]
        print(f"  {palindrome:<12} {c['n_gc_events']:>10} {c['n_mut_events']:>11} {c['n_gc_sites']:>9} {c['n_mut_sites']:>10} {c['n_both_sites']:>11}")
    print("  " + "-" * 65)
    print(f"  {'Sum':<12} {totals['n_gc_events']:>10} {totals['n_mut_events']:>11} {totals['n_gc_sites']:>9} {totals['n_mut_sites']:>10} {totals['n_both_sites']:>11}")

    # ── Write summary TSV ─────────────────────────────────────────────────────
    tsv_path = os.path.join(OUTPUT_DIR, f"summary_statistics_{TIMESTAMP}.tsv")
    cols = [
        "palindrome", "n_samples", "arm_length_bp",
        "gc_events", "gc_events_per_kb", "reversions", "fixations", "gc_bias_pct",
        "n_tracts", "median_tract_bp", "mean_tract_bp", "max_tract_bp",
        "gc_sites", "mut_events", "mut_sites", "sites_gc_and_mut",
        "at_to_gc", "gc_to_at",
    ]
    rows = []
    for palindrome in PALINDROMES:
        events     = all_events[palindrome]
        arm_length = all_arm_lengths[palindrome]
        n_samples  = sample_count_from_vcf(palindrome)
        n_gc       = sum(1 for e in events if e["type"] == "geneconversion")
        to_anc, to_der = gc_resolution_fate(events)
        at_gc, gc_at   = gc_bias(events)
        bias_pct   = 100 * at_gc / (at_gc + gc_at) if (at_gc + gc_at) else 0
        spans      = tract_lengths(events)
        c          = event_site_counts(events)
        rows.append({
            "palindrome":       palindrome,
            "n_samples":        n_samples,
            "arm_length_bp":    arm_length,
            "gc_events":        n_gc,
            "gc_events_per_kb": round(n_gc / (arm_length / 1000), 4),
            "reversions":       to_anc,
            "fixations":        to_der,
            "gc_bias_pct":      round(bias_pct, 2),
            "n_tracts":         len(spans),
            "median_tract_bp":  int(np.median(spans)) if spans else 0,
            "mean_tract_bp":    int(np.mean(spans))   if spans else 0,
            "max_tract_bp":     max(spans)             if spans else 0,
            "gc_sites":         c["n_gc_sites"],
            "mut_events":       c["n_mut_events"],
            "mut_sites":        c["n_mut_sites"],
            "sites_gc_and_mut": c["n_both_sites"],
            "at_to_gc":         at_gc,
            "gc_to_at":         gc_at,
        })

    with open(tsv_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"\nSummary TSV: {tsv_path}")


if __name__ == "__main__":
    main()
