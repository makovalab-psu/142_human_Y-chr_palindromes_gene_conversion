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
    ax.axvline(median, color="black", linestyle="--", linewidth=1, label=f"median={median}bp")
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
    fig_h   = fig_w * (1.5 * (n + 1) / 18)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = GridSpec(
        n + 1, 3, figure=fig,
        width_ratios=[4, 2, 1.5],
        height_ratios=[1] * n + [0.15],
        hspace=0.30,
    )

    bar_axes  = [fig.add_subplot(gs[i, 0]) for i in range(n)]
    viol_axes = [fig.add_subplot(gs[i, 1]) for i in range(n)]
    stat_axes = [fig.add_subplot(gs[i, 2]) for i in range(n)]
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

    for i, palindrome in enumerate(PALINDROMES):
        ax_bar  = bar_axes[i]
        ax_viol = viol_axes[i]
        ax_stat = stat_axes[i]

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
            ax_viol.text(med_v, 1.25, f"{med_v}bp",
                         ha="center", va="bottom", fontsize=4.5, color="#333333")
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

        ax_stat.axis("off")
        kw = dict(transform=ax_stat.transAxes, va="center", ha="left", family="monospace", fontsize=6)
        ax_stat.text(0.05, 0.78, f"events={n_gc}   events/kb={ev_per_kb:.2f}", **kw)
        ax_stat.text(0.05, 0.50, f"reversions={to_anc}   fixations={to_der}", **kw)
        ax_stat.text(0.05, 0.22, f"GC_bias={pct_bias:.0f}% ({at_to_gc}/{gc_to_at})   tracts={n_tracts}", **kw)

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

    viol_axes[-1].set_xlabel("Tract span")
    viol_axes[-1].xaxis.set_major_formatter(bp_fmt)

    # ── summary line: fixed position at figure bottom ────────────────────────
    pct_bias_all  = 100 * cum_at_gc / (cum_at_gc + cum_gc_at) if (cum_at_gc + cum_gc_at) else 0
    med_all       = int(np.median(cum_spans)) if cum_spans else 0
    max_all       = max(cum_spans) if cum_spans else 0
    total_arm_kb  = sum(all_arm_lengths[p] for p in PALINDROMES) / 1000
    ev_per_kb_all = cum_gc / total_arm_kb
    fig.text(0.5, 0.035,
             f"ALL — events={cum_gc}   events/kb={ev_per_kb_all:.2f}   reversions={cum_anc}   fixations={cum_der}",
             va="bottom", ha="center", family="monospace", fontsize=6)
    fig.text(0.5, 0.01,
             f"ALL — GC_bias={pct_bias_all:.0f}% ({cum_at_gc}/{cum_gc_at})   tracts={len(cum_spans)}",
             va="bottom", ha="center", family="monospace", fontsize=6)

    # ── column headers ────────────────────────────────────────────────────────
    bar_axes[0].set_title("GC event frequency", pad=3)
    viol_axes[0].set_title("Tract length distribution", pad=3)
    stat_axes[0].set_title("Summary statistics", pad=3)

    fig.suptitle("GC events along palindrome arms")
    plt.tight_layout(rect=[0, 0.07, 1.0, 0.97], h_pad=0.3)

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
        plot_gc_frequency(palindrome, gc_counts, arm_length)

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
        plot_tract_lengths(palindrome, spans)

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


if __name__ == "__main__":
    main()
