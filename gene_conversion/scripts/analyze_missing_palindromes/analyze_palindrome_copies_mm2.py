#!/usr/bin/env python3
"""
Palindrome copy number & indel analysis via minimap2.

Batch mode (default):
    python analyze_palindrome_copies_mm2.py \\
        --query QUERY.fa \\
        --assembly-dir /path/to/assemblies \\
        --assembly-ext fa \\
        --outdir results/ \\
        --large-indel-threshold 50 \\
        [--preset asm5] \\
        [--min-mapq 0]

Single-sample detailed mode:
    python analyze_palindrome_copies_mm2.py \\
        --query QUERY.fa \\
        --assembly-dir /path/to/assemblies \\
        --outdir results/ \\
        --sample HG01890 \\
        [--detail-qcov-min 20] \\
        [--detail-indel-min 1]
"""
import argparse
import re
import subprocess
from pathlib import Path


QCOV_FILTER  = 80.0   # batch mode: minimum qcov to keep a hit
QCOV_FLAG    = 90.0   # batch mode: flag if any hit below this
OVERLAP_FRAC = 0.50   # flag same-contig hits overlapping > this fraction of qlen


# ── cs tag parsing ────────────────────────────────────────────────────────────

def parse_cs_indels(cs: str, threshold: int) -> list[dict]:
    """
    Parse minimap2 cs tag for indels >= threshold bp.
      +seq — insertion in query (gap in reference)
      -seq — deletion from query (gap in query)
    """
    indels = []
    for op, seq in re.findall(r"([:\*\+\-])([A-Za-z0-9]+)", cs):
        if op == "+" and len(seq) >= threshold:
            indels.append({"type": "insertion", "size": len(seq)})
        elif op == "-" and len(seq) >= threshold:
            indels.append({"type": "deletion", "size": len(seq)})
    return indels


def cs_substitution_count(cs: str) -> int:
    """Count substitutions (*xy) in a cs tag."""
    return len(re.findall(r"\*[a-z]{2}", cs))


# ── minimap2 ──────────────────────────────────────────────────────────────────

def run_minimap2_raw(
    query: Path,
    target: Path,
    preset: str,
    n_secondary: int,
    extra_flags: list[str] | None = None,
) -> str:
    """Return raw PAF output as a string."""
    cmd = [
        "minimap2",
        "-x", preset,
        "--cs",
        "--secondary=yes",
        "-N", str(n_secondary),
        *(extra_flags or []),
        str(target),
        str(query),
    ]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def parse_paf(paf: str, qcov_min: float, min_mapq: int) -> list[dict]:
    """Parse PAF text into list of hit dicts, applying qcov and mapq filters."""
    hits = []
    for line in paf.splitlines():
        if not line:
            continue
        f = line.split("\t")
        q_len   = int(f[1])
        q_start = int(f[2])
        q_end   = int(f[3])
        strand  = "plus" if f[4] == "+" else "minus"
        t_name  = f[5]
        t_start = int(f[7])
        t_end   = int(f[8])
        matches = int(f[9])
        aln_len = int(f[10])
        mapq    = int(f[11])

        if mapq < min_mapq:
            continue

        qcov   = (q_end - q_start) / q_len * 100
        pident = matches / aln_len * 100 if aln_len > 0 else 0.0

        if qcov < qcov_min:
            continue

        cs_tag = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")

        hits.append({
            "q_len":   q_len,
            "q_start": q_start,
            "q_end":   q_end,
            "t_name":  t_name,
            "t_start": t_start,
            "t_end":   t_end,
            "strand":  strand,
            "pident":  pident,
            "qcov":    qcov,
            "mapq":    mapq,
            "cs":      cs_tag,
            "indels":  parse_cs_indels(cs_tag, threshold=0),
        })
    return hits


# ── Flagging ──────────────────────────────────────────────────────────────────

def flag_assembly(hits: list[dict]) -> list[str]:
    reasons = []
    if not hits:
        return ["no_hits"]

    qlen = hits[0]["q_len"]

    if len(hits) != 2:
        reasons.append(f"hit_count={len(hits)}")

    for h in hits:
        if h["qcov"] < QCOV_FLAG:
            reasons.append(f"low_qcov={h['qcov']:.1f}%_on_{h['t_name']}")

    if len(hits) == 2:
        if hits[0]["strand"] == hits[1]["strand"]:
            reasons.append(f"same_strand={hits[0]['strand']}")

        if hits[0]["t_name"] == hits[1]["t_name"]:
            s0 = (min(hits[0]["t_start"], hits[0]["t_end"]),
                  max(hits[0]["t_start"], hits[0]["t_end"]))
            s1 = (min(hits[1]["t_start"], hits[1]["t_end"]),
                  max(hits[1]["t_start"], hits[1]["t_end"]))
            overlap = max(0, min(s0[1], s1[1]) - max(s0[0], s1[0]))
            if overlap > OVERLAP_FRAC * qlen:
                reasons.append(f"overlapping_hits_on_{hits[0]['t_name']}")

    return reasons


# ── Batch mode ────────────────────────────────────────────────────────────────

def run_batch(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    assemblies = sorted(args.assembly_dir.glob(f"*.{args.assembly_ext}"))
    if not assemblies:
        raise SystemExit(f"No *.{args.assembly_ext} files in {args.assembly_dir}")

    print(f"Query:      {args.query}")
    print(f"Assemblies: {len(assemblies)}")
    print(f"Preset:     {args.preset}")
    print(f"Outdir:     {args.outdir}")
    print(f"Indel threshold: {args.large_indel_threshold} bp\n")

    summary_rows: list[dict] = []
    count_rows:   list[dict] = []
    flagged:      dict[str, list[str]] = {}

    for assembly_fa in assemblies:
        name = assembly_fa.stem
        paf  = run_minimap2_raw(args.query, assembly_fa, args.preset, n_secondary=10)
        hits = parse_paf(paf, qcov_min=QCOV_FILTER, min_mapq=args.min_mapq)

        if not hits:
            summary_rows.append({
                "assembly": name, "hit_num": 0,
                "subject_seq": ".", "t_start": ".", "t_end": ".", "strand": ".",
                "pident": ".", "qcov": ".", "mapq": ".",
                "n_large_indels": 0, "indel_details": ".",
            })
        else:
            for k, h in enumerate(hits, 1):
                large = [d for d in h["indels"] if d["size"] >= args.large_indel_threshold]
                summary_rows.append({
                    "assembly":       name,
                    "hit_num":        k,
                    "subject_seq":    h["t_name"],
                    "t_start":        h["t_start"],
                    "t_end":          h["t_end"],
                    "strand":         h["strand"],
                    "pident":         f"{h['pident']:.2f}",
                    "qcov":           f"{h['qcov']:.2f}",
                    "mapq":           h["mapq"],
                    "n_large_indels": len(large),
                    "indel_details":  ";".join(f"{d['type']}:{d['size']}" for d in large) or ".",
                })

        n_with_indels = sum(
            1 for h in hits
            if any(d["size"] >= args.large_indel_threshold for d in h["indels"])
        )
        count_rows.append({
            "assembly":               name,
            "total_hits":             len(hits),
            "hits_with_large_indels": n_with_indels,
        })

        reasons = flag_assembly(hits)
        if reasons:
            flagged[name] = reasons

        status = f"{len(hits)} hits" + (f"  !! {'; '.join(reasons)}" if reasons else "")
        print(f"  {name}: {status}")

    summary_cols = ["assembly", "hit_num", "subject_seq", "t_start", "t_end",
                    "strand", "pident", "qcov", "mapq", "n_large_indels", "indel_details"]
    with open(args.outdir / "summary.tsv", "w") as f:
        f.write("\t".join(summary_cols) + "\n")
        for row in summary_rows:
            f.write("\t".join(str(row[c]) for c in summary_cols) + "\n")

    with open(args.outdir / "per_assembly_counts.tsv", "w") as f:
        f.write("assembly\ttotal_hits\thits_with_large_indels\n")
        for row in count_rows:
            f.write(f"{row['assembly']}\t{row['total_hits']}\t{row['hits_with_large_indels']}\n")

    with open(args.outdir / "flagged_assemblies.txt", "w") as f:
        for name, reasons in sorted(flagged.items()):
            f.write(f"{name}\t{'; '.join(reasons)}\n")

    print(f"\nSummary: {args.outdir}/summary.tsv")
    print(f"Counts:  {args.outdir}/per_assembly_counts.tsv")
    print(f"Flagged: {len(flagged)} assemblies → {args.outdir}/flagged_assemblies.txt")


# ── Single-sample detailed mode ───────────────────────────────────────────────

def run_detailed(args: argparse.Namespace) -> None:
    candidates = sorted(args.assembly_dir.glob(f"{args.sample}*.{args.assembly_ext}"))
    if not candidates:
        raise SystemExit(
            f"No assembly file matching '{args.sample}*.{args.assembly_ext}' "
            f"in {args.assembly_dir}"
        )
    assembly_fa = candidates[0]

    args.outdir.mkdir(parents=True, exist_ok=True)
    paf_path = args.outdir / f"{args.sample}_detailed.paf"

    print(f"Sample:          {args.sample}")
    print(f"Assembly:        {assembly_fa}")
    print(f"Query:           {args.query}")
    print(f"Preset:          {args.preset}")
    print(f"qcov min:        {args.detail_qcov_min}%")
    print(f"Indel min:       {args.detail_indel_min} bp")
    print(f"Outdir:          {args.outdir}\n")

    # Relaxed alignment: more secondaries, no min-score filtering
    paf = run_minimap2_raw(
        args.query, assembly_fa, args.preset,
        n_secondary=50,
        extra_flags=["-p", "0.1"],   # retain alignments down to 10% of best score
    )

    # Save raw PAF
    paf_path.write_text(paf)
    print(f"Raw PAF saved: {paf_path}\n")

    hits = parse_paf(paf, qcov_min=args.detail_qcov_min, min_mapq=0)

    if not hits:
        print("No hits found.")
        return

    # Sort by qcov descending
    hits.sort(key=lambda h: -h["qcov"])

    report_lines = []
    report_lines.append(f"Detailed alignment report — {args.sample}")
    report_lines.append(f"Query: {args.query}  ({hits[0]['q_len']:,} bp)")
    report_lines.append(f"Assembly: {assembly_fa}")
    report_lines.append(f"{'─' * 80}")
    report_lines.append(
        f"{'#':<4} {'contig':<30} {'strand':<7} {'t_start':>10} {'t_end':>10} "
        f"{'qcov':>7} {'pident':>7} {'mapq':>5} {'indels':>7}"
    )
    report_lines.append(f"{'─' * 80}")

    for k, h in enumerate(hits, 1):
        large = [d for d in h["indels"] if d["size"] >= args.detail_indel_min]
        report_lines.append(
            f"{k:<4} {h['t_name']:<30} {h['strand']:<7} "
            f"{h['t_start']:>10,} {h['t_end']:>10,} "
            f"{h['qcov']:>6.1f}% {h['pident']:>6.1f}% {h['mapq']:>5} "
            f"{len(large):>7}"
        )
        if large:
            for d in large:
                report_lines.append(
                    f"     {'':30}   {d['type']:>10} {d['size']:>6} bp"
                )

    report_lines.append(f"{'─' * 80}")
    report_lines.append(f"Total hits reported: {len(hits)}")
    flags = flag_assembly([h for h in hits if h["qcov"] >= QCOV_FILTER])
    if flags:
        report_lines.append(f"Flags (standard thresholds): {'; '.join(flags)}")
    else:
        report_lines.append("Flags (standard thresholds): none")

    report = "\n".join(report_lines)
    print(report)

    report_path = args.outdir / f"{args.sample}_detailed_report.txt"
    report_path.write_text(report + "\n")
    print(f"\nReport saved: {report_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # shared
    parser.add_argument("--query",                 required=True, type=Path)
    parser.add_argument("--assembly-dir",          required=True, type=Path)
    parser.add_argument("--assembly-ext",          default="fa")
    parser.add_argument("--outdir",                required=True, type=Path)
    parser.add_argument("--preset",                default="asm5")
    # batch mode
    parser.add_argument("--large-indel-threshold", default=50,  type=int,
                        help="Minimum indel size to report in batch mode (default: 50)")
    parser.add_argument("--min-mapq",              default=0,   type=int,
                        help="Minimum mapping quality in batch mode (default: 0)")
    # detailed mode
    parser.add_argument("--sample",                default=None,
                        help="Run detailed single-sample report for this sample ID")
    parser.add_argument("--detail-qcov-min",       default=20.0, type=float,
                        help="Minimum qcov%% in detailed mode (default: 20)")
    parser.add_argument("--detail-indel-min",      default=1,   type=int,
                        help="Minimum indel size to report in detailed mode (default: 1)")
    args = parser.parse_args()

    if args.sample:
        run_detailed(args)
    else:
        run_batch(args)


if __name__ == "__main__":
    main()
