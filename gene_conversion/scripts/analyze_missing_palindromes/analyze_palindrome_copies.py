#!/usr/bin/env python3
"""
Palindrome copy number & indel analysis via BLAST.
Usage:
    python analyze_palindrome_copies.py \\
        --query QUERY.fa \\
        --assembly-dir /path/to/assemblies \\
        --assembly-ext fa \\
        --outdir results/ \\
        --large-indel-threshold 50 \\
        [--keep-db]
"""
import argparse
import os
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


BLAST_FMT = "6 qseqid sseqid pident length qlen qstart qend sstart send gapopen gaps btop"
QCOV_FILTER   = 80.0   # minimum qcov to keep a hit
QCOV_FLAG     = 90.0   # flag if any hit below this
OVERLAP_FRAC  = 0.50   # flag if same-contig hits overlap > this fraction of qlen


# ── BTOP parsing ─────────────────────────────────────────────────────────────

def parse_btop_indels(btop: str, threshold: int) -> list[dict]:
    """
    Extract large indels from a BTOP string.
    Integers = matched bases; letter pairs = substitutions;
    '-' runs = gaps (gap in subject = insertion in query; gap in query = deletion from query).
    Returns list of {type, size} for gaps >= threshold.
    """
    indels = []
    # tokenise: integers, letter pairs, or dash-runs
    for token in re.findall(r"\d+|[A-Za-z]{2}|-+", btop):
        if token.startswith("-"):
            size = len(token)
            if size >= threshold:
                # even length: alternating subject/query gaps per BLAST BTOP encoding
                # odd length should not occur; treat conservatively as mixed
                if size % 2 == 0:
                    # gaps in subject (insertions in query) and gaps in query (deletions)
                    # BTOP encodes subject gap as '-X' and query gap as 'X-'
                    # but in practice a run of dashes means consecutive gaps
                    indels.append({"type": "gap_run", "size": size})
                else:
                    indels.append({"type": "gap_run", "size": size})
    return indels


def parse_btop_indels_detailed(btop: str, threshold: int) -> list[dict]:
    """
    Parse BTOP token by token to distinguish insertions (gap in subject, '-N')
    from deletions (gap in query, 'N-').
    Each substitution/gap pair is two characters; integers are match runs.
    """
    indels = []
    i = 0
    while i < len(btop):
        # integer run
        m = re.match(r"\d+", btop[i:])
        if m:
            i += len(m.group())
            continue
        # two-character token
        if i + 1 < len(btop):
            pair = btop[i:i+2]
            if pair[0] == "-":          # gap in subject = insertion in query
                # count consecutive subject-gap pairs
                run = 0
                j = i
                while j + 1 < len(btop) and btop[j] == "-" and btop[j+1] != "-":
                    run += 1
                    j += 2
                if run >= threshold:
                    indels.append({"type": "insertion", "size": run})
                i = j
                continue
            elif pair[1] == "-":        # gap in query = deletion from query
                run = 0
                j = i
                while j + 1 < len(btop) and btop[j] != "-" and btop[j+1] == "-":
                    run += 1
                    j += 2
                if run >= threshold:
                    indels.append({"type": "deletion", "size": run})
                i = j
                continue
            else:
                i += 2  # substitution pair
                continue
        i += 1
    return indels


# ── BLAST ─────────────────────────────────────────────────────────────────────

def make_blast_db(assembly_fa: Path, db_dir: Path) -> Path:
    db_prefix = db_dir / assembly_fa.stem
    subprocess.run(
        ["makeblastdb", "-in", str(assembly_fa), "-dbtype", "nucl",
         "-out", str(db_prefix)],
        check=True, capture_output=True,
    )
    return db_prefix


def run_blastn(query: Path, db_prefix: Path, out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "blastn", "-task", "blastn",
            "-query", str(query),
            "-db", str(db_prefix),
            "-dust", "no",
            "-perc_identity", "95",
            "-word_size", "11",
            "-outfmt", BLAST_FMT,
            "-out", str(out_tsv),
        ],
        check=True, capture_output=True,
    )


def remove_blast_db(db_prefix: Path) -> None:
    for ext in (".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto"):
        p = Path(str(db_prefix) + ext)
        if p.exists():
            p.unlink()


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_blast_hits(tsv: Path, threshold: int, qcov_min: float) -> list[dict]:
    hits = []
    if not tsv.exists() or tsv.stat().st_size == 0:
        return hits
    with open(tsv) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            (qseqid, sseqid, pident, length, qlen,
             qstart, qend, sstart, send, gapopen, gaps, btop) = parts

            pident  = float(pident)
            qlen    = int(qlen)
            qstart  = int(qstart)
            qend    = int(qend)
            sstart  = int(sstart)
            send    = int(send)
            qcov    = (qend - qstart + 1) / qlen * 100
            strand  = "minus" if sstart > send else "plus"

            if qcov < qcov_min:
                continue

            indels = parse_btop_indels_detailed(btop, threshold)
            hits.append({
                "sseqid":         sseqid,
                "sstart":         sstart,
                "send":           send,
                "strand":         strand,
                "pident":         pident,
                "qcov":           qcov,
                "qlen":           qlen,
                "n_large_indels": len(indels),
                "indel_details":  ";".join(f"{d['type']}:{d['size']}" for d in indels) or ".",
            })
    return hits


# ── Flagging ──────────────────────────────────────────────────────────────────

def flag_assembly(assembly: str, hits: list[dict], qlen: int) -> list[str]:
    reasons = []

    if len(hits) != 2:
        reasons.append(f"hit_count={len(hits)}")

    for h in hits:
        if h["qcov"] < QCOV_FLAG:
            reasons.append(f"low_qcov={h['qcov']:.1f}%_on_{h['sseqid']}")

    if len(hits) == 2:
        strands = [h["strand"] for h in hits]
        if strands[0] == strands[1]:
            reasons.append(f"same_strand={strands[0]}")

        if hits[0]["sseqid"] == hits[1]["sseqid"]:
            s0 = (min(hits[0]["sstart"], hits[0]["send"]),
                  max(hits[0]["sstart"], hits[0]["send"]))
            s1 = (min(hits[1]["sstart"], hits[1]["send"]),
                  max(hits[1]["sstart"], hits[1]["send"]))
            overlap = max(0, min(s0[1], s1[1]) - max(s0[0], s1[0]))
            if overlap > OVERLAP_FRAC * qlen:
                reasons.append(f"overlapping_hits_on_{hits[0]['sseqid']}")

    return reasons


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query",                required=True, type=Path)
    parser.add_argument("--assembly-dir",         required=True, type=Path)
    parser.add_argument("--assembly-ext",         default="fa")
    parser.add_argument("--outdir",               required=True, type=Path)
    parser.add_argument("--large-indel-threshold",default=50, type=int)
    parser.add_argument("--keep-db",              action="store_true")
    args = parser.parse_args()

    blast_raw_dir = args.outdir / "blast_raw"
    db_dir        = args.outdir / "blast_db"
    blast_raw_dir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)

    assemblies = sorted(args.assembly_dir.glob(f"*.{args.assembly_ext}"))
    if not assemblies:
        raise SystemExit(f"No *.{args.assembly_ext} files found in {args.assembly_dir}")

    print(f"Query:      {args.query}")
    print(f"Assemblies: {len(assemblies)}")
    print(f"Outdir:     {args.outdir}")
    print(f"Indel threshold: {args.large_indel_threshold} bp\n")

    summary_rows:  list[dict] = []
    count_rows:    list[dict] = []
    flagged:       dict[str, list[str]] = {}
    db_prefixes:   list[Path] = []

    for assembly_fa in assemblies:
        name = assembly_fa.stem
        raw_tsv    = blast_raw_dir / f"{name}.tsv"
        db_prefix  = make_blast_db(assembly_fa, db_dir)
        db_prefixes.append(db_prefix)

        run_blastn(args.query, db_prefix, raw_tsv)
        hits = parse_blast_hits(raw_tsv, args.large_indel_threshold, QCOV_FILTER)

        qlen = hits[0]["qlen"] if hits else 0

        if not hits:
            summary_rows.append({
                "assembly": name, "hit_num": 0,
                "subject_seq": ".", "sstart": ".", "send": ".", "strand": ".",
                "pident": ".", "qcov": ".", "n_large_indels": 0, "indel_details": ".",
            })
        else:
            for k, h in enumerate(hits, 1):
                summary_rows.append({
                    "assembly":      name,
                    "hit_num":       k,
                    "subject_seq":   h["sseqid"],
                    "sstart":        h["sstart"],
                    "send":          h["send"],
                    "strand":        h["strand"],
                    "pident":        f"{h['pident']:.2f}",
                    "qcov":          f"{h['qcov']:.2f}",
                    "n_large_indels":h["n_large_indels"],
                    "indel_details": h["indel_details"],
                })

        n_with_indels = sum(1 for h in hits if h["n_large_indels"] > 0)
        count_rows.append({
            "assembly":             name,
            "total_hits":           len(hits),
            "hits_with_large_indels": n_with_indels,
        })

        reasons = flag_assembly(name, hits, qlen)
        if reasons:
            flagged[name] = reasons

        status = f"{len(hits)} hits" + (f"  !! {'; '.join(reasons)}" if reasons else "")
        print(f"  {name}: {status}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    summary_cols = ["assembly", "hit_num", "subject_seq", "sstart", "send",
                    "strand", "pident", "qcov", "n_large_indels", "indel_details"]
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

    print(f"\nSummary:   {args.outdir}/summary.tsv")
    print(f"Counts:    {args.outdir}/per_assembly_counts.tsv")
    print(f"Flagged:   {len(flagged)} assemblies → {args.outdir}/flagged_assemblies.txt")

    # ── Cleanup ───────────────────────────────────────────────────────────────
    if not args.keep_db:
        for db_prefix in db_prefixes:
            remove_blast_db(db_prefix)
        shutil.rmtree(db_dir, ignore_errors=True)
        print("BLAST db files removed.")


if __name__ == "__main__":
    main()
