#!/usr/bin/env python3
"""Identify putative conversion tracts: co-occurring GC events on the same node within 1kb."""
import os
from collections import defaultdict

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
WORK_DIR = "work"
MAX_DIST = 1000  # bp


def load_gc_events(events_file: str) -> dict[str, list[int]]:
    """Return {node: [sorted positions]} for geneconversion events only."""
    node_positions: dict[str, list[int]] = defaultdict(list)
    with open(events_file) as f:
        next(f)  # header
        for line in f:
            fields = line.strip().split("\t")
            if fields[6] != "geneconversion":
                continue
            pos = int(fields[1])
            node = fields[7]
            node_positions[node].append(pos)
    return {node: sorted(positions) for node, positions in node_positions.items()}


def cluster_positions(positions: list[int], max_dist: int) -> list[list[int]]:
    """Group sorted positions into clusters where consecutive gap <= max_dist."""
    if not positions:
        return []
    clusters = [[positions[0]]]
    for pos in positions[1:]:
        if pos - clusters[-1][-1] <= max_dist:
            clusters[-1].append(pos)
        else:
            clusters.append([pos])
    return [c for c in clusters if len(c) >= 2]


def get_individuals(events_file: str, palindrome: str, node: str, positions: list[int]) -> str:
    """Get individuals for the first position in the tract (all should be same node)."""
    pos_set = set(positions)
    with open(events_file) as f:
        next(f)
        for line in f:
            fields = line.strip().split("\t")
            if fields[6] == "geneconversion" and fields[7] == node and int(fields[1]) in pos_set:
                return fields[8]
    return ""


def main() -> None:
    out_file = os.path.join(WORK_DIR, "conversion_tracts.tsv")

    all_tracts = []
    for palindrome in PALINDROMES:
        events_file = os.path.join(WORK_DIR, palindrome, "events", f"{palindrome}.events.txt")
        node_positions = load_gc_events(events_file)

        for node, positions in node_positions.items():
            tracts = cluster_positions(positions, MAX_DIST)
            for tract in tracts:
                individuals = get_individuals(events_file, palindrome, node, tract)
                n_individuals = len(individuals.split(",")) if individuals else 0
                all_tracts.append({
                    "palindrome": palindrome,
                    "node": node,
                    "n_individuals": n_individuals,
                    "individuals": individuals,
                    "n_sites": len(tract),
                    "start": tract[0],
                    "end": tract[-1],
                    "span_bp": tract[-1] - tract[0],
                    "positions": ",".join(map(str, tract)),
                })

    # Sort by span descending
    all_tracts.sort(key=lambda x: (-x["n_sites"], -x["span_bp"]))

    with open(out_file, "w") as f:
        f.write("palindrome\tnode\tn_individuals\tn_sites\tstart\tend\tspan_bp\tpositions\tindividuals\n")
        for t in all_tracts:
            f.write(
                f"{t['palindrome']}\t{t['node']}\t{t['n_individuals']}\t{t['n_sites']}\t"
                f"{t['start']}\t{t['end']}\t{t['span_bp']}\t{t['positions']}\t{t['individuals']}\n"
            )

    print(f"Found {len(all_tracts)} putative conversion tracts (max_dist={MAX_DIST}bp)")
    print(f"Written to {out_file}\n")

    # Summary per palindrome
    print(f"{'palindrome':<12} {'tracts':>7} {'max_sites':>10} {'max_span_bp':>12}")
    print("-" * 45)
    for pal in PALINDROMES:
        pal_tracts = [t for t in all_tracts if t["palindrome"] == pal]
        if not pal_tracts:
            print(f"{pal:<12} {'0':>7}")
            continue
        max_sites = max(t["n_sites"] for t in pal_tracts)
        max_span = max(t["span_bp"] for t in pal_tracts)
        print(f"{pal:<12} {len(pal_tracts):>7} {max_sites:>10} {max_span:>12}")


if __name__ == "__main__":
    main()
