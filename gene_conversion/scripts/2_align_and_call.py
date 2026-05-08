#!/usr/bin/env python3
"""Align palindrome arms to HG002 arm A reference and call SNP variants via CS tag."""
import os
import re
import pickle
import subprocess
from multiprocessing import Pool
from collections import defaultdict

PALINDROMES = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
EXCLUDE_SAMPLES = {"HG03456"}
WORK_DIR = "work"
REF_SAMPLE = "HG002"
NCPUS = 10

CS_OP = re.compile(r":[0-9]+|\*[a-z]{2}|\+[a-z]+|-[a-z]+")


def parse_cs(cs_str: str, ref_start: int) -> dict[int, tuple[str, str]]:
    """Return {ref_pos (0-based): (ref_allele, query_allele)} for SNPs only."""
    variants: dict[int, tuple[str, str]] = {}
    pos = ref_start
    for op in CS_OP.findall(cs_str):
        if op[0] == ":":
            pos += int(op[1:])
        elif op[0] == "*":
            variants[pos] = (op[1].upper(), op[2].upper())
            pos += 1
        elif op[0] == "+":
            pass  # insertion in query, ref pos unchanged
        elif op[0] == "-":
            pos += len(op) - 1  # deletion from query, ref advances
    return variants


def align_arm(args: tuple) -> tuple[str, str, str, dict]:
    """Align one arm to reference. Returns (palindrome, sample, arm, {pos: (ref, query)})."""
    palindrome, sample, arm, arm_fa, ref_fa = args

    result = subprocess.run(
        ["minimap2", "-x", "asm5", "--cs", "-c", ref_fa, arm_fa],
        capture_output=True, text=True,
    )

    variants: dict[int, tuple[str, str]] = {}
    for line in result.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        ref_start = int(fields[7])
        cs_tag = next((f for f in fields if f.startswith("cs:Z:")), None)
        if cs_tag is None:
            continue
        variants.update(parse_cs(cs_tag[5:], ref_start))

    return palindrome, sample, arm, variants


def build_tasks(palindrome: str) -> list[tuple]:
    arms_dir = os.path.join(WORK_DIR, palindrome, "arms")
    ref_fa = os.path.join(arms_dir, f"{REF_SAMPLE}_A.fa")
    tasks = []
    for fname in sorted(os.listdir(arms_dir)):
        if not fname.endswith(".fa"):
            continue
        sample, arm = fname[:-3].rsplit("_", 1)
        if sample in EXCLUDE_SAMPLES:
            continue
        arm_fa = os.path.join(arms_dir, fname)
        tasks.append((palindrome, sample, arm, arm_fa, ref_fa))
    return tasks


def main() -> None:
    all_tasks: list[tuple] = []
    for palindrome in PALINDROMES:
        all_tasks.extend(build_tasks(palindrome))

    print(f"Running {len(all_tasks)} alignments with {NCPUS} workers...")
    with Pool(NCPUS) as pool:
        results = pool.map(align_arm, all_tasks)

    # Group by palindrome → {sample: {arm: variants}}
    by_palindrome: dict[str, dict] = defaultdict(lambda: defaultdict(dict))
    for palindrome, sample, arm, variants in results:
        by_palindrome[palindrome][sample][arm] = variants

    for palindrome, samples in by_palindrome.items():
        out = os.path.join(WORK_DIR, palindrome, "variants", "variants.pkl")
        with open(out, "wb") as f:
            pickle.dump(dict(samples), f)
        print(f"{palindrome}: {len(samples)} samples saved to {out}")

    print("Done.")


if __name__ == "__main__":
    main()
