# Human Y Chromosome Palindromes & Gene Conversion

Computational pipeline to identify palindromic regions on the human Y chromosome and characterize palindrome aided gene conversion events within.

## Overview

Two pipelines:

**`palindromes/`** — Discovers and catalogs Y chromosome palindromes across samples.
1. Self-aligns masked Y chromosome sequences with LASTZ to find inverted repeats.
2. Parses palindrome detector output to extract arm A/B boundaries.
3. Clusters palindromes into homologous groups via alignment-coverage graphs.
4. Computes pairwise sequence identities and renders heatmaps.

**`gene_conversion/`** — Detects and analyzes gene conversion events within palindrome arms.
1. Extracts per-arm sequences from masked FASTAs using BED coordinates.
2. Aligns arms to a reference (HG002) and calls variants to produce a multi-sample VCF.
3. Applies phylogenetic parsimony to classify each variant as a point mutation or gene conversion event.
4. Clusters conversion events into tracts, estimates tract lengths, and quantifies AT→GC bias.

## Dependencies

**External tools:** LASTZ, samtools

**Python packages** (see `gene_conversion/requirements.txt`):
```
biopython, networkx, pandas, matplotlib, seaborn, ete3
```

## Usage

```bash
# Run the full gene conversion pipeline
bash gene_conversion/run_pipeline.sh
```

Individual scripts are numbered sequentially within each sub-pipeline and can be run step-by-step.

## Input Data

- Masked Y chromosome FASTA files per sample (`{sample}_masked.fa`)
- BED file of palindrome arm coordinates
- Species/population tree topology (used for parsimony-based event calling)
