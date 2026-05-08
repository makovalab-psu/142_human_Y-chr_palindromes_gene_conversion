#!/usr/bin/env bash
# Full HPRC-Y gene conversion pipeline.
# All parameters are set inside each script; this wrapper just runs them in order.
set -euo pipefail

PYTHON=".venv/bin/python"
SCRIPTS="scripts"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Step 1: Extract palindrome arm sequences"
$PYTHON $SCRIPTS/1_extract_arms.py

log "Step 2: Align arms to HG002 reference and call variants"
$PYTHON $SCRIPTS/2_align_and_call.py

log "Step 3: Build VCF (arm A vs arm B genotypes, min 80% sample coverage)"
$PYTHON $SCRIPTS/3_build_vcf.py

log "Step 4: Prune species tree to per-palindrome sample sets"
$PYTHON $SCRIPTS/4_prep_tree.py

log "Step 5: Parsimony-based classification of mutations vs gene conversion events"
$PYTHON $SCRIPTS/5_find_events.py

log "Step 6: Per-site phylogeny plots (Graphviz)"
$PYTHON $SCRIPTS/6_plot.py

log "Step 7: Identify putative conversion tracts (≥2 GC events on same node within 1kb)"
$PYTHON $SCRIPTS/7_conversion_tracts.py

log "Step 8: GC analysis and summary figure"
$PYTHON $SCRIPTS/8_gc_analysis.py

log "Step 9: Per-site parsimony tree plots (PDF)"
$PYTHON $SCRIPTS/9_plot_trees.py

log "Pipeline complete."
