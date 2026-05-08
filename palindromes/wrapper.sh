#!/usr/bin/env bash

set -x 

DATA_DIR="data/palindromes_A8K_S500K"

# ## run self alignment with lastZ
# 
# bash src/palindromes/01_self_align.sh

# ## run palindrover
# bash src/palindromes/03_run_palindrome.sh

# ## collect palindromes into a BED file
# BED_DIR="${DATA_DIR}/BED"
# 
# mkdir -p "$BED_DIR"
# python src/palindromes/04_collect.py -p "$DATA_DIR" | bedtools sort -i - > "${BED_DIR}/HPRC_ALL.bed"

# ## extract left arms of palindromes
# bash src/palindromes/06_extract_fasta.sh

# ## aling all left arms to one "reference" assembly tried with HG002 but that has the less common
# ## inversion where palindrome 1 doesn't really exist and P2 is long using HG00621 as reference with long P1
# bash src/palindromes/07_align_left_pals_HG002.sh

ref="HG02071"
match_table="output/output.table_all_v_${ref}_w_unmatched"

# ## create table output from single ref aligned left palindromes
# python "src/palindromes/08b_create_table.py" \
#     "data/palindromes_A8K_S500K" \
#     $ref > $match_table

## create stats from alignments to ref
# python src/palindromes/10.1_extract_metrics_hg002vall_split_multicopy.py \
#    -i $match_table \
#    -p $DATA_DIR \
#    -o out_stats.png

# ## align all vs all for better extraction of homologous arms
# sbatch src/palindromes/09_merge_fasta_algn_all_v_all_sbatch.sh


python src/palindromes/13d_extract_metrics.py -i output_elements_table.csv -p data/palindromes_A8K_S500K -o output.pdf
