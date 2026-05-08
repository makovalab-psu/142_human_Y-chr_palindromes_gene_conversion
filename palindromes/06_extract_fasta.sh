#!/usr/bin/env bash

REPEAT_FILES="/storage/home/kxp5629/group_storage/shared/HPRC_R2_Y/v2/verkko-v2.2.1/annotations/repeatmasker"
MASKED_FILES="/storage/home/kxp5629/group_storage/shared/HPRC_R2_Y/masked/masked_sequences"
SELF_ALIGN_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/self_align"

TEMP_SEQLIST_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/self_align/SEQLIST_DIR"
PALINDROME_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/palindromes"

BLACKLIST_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/blacklist"



source /storage/home/kxp5629/miniconda3/etc/profile.d/conda.sh
# conda activate 10_Y

max_jobs=10
job_count=0

ARM_LEN="8K"
SPACER_LEN="500K"

PALINDROME_DIR="${PALINDROME_DIR}_A${ARM_LEN}_S${SPACER_LEN}"


for palindrome_file in "$PALINDROME_DIR"/*.pal; do
    echo $palindrome_file
    contig=$(basename "$palindrome_file")
    sample="${contig%%_*}"
    echo $sample
    bed_file=${palindrome_file%.pal}.left.bed
    cut -f1,2,3,10 $palindrome_file > $bed_file
    fasta_file="${MASKED_FILES}/${sample}_masked.fa"
    bedtools getfasta -fi $fasta_file -bed $bed_file -name >${palindrome_file%.pal}.left_pal.fasta
done