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

mkdir -p "$PALINDROME_DIR"
mkdir -p "$BLACKLIST_DIR"

for fasta_file in "$SELF_ALIGN_DIR"/*.fasta; do

    file_name=$(basename "$fasta_file")
    contig_name=${file_name%.self_aling.fasta}
    sample_name=${contig_name%_random*}
    sample_name=${sample_name%_chrY}

    echo $sample_name

    output="${PALINDROME_DIR}/${contig_name}.pal"

    echo "${contig_name} found"

    if [ -f "${output}" ]; then

        echo "${contig_name} found"
        continue
    fi

    echo "run"

    blacklist_file="${BLACKLIST_DIR}/${sample_name}.bl.dat"

    if [ ! -f "$blacklist_file" ]; then
        zcat "${REPEAT_FILES}/${sample_name}/${sample_name}"*repmask-tblout-norm.bed.gz | \
            awk '{chrom=$1; start=$2; end=$3; class=$7; pattern="Sattellite/";
                 if ((class=="Low_complexity")||(class=="Simple_repeat")||(class=="Satellite")||(class~pattern))
                 print chrom,start,end,class;}' \
            > "$blacklist_file"
    fi

    python3 src/palindrover/palindrover.py \
            --minlength="${ARM_LEN}" \
            --minidentity=98% \
            --maxspacer="${SPACER_LEN}" \
            --blacklist:80%=${blacklist_file} \
            --column:palname \
            --group:overlaps \
            --report:blacklisted \
            --debug=blacklisted% \
            --blacklisted=${PALINDROME_DIR}/${contig_name}.bl \
            < $fasta_file > \
        "${output}"


    ((job_count++))
    if (( job_count % max_jobs == 0 )); then
        wait
    fi


done

wait 

echo "done"