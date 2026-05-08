#!/usr/bin/env bash
PALINDROME_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/palindromes"


source /storage/home/kxp5629/miniconda3/etc/profile.d/conda.sh
conda activate 10_Y

max_jobs=10
job_count=0

ARM_LEN="8K"
SPACER_LEN="500K"

PALINDROME_DIR="${PALINDROME_DIR}_A${ARM_LEN}_S${SPACER_LEN}"


compare_to="HG002"
# compare_to = "HG00261"
compare_to="HG02071"

for fasta_file in "$PALINDROME_DIR"/*.fasta; do
    echo $fasta_file
    contig=$(basename $fasta_file)
    contig=${contig%.left_pal.fasta}
    echo $contig

    #         	

    command="lastz --scores=lastz_scoring \
		    --allocate:traceback=2048M \
            --filter=coverage:80 \
             --format=general:name1,size1,length1,zstart1,end1,name2,strand2,size2,length2,zstart2+,end2+,id% \
           	$fasta_file[multiple] \
           	/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/${compare_to}_chrY.left_pal.fasta |
            ./src/add_cov_to_lastz.py - \
            > $PALINDROME_DIR/$contig-${compare_to}.lastz"
 
    while [ $(jobs -r | wc -l) -ge 16 ]; do
             wait -n
    done

    eval  "$command" &
done

wait

