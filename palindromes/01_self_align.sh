MASKED_FILES="/storage/home/kxp5629/group_storage/shared/HPRC_R2_Y/masked/masked_sequences"
SELF_ALIGN_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/self_align"

TEMP_SEQLIST_DIR="/storage/group/kdm16/default/kxp5629/proj/10_HPRC_R2_Y/data/self_align/SEQLIST_DIR"

mkdir -p "$TEMP_SEQLIST_DIR"

source /storage/home/kxp5629/miniconda3/etc/profile.d/conda.sh
conda activate 10_Y

for fasta_file in "$MASKED_FILES"/*.fa; do

    file_name=$(basename "$fasta_file")
    sample_name=${file_name%_masked.fa}
    echo $sample_name

    mapfile -t results < <(grep ">" "$fasta_file")

    for seq_id in "${results[@]}"; do
        echo "    Processing sequence: ${seq_id#>}"
        # Add your processing logic here

        seq_file="$TEMP_SEQLIST_DIR/${seq_id#>}.fasta"
        echo "${seq_id#>}" > "$seq_file"

        output_file="${SELF_ALIGN_DIR}/${seq_id#>}.self_aling.fasta"

        echo "$output_file"

        if [ -f $output_file ]; then
          continue
        fi

        command="lastz --scores=lastz_scoring \
		--allocate:traceback=2048M \
           	--format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
           	${fasta_file}[subset=${seq_file}] \
           	${fasta_file}[subset=${seq_file}] \
           	--strand=minus > ${output_file}"

        echo $command

        eval  "$command &"
         while [ $(jobs -r | wc -l) -ge 8 ]; do
              wait -n
         done
    done



done

wait

echo "done"
