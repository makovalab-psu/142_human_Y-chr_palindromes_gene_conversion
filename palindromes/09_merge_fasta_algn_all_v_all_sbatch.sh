#!/bin/bash
# filepath: /storage/home/kxp5629/proj/10_HPRC_R2_Y/09_merge_fasta_algn_all_v_all_array.sbatch

#SBATCH --job-name=lastz_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=open
#SBATCH --array=1-795%20

# Set up directories
FASTAS_DIR="/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K"
ALL_FASTA_DIR="${FASTAS_DIR}/all_v_all"
ALL_FASTA="$ALL_FASTA_DIR/all.left_pal.fasta"
IND_ALG="${ALL_FASTA_DIR}/individual_alignments"

mkdir -p "${ALL_FASTA_DIR}"
mkdir -p "$IND_ALG"

# Activate conda environment
source /storage/home/kxp5629/miniconda3/etc/profile.d/conda.sh
conda activate 10_Y

if [ ! -f "$ALL_FASTA" ]; then
    cat "$FASTAS_DIR"/*left_pal.fasta > "$ALL_FASTA"
fi


# Get list of fasta files and select the one for this array task
FASTA_FILES=("$FASTAS_DIR"/*left_pal.fasta)
fasta="${FASTA_FILES[$((SLURM_ARRAY_TASK_ID-1))]}"

if [ ! -f "$fasta" ]; then
    echo "No fasta file for task ID $SLURM_ARRAY_TASK_ID"
    exit 0
fi

sample=$(basename "$fasta")
sample="${sample%.left_pal.fasta}"
output="$IND_ALG/$sample.to_all.lastz"

if [ -f "$output" ]; then
    echo "Output already exists for $sample"
    exit 0
fi

echo "Processing $sample at $(date)"

lastz "$fasta"[multiple] "$ALL_FASTA" \
     --identity=90 \
     --hspthresh=5000 \
     --coverage=70 \
     --step=10 \
     --allocate:traceback=2048M \
     --format=general:name1,size1,length1,zstart1,end1,name2,strand2,size2,length2,zstart2+,end2+,id% \
     --output="$output"

echo "Completed $sample at $(date)"