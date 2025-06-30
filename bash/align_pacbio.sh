#!/bin/bash

# ========================================
# SLURM OPTIONS
# Edit these to match your HPC cluster
# ========================================

#SBATCH --mem=32g
#SBATCH --cpus-per-task=8
#SBATCH --array=1-138 # this should match number of your samples
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your-email@gmail.com

set -e

# ========================================
# align_pacbio.sh
#
# SLURM job script to align PacBio reads
# using minimap2 and convert to PSL format.
#
# Usage:
# - Edit the USER-DEFINED VARIABLES section
# - Submit with:
#     sbatch align_pacbio.sh
# ========================================

# === USER-DEFINED VARIABLES ===

OUTDIR="/path/to/output"
REFERENCE="/path/to/genome.fna"
GTF="/path/to/annotation.gtf"
CHR="/path/to/chrNameLength.txt"
FASTQDIR="/path/to/fastq"
FASTQLIST="/path/to/all_encode_pacbio.txt"

FLAIR_DIR="/path/to/flair/bin"

# ========================================

# Check required tools
command -v minimap2 >/dev/null 2>&1 || { echo >&2 "minimap2 not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools not found"; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "python not found"; exit 1; }

FASTQ=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FASTQLIST)
echo "Processing FASTQ: $FASTQ"

OUTFILE=$(echo $FASTQ | cut -d . -f 1)

module load minimap2/2.17
module load samtools/1.9

minimap2 -t 8 -ax splice:hq -uf $REFERENCE $FASTQDIR/$FASTQ > ${OUTDIR}/${OUTFILE}.sam

samtools view -hSb ${OUTDIR}/${OUTFILE}.sam > ${OUTDIR}/${OUTFILE}.all.bam
samtools view -q 60 -F 2304 -hb ${OUTDIR}/${OUTFILE}.all.bam > ${OUTDIR}/${OUTFILE}.bam
samtools sort ${OUTDIR}/${OUTFILE}.bam -o ${OUTDIR}/${OUTFILE}.sorted.bam
samtools index ${OUTDIR}/${OUTFILE}.sorted.bam

python ${FLAIR_DIR}/bam2Bed12 -i ${OUTDIR}/${OUTFILE}.sorted.bam > ${OUTDIR}/${OUTFILE}.sorted.bed
python ${FLAIR_DIR}/bed_to_psl $CHR ${OUTDIR}/${OUTFILE}.sorted.bed ${OUTDIR}/${OUTFILE}.psl
