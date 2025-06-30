#!/bin/bash

# ========================================
# collapse_pacbio.sh
#
# SLURM job script to:
# - Run isoform collapsing per chromosome
# - Add UTRs to collapsed isoforms
#
# Usage:
# - Edit the USER-DEFINED VARIABLES section
# - Submit with:
#     sbatch collapse_pacbio.sh
#
# Requirements:
# - Python scripts:
#     - sp_collapse_isoforms.py
#     - add_utr.py
# - bedtools
#
# Note:
# SLURM array indices:
#   1-22 → chromosomes 1-22
#   23   → chromosome X
#   24   → chromosome Y
# ========================================

#SBATCH --mem=128g
#SBATCH --cpus-per-task=12
#SBATCH --array=1-24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@domain.org

set -e

# ========================================
# USER-DEFINED VARIABLES
# ========================================

# Input and output directories
INDIR="/path/to/aligned_no_treatment"
OUTDIR="/path/to/collapsed_no_treatment"

# Python scripts
SP_COLLAPSE_SCRIPT="/path/to/sp_collapse_isoforms.py"
ADD_UTR_SCRIPT="/path/to/add_utr.py"

# Additional inputs
MANIFEST="/path/to/encode_pacbio_reads_manifest.tsv"
TSS_BED="/path/to/gencode.first_exon.bed"
TES_BED="/path/to/gencode.last_exon.bed"
GTF_DIR="/path/to/gencode_gtf_by_chr"

# ========================================

# Save original IFS
SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

# Determine chromosome name
CHR=$SLURM_ARRAY_TASK_ID

if [ "$CHR" == "23" ]; then
    CHR="X"
fi

if [ "$CHR" == "24" ]; then
    CHR="Y"
fi

echo "Processing chromosome: ${CHR}"

# Load bedtools (only if your cluster uses modules)
module load bedtools || true

# Create chromosome-specific output dir
mkdir -p "${OUTDIR}/chr${CHR}"

# Check required input files
IN_PSL="${INDIR}/all_pacbio_samples_chr${CHR}.psl"
GTF_FILE="${GTF_DIR}/gencode.v46.MScustom_chr${CHR}.gtf"

if [ ! -f "$IN_PSL" ]; then
    echo "ERROR: Input PSL file not found: $IN_PSL"
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file not found: $GTF_FILE"
    exit 1
fi

# Run sp_collapse_isoforms.py
python $SP_COLLAPSE_SCRIPT \
    $IN_PSL \
    --quantify \
    -tsv $MANIFEST \
    -tss $TSS_BED \
    -tes $TES_BED \
    -id \
    -n 10 \
    -gtf $GTF_FILE \
    -o ${OUTDIR}/chr${CHR}/all_samples_sp_collapse_chr${CHR}

# Run add_utr.py
python $ADD_UTR_SCRIPT \
    ${OUTDIR}/chr${CHR}/all_samples_sp_collapse_chr${CHR}_gid_read_map.txt \
    $IN_PSL \
    ${OUTDIR}/chr${CHR}/all_samples_sp_collapse_chr${CHR}_isoform_gid.psl \
    -o ${OUTDIR}/chr${CHR}/all_samples_sp_collapse_chr${CHR}

IFS=$SAVEIFS

echo "Done processing chromosome ${CHR}."

