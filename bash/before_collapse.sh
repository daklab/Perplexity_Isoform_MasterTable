#!/bin/bash

# ========================================
# before_collapse.sh
#
# SLURM job script to:
# - Split a GTF file into separate files per chromosome
# - Merge PSL files across samples
# - Split merged PSL file into separate files per chromosome
#
# Usage:
# - Edit the USER-DEFINED VARIABLES section
# - Submit with:
#     sbatch before_collapse.sh
#
# Requirements:
# - awk
# ========================================

#SBATCH --mem=32g
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@domain.org

set -e

# ========================================
# USER-DEFINED VARIABLES
# ========================================

INDIR="/path/to/aligned_output"
GTF="/path/to/gencode.v46.MScustom.gtf"
OUTDIR="/path/to/output"

# ========================================

# Check input GTF exists
if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF file not found at $GTF"
    exit 1
fi

# Create output directory if needed
mkdir -p ${OUTDIR}

echo "Splitting GTF file by chromosome..."
awk -F'\t' '{print > "'${OUTDIR}'/"$1".gtf"}' ${GTF}

# Rename files to follow chr naming
for i in {1..22} X Y; do
    if [ -f ${OUTDIR}/chr${i}.gtf ]; then
        mv ${OUTDIR}/chr${i}.gtf ${OUTDIR}/gencode.v46.MScustom_chr${i}.gtf
        echo "Created: ${OUTDIR}/gencode.v46.MScustom_chr${i}.gtf"
    fi
done

echo "Merging PSL files..."
cat ${INDIR}/ENCFF*.psl > ${OUTDIR}/all_samples_all_reads.psl

echo "Splitting merged PSL by chromosome..."
awk -F'\t' '{print > "'${OUTDIR}'/"$14".psl"}' ${OUTDIR}/all_samples_all_reads.psl

for i in {1..22} X Y; do
    if [ -f ${OUTDIR}/chr${i}.psl ]; then
        mv ${OUTDIR}/chr${i}.psl ${OUTDIR}/all_pacbio_samples_chr${i}.psl
        echo "Created: ${OUTDIR}/all_pacbio_samples_chr${i}.psl"
    fi
done

echo "DONE."

