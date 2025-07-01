#!/bin/bash

# ========================================
# sqanti_IR.sh
#
# SLURM job script to:
# - Concatenate PSL files from multiple chromosomes
# - Convert PSL to GTF
# - Run SQANTI3 QC
# - Filter intron retention events
#
# Usage:
# - Edit the USER-DEFINED VARIABLES section
# - Submit with:
#     sbatch sqanti_IR.sh
#
# Requirements:
# - conda environment with SQANTI3 installed
# - Python scripts:
#     - psl_to_gtf_sp.py
#     - filter_IR.py
#
# Note:
# - Adjust SLURM directives as needed for your cluster.
# ========================================

#SBATCH --mem=16g
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=youremail@domain.org

set -e

# ========================================
# USER-DEFINED VARIABLES
# ========================================

# Input / output directories
COLDIR="/path/to/collapsed_no_treatment"
OUTDIR="/path/to/sqanti_IR_no_treatment"

# Reference genome and annotation
REFERENCE="/path/to/GRCh38_no_alt_analysis_set.fna"
GTF="/path/to/gencode.v46.MScustom.gtf"

# Scripts and environments
PSL_TO_GTF_SCRIPT="/path/to/psl_to_gtf_sp.py"
FILTER_IR_SCRIPT="/path/to/filter_IR.py"
SQANTI3_ENV="SQANTI3.env"
PACBIO_ENV="encode_pacbio"

# SQANTI3 executable
SQANTI3_QC="/path/to/sqanti3_qc.py"

# ========================================

# Check required files
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "ERROR: Annotation GTF not found: $GTF"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p ${OUTDIR}

echo "Concatenating PSL files..."
for i in {1..22} X Y; do
    if [ -f ${COLDIR}/chr${i}/all_samples_sp_collapse_chr${i}_no_treatment_isoform_full.psl ]; then
        cat ${COLDIR}/chr${i}/all_samples_sp_collapse_chr${i}_no_treatment_isoform_full.psl
    else
        echo "Warning: PSL file missing for chr${i}, skipping."
    fi
done > ${COLDIR}/all_samples_sp_collapse_all_chr_no_treatment_full.psl

echo "Converting PSL to GTF..."
python ${PSL_TO_GTF_SCRIPT} \
    ${COLDIR}/all_samples_sp_collapse_all_chr_no_treatment_full.psl \
> ${COLDIR}/all_samples_sp_collapse_all_chr_no_treatment_full.gtf

echo "Activating SQANTI3 environment..."
eval "$(conda shell.bash hook)"
conda activate ${SQANTI3_ENV}

echo "Running SQANTI3 QC..."
python ${SQANTI3_QC} \
    ${COLDIR}/all_samples_sp_collapse_all_chr_no_treatment_full.gtf \
    ${GTF} \
    ${REFERENCE} \
    --min_ref_len 0 \
    --force_id_ignore \
    --aligner_choice minimap2 \
    --skipORF \
    -o all_samples_all_chr_no_treatment_sqanti3_qc \
    -d ${OUTDIR} \
    --report pdf

echo "Switching to PacBio environment..."
conda activate ${PACBIO_ENV}

echo "Filtering intron retention events..."
python ${FILTER_IR_SCRIPT} \
    ${OUTDIR}/all_samples_all_chr_no_treatment_sqanti3_qc_classification.txt \
    ${COLDIR}/all_samples_sp_collapse_all_chr_no_treatment_full.psl \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_no_treatment_noIR.psl

echo "DONE."
