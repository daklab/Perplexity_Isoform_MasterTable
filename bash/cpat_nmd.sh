#!/bin/bash

# ========================================
# cpat_nmd.sh
#
# SLURM job script to:
# - Convert PSL to FASTA
# - Run CPAT to assess coding potential
# - Translate ORFs to protein sequences
# - Clean FASTA headers
# - Filter NMD candidates
#
# Usage:
# - Edit the USER-DEFINED VARIABLES section
# - Submit with:
#     sbatch cpat_nmd.sh
#
# Requirements:
# - cpat.py
# - faTrans
# - Python scripts:
#     - psl_to_sequence.py
#     - filter_NMD.py
# - awk, tr, sed
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

# Input and output directories
INDIR="/path/to/sqanti_IR"
OUTDIR="/path/to/cpat_NMD"
CPATDIR="/path/to/cpat"

# Reference genome
REFERENCE="/path/to/GRCh38_no_alt_analysis_set.fna"

# Script paths
PSL_TO_SEQUENCE_SCRIPT="/path/to/psl_to_sequence.py"
FILTER_NMD_SCRIPT="/path/to/filter_NMD.py"

# ========================================

# Check required files
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE"
    exit 1
fi

mkdir -p ${OUTDIR}

# ========================================
# Convert PSL to FASTA
# ========================================

echo "Converting PSL to FASTA..."
python ${PSL_TO_SEQUENCE_SCRIPT} \
    ${INDIR}/all_samples_sp_collapse_all_chr_noIR.psl \
    ${REFERENCE} \
    ${INDIR}/all_samples_sp_collapse_all_chr_noIR.fa

# ========================================
# Run CPAT
# ========================================

echo "Running CPAT..."
cpat.py \
    -x ${CPATDIR}/Human_Hexamer.tsv \
    -d ${CPATDIR}/Human_logitModel.RData \
    -g ${INDIR}/all_samples_sp_collapse_all_chr_noIR.fa \
    --min-orf=90 \
    -o ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat \
    --top-orf=1

# ========================================
# Translate ORFs
# ========================================

echo "Translating ORFs..."
faTrans \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs.fa \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn.fa \
    -stop

# ========================================
# Clean FASTA headers
# ========================================

echo "Cleaning FASTA headers..."

# Collapse multi-line FASTA into one line per record
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' \
    < ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn.fa \
    | tr "\t" "\n" \
    > ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_linear.fa

# Modify FASTA headers: remove _ORF_x suffix and lowercase to uppercase
awk -F'_' 'BEGIN {OFS="_"} /^>/ {sub(/_ORF_[0-9]+$/, "", $0); print toupper($1), $2; next} {print}' \
    < ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_linear.fa \
    > ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_clean.fa

# Ensure consistent ENST capitalization
sed -i 's/>enst/>ENST/g' ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_clean.fa
sed -i 's/&enst/&ENST/g' ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_clean.fa

# ========================================
# Filter NMD
# ========================================

echo "Filtering NMD candidates..."

python ${FILTER_NMD_SCRIPT} \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_seqs_ptn_clean.fa \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_noNMD.ORF_seqs_ptn.fa \
    ${INDIR}/all_samples_sp_collapse_all_chr_noIR.psl \
    ${OUTDIR}/all_samples_sp_collapse_all_chr_noIR_cpat.ORF_prob.best.tsv

echo "DONE."

