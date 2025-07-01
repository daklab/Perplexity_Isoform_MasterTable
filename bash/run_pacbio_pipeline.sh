#!/bin/bash

###############################################################
# run_pacbio_pipeline.sh
#
# Master pipeline runner for Stella Park's PacBio pipeline
#
# Runs in sequence:
#   1. align_pacbio.sh
#   2. before_collapse.sh
#   3. collapse_pacbio.sh
#   4. sqanti_IR.sh
#   5. cpat_nmd.sh
#   6. generate_master_table.R
#
# User controls:
#   SUBMIT="sbatch"   # run via SLURM
#   SUBMIT="bash"     # run locally
#
# Author: Stella Park (script generated)
###############################################################

# -------------------------------
# User Settings
# -------------------------------

# Choose how to run the scripts:
#   sbatch → submit jobs to Slurm
#   bash   → run locally (no job scheduler)
SUBMIT="sbatch"

# Path to scripts
SCRIPT_DIR="/path/to/your/scripts"

# Path to Rscript binary
RSCRIPT_BIN="/usr/bin/Rscript"

# Path to generate_master_table.R
R_SCRIPT="${SCRIPT_DIR}/generate_master_table.R"

# -------------------------------
# Run functions
# -------------------------------

run_step() {
    script_name="$1"

    echo "============================================================"
    echo "Running step: $script_name"
    echo "============================================================"

    if [[ "$SUBMIT" == "sbatch" ]]; then
        JOB_ID=$($SUBMIT ${SCRIPT_DIR}/${script_name})
        echo "Submitted Slurm job $JOB_ID for $script_name"
    else
        bash ${SCRIPT_DIR}/${script_name}
    fi

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Script $script_name failed. Exiting pipeline."
        exit 1
    fi

    echo "Step $script_name completed successfully."
    echo
}

# -------------------------------
# Run pipeline steps
# -------------------------------

# 1. align_pacbio.sh
run_step "align_pacbio.sh"

# 2. before_collapse.sh
run_step "before_collapse.sh"

# 3. collapse_pacbio.sh
run_step "collapse_pacbio.sh"

# 4. sqanti_IR.sh
run_step "sqanti_IR.sh"

# 5. cpat_nmd.sh
run_step "cpat_nmd.sh"

# 6. generate_master_table.R
echo "============================================================"
echo "Running final R script: generate_master_table.R"
echo "============================================================"

if [[ "$SUBMIT" == "sbatch" ]]; then
    # Wrap Rscript into a Slurm job submission
    sbatch --wrap="${RSCRIPT_BIN} ${R_SCRIPT}"
else
    ${RSCRIPT_BIN} ${R_SCRIPT}
fi

if [[ $? -ne 0 ]]; then
    echo "ERROR: R script generate_master_table.R failed. Exiting pipeline."
    exit 1
fi

echo
echo "============================================================"
echo "All pipeline steps completed successfully!"
echo "Master table generated."
echo "============================================================"

