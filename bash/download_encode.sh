#!/bin/bash

# Download ENCODE4 PacBio FASTQ files
# Usage: bash download_encode_data.sh /path/to/output_dir

set -euo pipefail

OUTDIR="$1"
if [[ -z "$OUTDIR" ]]; then
    echo "Usage: bash download_encode_data.sh /path/to/fastq_dir"
    exit 1
fi

ACCESSIONS="$(dirname "${BASH_SOURCE[0]}")/../reference_files/encode4_lrs_accessions.txt"

mkdir -p "$OUTDIR"

while IFS= read -r ACC; do
    [[ -z "$ACC" ]] && continue
    if [[ -f "${OUTDIR}/${ACC}.fastq.gz" ]]; then
        echo "SKIP: $ACC"
    else
        echo "Downloading: $ACC"
        wget -q -O "${OUTDIR}/${ACC}.fastq.gz" \
            "https://www.encodeproject.org/files/${ACC}/@@download/${ACC}.fastq.gz"
    fi
done < "$ACCESSIONS"

echo "Done. Downloaded to: $OUTDIR"