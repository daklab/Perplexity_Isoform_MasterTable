# Perplexity_Isoform_MasterTable

Bash scripts and codes are provided to process PacBio long-read RNA-seq data and generate isoform-level master table.
The SLURM settings in all bash scripts are cluster-specific. Adjust memory, CPUs, job array size, and email before submitting on your system.

The pipeline executes the following steps in order:

1. **align_pacbio.sh**  
   - Aligns PacBio reads to the reference genome using minimap2.

2. **before_collapse.sh**  
   - Prepares GTFs and PSL files for downstream collapse steps.

3. **collapse_pacbio.sh**  
   - Collapses isoforms using custom Python scripts to produce a non-redundant isoform set.

4. **sqanti_IR.sh**  
   - Runs SQANTI3 for isoform QC and filters intron retention events.

5. **cpat_nmd.sh**  
   - Runs CPAT to analyze coding potential and filters nonsense-mediated decay (NMD) transcripts.

6. **generate_master_table.R**  
   - Generates a comprehensive master table of all isoforms and genes, including:
     - isoform categories (protein coding, NMD, noORF, RI)
     - canonical/annotated transcript annotation
     - entropy/perplexity
     - TPM expression levels
     - exon counts
     - splicing event integration
     - ORF-level summaries

## How to Run

```bash/run_pacbio_pipeline.sh``` directs to all bash and R scripts needed to run the entire pipeline. For each file, edit user-defined variables to run on your data.
