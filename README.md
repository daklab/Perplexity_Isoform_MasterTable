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

Download **environment.yaml** file to figure out all dependencies and avoid manual installation of each package.
To create a Conda environment from the file, use:
```
conda env create -f encode_pacbio.yml
conda activate encode_pacbio
```


### `align_pacbio.sh`
Before running this script, you should confirm that all `fastq` files you wish to process are in a single directory. Then, create a textfile containing all filenames of `fastq` files you wish to process.
```
first.fastq.gz
second.fastq.gz
third.fastq.gz
``` 
You should also confirm that you have a file containing chromosome sizes for all chromosomes. You can generate one using pyfaidx as follows:
```
faidx /path/to/genome/reference.fa -i chromsizes > /path/to/output/chrNameLength.txt
```

In the bash script, insert correct user-defined variables for 1) a path to output directory 2) a path to genome reference 3) genome annotationg gtf 4) a chromosome length file 5) path to the directory where all raw `fastq` files are stored, and 6) textfile containing filenames that we craeted above.

For 6) path to FLAIR helper scripts, insert a path to where all FLAIR helper scripts are stored. You can download from FLAIR github directly, or just download files shared in this github repo.

### `before_collapse.sh`
In the bash script, insert correct user-defined variables for 1) path to input directory where all output files of `align_pacbio.sh` are stored in 2) genome annotation gtf 3) a path to ouput directory 4) a prefix for subset gtf files that will be generated from the input genome annotation gtf. 

Also note that this script assumes the data are human data with chromosomes 1-22, X, and Y. Since this is hard-coded, if you wish to run this script on different species, please edit the numbers in the loop
```
for i in {1..22} X Y; do
    ...
done
```

### `collapse_pacbio.sh`
Before running this script, you should generate a `manifest.tsv` file that contains relevant information for each aligned fastq file. This file should only contain three columns, being (1) sample accession (2) tissue information, and (3) path to fastq file.
```
ENCFF011DHZ prefrontal_cortex /path/to/ENCFF011DHZ.fastq.gz
ENCFF012NOH adrenal_gland /path/to/ENCFF012NOH.fastq.gz
...
```

You should also have two bed files that contain relevant information about acceptable TSS and TES range (one bed file each). In our paper, ...

`bash/run_pacbio_pipeline.sh` directs to all bash and R scripts needed to run the entire pipeline. For each file, edit user-defined variables to run on your data.
