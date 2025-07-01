#!/usr/bin/env Rscript

# ==========================================
# generate_master_table.R
#
# Recreates the full analysis pipeline
# from Stella Park's Rmd notebook for
# PacBio transcript analysis
#
# Usage:
#   - Edit the USER VARIABLES section below
#   - Run:
#       Rscript generate_master_table.R
#
# Outputs:
#   - master_table.tsv
#   - multiple intermediate RDS files
#
# Author: Stella Park (converted to script)
# ==========================================

# ---------------------------
# USER VARIABLES
# ---------------------------

gtf_path <- "/path/to/gencode.v46.basic.annotation.gtf"
msc_gtf_path <- "/path/to/gencode.v46.MScustom.gtf"
quant_tsv_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_quantify.tsv"
after_ir_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_noIR.txt"
no_orf_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_noIR_cpat.no_ORF.txt"
after_nmd_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_noIR_noNMD_noshort.txt"
canonical_map_path <- "/path/to/gencode.v46.MScustom.canonical.txt"
ptc_whole_sequence_path <- "/path/to/ptc_whole_sequence.tsv"
novel_junction_list_path <- "/path/to/NNC_isoforms_250610.txt"
gene_transcript_convert_path <- "/path/to/gencode_v46_MScustomgene_transcript_convert.txt"
tpm_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_quantify_rpm.tsv"
psl_path <- "/path/to/all_samples_sp_collapse_all_chr_no_treatment_full.psl"
suppa_ioe_path <- "/path/to/suppa_encode_collapse_all_strict.ioe"
event_summary_rds_path <- "/path/to/splicing_event_summary_gene_level_2.rds"
proteinatlas_path <- "/path/to/proteinatlas.tsv"
rbp_data_path <- "/path/to/VanNostrand_RBP_list.csv"
housekeeping_data_path <- "/path/to/Housekeeping_GenesHuman.csv"

# Output
master_table_output_path <- "master_table.tsv"

# ---------------------------
# Load Libraries
# ---------------------------

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(plotly)
  library(reshape2)
  library(rtracklayer)
  library(IRanges)
  library(data.table)
  library(purrr)
})

cat("Libraries loaded.\n")

# ---------------------------
# File Checks
# ---------------------------

paths <- c(
  gtf_path, msc_gtf_path, quant_tsv_path, after_ir_path,
  no_orf_path, after_nmd_path, canonical_map_path,
  ptc_whole_sequence_path, novel_junction_list_path,
  gene_transcript_convert_path, tpm_path, psl_path,
  suppa_ioe_path, event_summary_rds_path,
  proteinatlas_path, rbp_data_path, housekeeping_data_path
)

for (p in paths) {
  if (!file.exists(p)) {
    stop(paste("ERROR: Required file does not exist:", p))
  }
}

cat("All input files verified.\n")

# ---------------------------
# GTF Processing
# ---------------------------

cat("Reading GTF file...\n")

gtf <- read_tsv(gtf_path, col_names = FALSE, skip = 5)

gtf <- gtf %>%
  mutate(
    transcript_id = str_extract(X9, "transcript_id \"[^\"]+\""),
    transcript_id = str_replace_all(transcript_id, 'transcript_id|[" ]', ""),
    gene_id = str_extract(X9, "gene_id \"[^\"]+\""),
    gene_id = str_replace_all(gene_id, 'gene_id|[" ]', ""),
    gene_type = str_extract(X9, "gene_type \"[^\"]+\""),
    gene_type = str_replace_all(gene_type, 'gene_type|[" ]', "")
  )

gencode_genes <- gtf %>%
  filter(X3 == "gene") %>%
  select(gene_id, gene_type) %>%
  mutate(gene_id = sub("\\..*", "", gene_id))

protein_coding_genes <- filter(gencode_genes, gene_type == "protein_coding")

gencode_basic_transcripts <- gtf %>%
  filter(X3 == "transcript") %>%
  select(gene_id, transcript_id) %>%
  mutate(
    gene_id = sub("\\..*", "", gene_id),
    transcript_id = sub("\\..*", "", transcript_id)
  )

cat("Finished reading GTF.\n")

# ---------------------------
# MSCUSTOM GTF
# ---------------------------

cat("Reading MSCustom GTF...\n")

msc_gtf <- import(msc_gtf_path) %>% as.data.frame()

msc_gtf_nmd_filtered <- msc_gtf %>% filter(transcript_type != "nonsense_mediated_decay")
msc_gtf_nmd_only <- msc_gtf %>% filter(transcript_type == "nonsense_mediated_decay")
msc_gtf_ptc_only <- msc_gtf %>% filter(transcript_type == "protein_coding")

msc_nmd_filtered_tx <- unique(sub("\\..*", "", msc_gtf_nmd_filtered$transcript_id))
msc_nmd_only_tx <- unique(sub("\\..*", "", msc_gtf_nmd_only$transcript_id))
msc_ptc_only_tx <- unique(sub("\\..*", "", msc_gtf_ptc_only$transcript_id))

cat("Finished reading MSCustom GTF.\n")

# ---------------------------
# Load Quantification Data
# ---------------------------

cat("Loading quantification data...\n")

all_quant <- read.table(quant_tsv_path, row.names = 1, header = TRUE)

all_quant <- all_quant %>%
  mutate(
    transcript = sapply(strsplit(rownames(.), "_"), function(x) sub("\\..*", "", x[[1]])),
    gene = sapply(strsplit(rownames(.), "_"), function(x) sub("\\..*", "", x[[2]]))
  )

all_transcripts <- select(all_quant, gene, transcript)

cat("Loaded quantification data.\n")

# ---------------------------
# Isoform Categories
# ---------------------------

cat("Loading isoform classification files...\n")

after_ir <- read.table(after_ir_path, header = FALSE)
no_orf <- read.table(no_orf_path, header = FALSE)
after_nmd <- read.table(after_nmd_path, header = FALSE)

cat("Assigning isoform categories...\n")

all_transcripts <- all_transcripts %>%
  mutate(isoform_category = case_when(
    !(tolower(rownames(.)) %in% tolower(after_ir$V1)) ~ "RI",
    tolower(rownames(.)) %in% tolower(no_orf$V1) ~ "noORF",
    !(tolower(rownames(.)) %in% tolower(after_nmd$V1)) ~ "NMD",
    TRUE ~ "protein_coding"
  ))

# Annotated NMD override
all_transcripts[all_transcripts$transcript %in% msc_nmd_only_tx, "isoform_category"] <- "NMD"

# Rescue known protein-coding
rescued_rows <- all_transcripts %>%
  filter(isoform_category == "NMD") %>%
  mutate(rowname = rownames(.)) %>%
  rowwise() %>%
  filter(any(str_extract_all(rowname, "ENST[0-9]+")[[1]] %in% msc_ptc_only_tx)) %>%
  pull(rowname)

all_transcripts[rescued_rows, "isoform_category"] <- "protein_coding"

all_transcripts <- all_transcripts %>%
  mutate(isoform_category = case_when(
    !(gene %in% protein_coding_genes$gene_id) & isoform_category != "RI" ~ "noORF",
    TRUE ~ isoform_category
  ))

cat("Isoform categories assigned:\n")
print(all_transcripts %>% count(isoform_category))

# ---------------------------
# Canonical Transcript Mapping
# ---------------------------

cat("Loading canonical transcript map...\n")

canonical_map <- read.table(
  canonical_map_path,
  sep = ":", col.names = c("canonical_transcript", "gene")
)
canonical_map$canonical_transcript <- sub("\\..*", "", canonical_map$canonical_transcript)
canonical_map$gene <- sub("\\..*", "", canonical_map$gene)

multi_enst_df <- all_transcripts %>%
  filter(str_detect(rownames(.), "&")) %>%
  mutate(row_id = rownames(.),
         enst_ids = str_extract_all(row_id, "ENST[0-9]+")) %>%
  left_join(canonical_map, by = c("gene" = "gene"))

choose_transcript <- function(ensts, canonical, pc_list, fallback) {
  if (!is.na(canonical) && canonical %in% ensts) {
    canonical
  } else if (any(ensts %in% pc_list)) {
    ensts[ensts %in% pc_list][1]
  } else {
    fallback
  }
}

multi_enst_df$transcript <- mapply(
  choose_transcript,
  ensts = multi_enst_df$enst_ids,
  canonical = multi_enst_df$canonical_transcript,
  fallback = multi_enst_df$transcript,
  MoreArgs = list(pc_list = msc_ptc_only_tx)
)

rownames(multi_enst_df) <- multi_enst_df$row_id
all_transcripts[rownames(multi_enst_df), "transcript"] <- multi_enst_df$transcript
all_quant[rownames(all_transcripts), "transcript"] <- all_transcripts$transcript

cat("Canonical mapping completed.\n")

# ---------------------------
# Protein Category Annotation
# ---------------------------

cat("Annotating protein categories...\n")

proteinatlas <- read_tsv(proteinatlas_path, col_names = TRUE)
rbp_data <- read.csv(rbp_data_path)
housekeeping_data <- read.table(housekeeping_data_path, sep = ";", header = TRUE)

protein_class <- proteinatlas %>%
  select(c("Ensembl", "Protein class", "Molecular function")) %>%
  separate_longer_delim("Protein class", delim = ", ") %>%
  separate_longer_delim("Molecular function", delim = ", ")

rbps <- rbp_data %>%
  select("GeneID") %>%
  mutate("Category" = "RNA binding proteins")
colnames(rbps) <- c("Ensembl", "Category")

tfs <- protein_class %>%
  filter(`Protein class` == "Transcription factors") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Transcription factors")

enzymes <- protein_class %>%
  filter(`Protein class` == "Enzymes") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Enzymes")

chaperones <- protein_class %>%
  filter(`Molecular function` == "Chaperone") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Chaperones")

transporters <- protein_class %>%
  filter(`Protein class` == "Transporters") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Transporters")

receptors <- protein_class %>%
  filter(`Molecular function` == "Receptor") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Receptors")

chrom_regulators <- protein_class %>%
  filter(`Molecular function` == "Chromatin regulator") %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Chromatin regulators")

structural <- c("Actin capping", "Actin-binding", "Myosin")
actin_myosin <- protein_class %>%
  filter(`Molecular function` %in% structural) %>%
  select("Ensembl") %>% distinct() %>%
  mutate(Category = "Actin binding/capping and Myosin")

housekeeping <- housekeeping_data %>%
  inner_join(gencode_basic_transcripts, by = c("Ensembl" = "transcript_id")) %>%
  select(gene_id) %>%
  mutate(Ensembl = gene_id, Category = "Housekeeping") %>%
  select(-gene_id) %>%
  distinct()

plot_proteins <- bind_rows(
  rbps, tfs, enzymes, chaperones, transporters,
  receptors, chrom_regulators, actin_myosin, housekeeping
)

entropy_gene_all_class <- merge(
  all_transcripts, plot_proteins,
  by.x = "gene", by.y = "Ensembl", all.x = TRUE
)

entropy_gene_all_class$Category <- entropy_gene_all_class$Category %>% replace_na("Other")

collapsed_categories <- entropy_gene_all_class %>%
  select(gene, Category) %>%
  group_by(gene) %>%
  summarise(Category = paste(Category, collapse = ";"), .groups = "drop")

all_transcripts <- all_transcripts %>%
  left_join(collapsed_categories, by = "gene") %>%
  mutate(gene_protein_category = Category) %>%
  select(-Category)

cat("Protein category annotation complete.\n")

# ---------------------------
# Load TPM
# ---------------------------

cat("Loading TPM data...\n")

all_tpm <- read.table(tpm_path, row.names = 1, header = TRUE)

all_tpm <- all_tpm %>%
  mutate(transcript = all_quant[rownames(.),]$transcript,
         gene = all_quant[rownames(.),]$gene)

cat("TPM data loaded.\n")

# ---------------------------
# Calculate Ratios
# ---------------------------

cat("Calculating isoform ratios...\n")

all_quant_ratio <- all_quant %>%
  group_by(gene) %>%
  mutate(across(starts_with("ENCFF"), ~ .x / sum(.x))) %>%
  ungroup()

# Calculate max, min, mean prob
all_encff_mat <- all_quant_ratio %>% select(starts_with("ENCFF")) %>% as.matrix()
all_quant_ratio$max_ratio <- apply(all_encff_mat, 1, max, na.rm = TRUE)
all_quant_ratio$min_ratio <- apply(all_encff_mat, 1, min, na.rm = TRUE)
all_quant_ratio$prob <- rowMeans(all_encff_mat, na.rm = TRUE)

# Add ratio metrics to master table
all_transcripts <- all_transcripts %>%
  left_join(
    select(
      all_quant_ratio,
      transcript,
      isoform_max_ratio = max_ratio,
      isoform_min_ratio = min_ratio,
      isoform_prob = prob
    ),
    by = "transcript"
  )

cat("Ratios calculated.\n")

# ---------------------------
# Calculate Entropy & Perplexity
# ---------------------------

cat("Calculating gene-level perplexity...\n")

entropy_gene_all <- all_quant_ratio %>%
  filter(prob > 0) %>%
  group_by(gene) %>%
  summarise(
    isoform_potential = n(),
    entropy = -sum(prob * log2(prob)),
    perplexity = 2^entropy,
    .groups = "drop"
  )

all_transcripts <- all_transcripts %>%
  left_join(
    select(entropy_gene_all,
           gene,
           gene_potential = isoform_potential,
           gene_perplexity = perplexity),
    by = "gene"
  )

cat("Entropy and perplexity computed.\n")

# ---------------------------
# Exon Counting from PSL
# ---------------------------

cat("Reading PSL file for exon counting...\n")

all_pacbio_psl <- read.table(
  psl_path,
  sep = "\t",
  col.names = c(
    "matches", "misMatches", "repMatches", "nCount",
    "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert",
    "strand", "qName", "qSize", "qStart", "qEnd",
    "tName", "tSize", "tStart", "tEnd",
    "blockCount", "blockSizes", "qStarts", "tStarts"
  )
)

all_pacbio_psl <- all_pacbio_psl %>%
  mutate(transcript = all_quant[qName,]$transcript,
         gene = all_quant[qName,]$gene)

all_transcripts <- all_transcripts %>%
  left_join(
    select(all_pacbio_psl, transcript, isoform_number_of_exons = blockCount),
    by = "transcript"
  )

cat("Exon counts merged into master table.\n")

# ---------------------------
# Effective Isoform Logic
# ---------------------------

cat("Computing effective isoforms...\n")

all_quant_ratio_rank <- all_quant_ratio %>%
  group_by(gene) %>%
  mutate(rank = rank(-prob)) %>%
  ungroup()

rank_perplexity_merged <- merge(
  select(all_quant_ratio_rank, gene, transcript, rank),
  select(entropy_gene_all, gene, perplexity),
  by = "gene"
)

rank_perplexity_merged <- rank_perplexity_merged %>%
  mutate(effective = ifelse(rank <= round(perplexity), TRUE, FALSE))

all_transcripts <- all_transcripts %>%
  left_join(
    select(rank_perplexity_merged,
           transcript,
           isoform_ranking_within_gene = rank,
           isoform_effective = effective),
    by = "transcript"
  )

cat("Effective isoform logic complete.\n")

# ---------------------------
# Begin ORF ID Analysis
# ---------------------------

cat("Reading ORF-level data...\n")

protein <- fread(ptc_whole_sequence_path)

aa_collapsed <- protein %>%
  group_by(gene, sequence) %>%
  mutate(ORF_id = paste(isoform, collapse = ";")) %>%
  ungroup() %>%
  select(gene, transcript = isoform, ORF_id)

cat("ORF-level mappings prepared.\n")

# ---------------------------
# Finish ORF-level analysis
# ---------------------------

cat("Finishing ORF-level analysis...\n")

# Summarize ORF probabilities
ptc_quant_ratio <- all_quant_ratio %>%
  filter(transcript %in% aa_collapsed$transcript)

aa_collapsed_gene <- aa_collapsed %>%
  left_join(
    select(ptc_quant_ratio, transcript, prob),
    by = "transcript"
  ) %>%
  group_by(gene, ORF_id) %>%
  summarise(ORF_prob = sum(prob, na.rm = TRUE), .groups = "drop")

aa_collapsed_gene_perp <- aa_collapsed_gene %>%
  group_by(gene) %>%
  summarise(
    ORF_potential = n(),
    ORF_entropy = -sum(ORF_prob * log2(ORF_prob)),
    ORF_perplexity = 2^ORF_entropy,
    .groups = "drop"
  )

all_transcripts <- all_transcripts %>%
  left_join(
    unique(
      select(
        ptc_quant_ratio, gene,
        ptc_potential = gene_potential,
        ptc_perplexity = gene_perplexity
      )
    ),
    by = "gene"
  ) %>%
  left_join(
    select(aa_collapsed, transcript, ORF_id),
    by = "transcript"
  ) %>%
  left_join(
    select(aa_collapsed_gene, ORF_id, ORF_prob),
    by = "ORF_id"
  ) %>%
  left_join(
    select(aa_collapsed_gene_perp, gene, ORF_potential, ORF_perplexity),
    by = "gene"
  )

cat("ORF-level annotations added.\n")

# ---------------------------
# Effective Tissue Counts
# ---------------------------

cat("Calculating effective tissue counts...\n")

# Effective tissue logic
all_quant_ratio_rank_no_average <- all_quant_ratio %>%
  group_by(gene) %>%
  mutate(across(-c(transcript, max_ratio, min_ratio, prob),
                ~ rank(-.x, na.last = "keep"))) %>%
  ungroup()

perplexity_per_gene_per_sample <- all_quant_ratio %>%
  group_by(gene) %>%
  summarise(
    across(
      .cols = -c(transcript, max_ratio, min_ratio, prob),
      .fns = ~ if (all(is.na(.x))) 0 else 2^(-sum(.x * log2(.x + 1e-10))),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

comparison_df <- left_join(
  all_quant_ratio_rank_no_average,
  perplexity_per_gene_per_sample,
  by = "gene"
)

result_df <- comparison_df %>%
  mutate(across(
    .cols = matches("^ENCFF.*\\.x$"),
    .fns = ~ .x <= round(get(str_replace(cur_column(), "\\.x$", ".y"))),
    .names = "result_{.col}"
  ))

final_df <- result_df %>%
  select(c("gene", "transcript", starts_with("result_")))

true_counts <- final_df %>%
  rowwise() %>%
  mutate(
    count_true = sum(c_across(starts_with("result")), na.rm = TRUE),
    gene_expressed_samples = sum(!is.na(c_across(starts_with("result")))),
    percent_count_true = count_true / gene_expressed_samples * 100
  ) %>%
  ungroup()

all_transcripts <- all_transcripts %>%
  left_join(
    select(true_counts, transcript,
           isoform_effective_tissue_count = count_true,
           gene_expressed_samples,
           isoform_percent_effective_tissue_count = percent_count_true),
    by = "transcript"
  )

# Calculate ratio standard deviation
all_quant_ratio_std <- all_quant_ratio %>%
  rowwise() %>%
  mutate(ratio_sd = ifelse(
    sum(!is.na(c_across(starts_with("ENCFF")))) <= 1,
    0,
    sd(c_across(starts_with("ENCFF")), na.rm = TRUE)
  )) %>%
  ungroup()

all_transcripts <- all_transcripts %>%
  left_join(
    select(all_quant_ratio_std, transcript, isoform_ratio_sd = ratio_sd),
    by = "transcript"
  )

cat("Effective tissue counts complete.\n")

# Save effective tissue count object
saveRDS(true_counts, file = "effective_tissue_count.rds")
saveRDS(perplexity_per_gene_per_sample, file = "perplexity_per_gene_per_sample.rds")

# ---------------------------
# Splicing Event Data
# ---------------------------

cat("Reading splicing event summaries...\n")

event_summary <- readRDS(event_summary_rds_path)

# NOTE: Insert your logic for merging event_summary if desired
# e.g.
# all_transcripts <- all_transcripts %>%
#     left_join(select(event_summary, transcript, splicing_event = event_type),
#               by = "transcript")

cat("Splicing event merging skipped (no explicit merge requested).\n")

# ---------------------------
# Novel Junctions
# ---------------------------

cat("Adding novel junction annotations...\n")

novel_junction_isoform_list <- readLines(novel_junction_list_path)

all_transcripts <- all_transcripts %>%
  mutate(novel_junction = transcript %in% novel_junction_isoform_list)

cat("Novel junction annotations complete.\n")

# ---------------------------
# Gene Names and Writing Master Table
# ---------------------------

cat("Adding gene names and writing master table...\n")

gene_transcript_convert <- read.table(
  gene_transcript_convert_path,
  header = TRUE
)

gene_transcript_convert$trans_id <- sub("\\..*", "", gene_transcript_convert$trans_id)
gene_transcript_convert$gene_id <- sub("\\..*", "", gene_transcript_convert$gene_id)

# Add gene names to TPM and ratios
all_tpm <- all_tpm %>%
  left_join(
    distinct(select(gene_transcript_convert, gene_id, gene_name)),
    by = c("gene" = "gene_id")
  )
all_quant_ratio <- all_quant_ratio %>%
  left_join(
    distinct(select(gene_transcript_convert, gene_id, gene_name)),
    by = c("gene" = "gene_id")
  )

saveRDS(all_tpm, file = "all_tpm.rds")
saveRDS(all_quant_ratio, file = "all_quant_ratio.rds")

# Add gene and transcript names to master table
all_transcripts <- all_transcripts %>%
  left_join(
    select(gene_transcript_convert, trans_id, gene_name, transcript_name),
    by = c("transcript" = "trans_id")
  ) %>%
  left_join(
    distinct(select(gene_transcript_convert, gene_id, gene_name_new = gene_name)),
    by = c("gene" = "gene_id")
  ) %>%
  mutate(
    gene_name = if_else(is.na(gene_name), gene_name_new, gene_name)
  ) %>%
  select(-gene_name_new)

write.table(
  all_transcripts,
  file = master_table_output_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(paste0("Master table written to ", master_table_output_path, "\n"))

