#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

# === Load list of target Reactome pathways ===
target_pathways <- fread("../data/pathways.txt", header = FALSE)$V1

# === Get list of GSEA result files ===
gsea_dir <- "../results/gsea/reactome"
gsea_files <- list.files(gsea_dir, pattern = "^Reactome_.*\\.tsv$", full.names = TRUE)

# === Initialize list to collect NES per file ===
all_results <- list()

for (file in gsea_files) {
  gsea <- fread(file)

  # Filter only pathways of interest
  gsea <- gsea[pathway %in% target_pathways]

  # Replace non-significant with 0, keep NES otherwise
  gsea[, NES_sig := ifelse(padj < 0.05, NES, 0)]

  # Store
  condition_name <- tools::file_path_sans_ext(basename(file)) %>% str_remove("^Reactome_")
  all_results[[condition_name]] <- gsea[, .(pathway, NES_sig)]
  setnames(all_results[[condition_name]], "NES_sig", condition_name)
}

# === Merge all NES columns by pathway ===
summary_table <- reduce(all_results, full_join, by = "pathway")

# === Ensure all target pathways are included ===
summary_table <- merge(data.table(pathway = target_pathways), summary_table, by = "pathway", all.x = TRUE)

# === Replace any remaining NA with 0 ===
summary_table[is.na(summary_table)] <- 0

# === Save the final table ===
fwrite(summary_table, file = "../results/gsea/summary_reactome_nes.tsv", sep = "\t")