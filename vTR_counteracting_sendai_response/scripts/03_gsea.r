#!/usr/bin/env Rscript

library(fgsea)
library(data.table)
library(msigdbr)
library(tidyverse)

# === Load desired Reactome pathway names ===
desired_pathways <- fread("../data/pathways.txt", header = FALSE)[[1]]

# === Load Reactome gene sets via MSigDB (v10+ API-compatible) ===
reactome_sets <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
) %>% as.data.table()

# Filter only the desired pathways
reactome_sets <- reactome_sets[gs_name %in% desired_pathways]

# Check for missing pathway names
missing <- setdiff(desired_pathways, unique(reactome_sets$gs_name))
if (length(missing) > 0) {
  warning("These pathway names were not found in MSigDB: ", paste(missing, collapse = ", "))
}

# Format pathways list
pathways <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)

# === DE result files ===
de_files <- list.files("../results/deseq2", pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# Output directory
gsea_dir <- "../results/gsea/reactome"
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# === Run GSEA on each DE result ===
for (file in de_files) {
  message("Processing ", basename(file))

  res <- fread(file)
  res <- res[!is.na(padj) & !is.na(log2FoldChange)]

  # Identify column with gene symbols
  symbol_col <- intersect(c("symbol", "gene_symbol", "label"), colnames(res))
  if (length(symbol_col) == 0) stop("No gene symbol column found in ", file)
  res[, symbol := get(symbol_col[1])]

  # Deduplicate
  res <- res[!duplicated(symbol)]

  # Rank genes
  ranked <- res$log2FoldChange
  names(ranked) <- res$symbol
  ranked <- sort(ranked, decreasing = TRUE)

  # Filter to genes in desired pathways
  ranked <- ranked[names(ranked) %in% unique(unlist(pathways))]

  # Run fgsea
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranked,
    nperm = 10000,
    minSize = 15,
    maxSize = 500
  ) %>% as.data.table() %>%
    arrange(padj)

  # Save
  out_file <- file.path(gsea_dir, paste0("Reactome_", tools::file_path_sans_ext(basename(file)), ".tsv"))
  fwrite(fgsea_res, out_file, sep = "\t")
}

message("All Reactome GSEA analyses complete.")
