#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrepel)
library(biomaRt)

# === Setup directories ===
results_dir <- "../results/deseq2"
volcano_dir <- "../results/volcano_plots"
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)

# === List DE result files ===
files <- list.files(results_dir, pattern = "^DE_.*\\.tsv$", full.names = TRUE)

# === Gather unique Ensembl gene IDs (strip version numbers) ===
all_genes <- unique(unlist(lapply(files, function(f) {
  if (!"gene" %in% names(fread(f, nrows = 0))) {
    message("⚠️ Skipping ", f, " — no 'gene' column found.")
    return(NULL)
  }
  df <- fread(f, select = "gene")
  return(sub("\\..*", "", df$gene))  # remove version suffix
})))

if (length(all_genes) == 0) stop("No gene IDs found. Aborting.")

# === Get HGNC symbols from biomaRt ===
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = ensembl
)
setDT(gene_map)
setnames(gene_map, c("ensembl_gene_id", "hgnc_symbol"), c("gene_base", "gene_symbol"))

# === Loop through each DE file ===
for (file in files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  res <- fread(file)

  if (!all(c("gene", "log2FoldChange", "padj", "pvalue") %in% names(res))) {
    message("Skipping ", file, " — missing required columns.")
    next
  }

  # Add base Ensembl ID
  res[, gene_base := sub("\\..*", "", gene)]

  # Drop conflicting columns if they exist
  cols_to_drop <- intersect(c("gene_symbol", "gene_symbol.x", "gene_symbol.y", "symbol", "label"), names(res))
  if (length(cols_to_drop) > 0) res[, (cols_to_drop) := NULL]

  # Merge with gene symbols
  res <- merge(res, gene_map, by = "gene_base", all.x = TRUE, sort = FALSE)

  # Label = gene symbol if available, else Ensembl ID
  res[, label := ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)]

  # Fix p-values for plotting
  res[, pvalue := pmax(pvalue, 1e-300, na.rm = TRUE)]

  # Set status for volcano
  res[, status := fifelse(padj < 0.05 & log2FoldChange > 1, "up",
                          fifelse(padj < 0.05 & log2FoldChange < -1, "down", "ns"))]

  # Top 40 by significance
  top_genes <- res[!is.na(padj)][order(padj)][1:min(.N, 40)]

  # Plot
  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = status)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c(up = "red", down = "blue", ns = "grey")) +
    geom_text_repel(data = top_genes, aes(label = label), size = 2.5, max.overlaps = 50) +
    theme_classic() +
    labs(title = base_name, x = "log2 Fold Change", y = "-log10(p-value)") +
    theme(legend.position = "none")

  # Save plot
  out_png <- file.path(volcano_dir, paste0(base_name, ".png"))
  ggsave(out_png, plot = p, width = 6, height = 5, units = "in", dpi = 150)
  message("Saved volcano: ", out_png)

  # Clean final table and write it back
  res_out <- res[, .SD, .SDcols = setdiff(names(res), c("gene_base"))]
  cols_final <- c("gene", setdiff(names(res_out), c("gene", "label")), "label")
  setcolorder(res_out, cols_final)
  fwrite(res_out, file, sep = "\t")
}