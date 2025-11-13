# === run_deseq2.R ===
library(data.table)
library(Matrix)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Define fallback operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# === Load input ===
counts <- fread("../results/matrix/umi.counts.txt")
metadata <- fread("../data/metadata.tsv")

# === Add full barcode names ===
metadata[, Barcode_full := paste0(Batch, ".", Plate, "_", Barcode)]

# === Prepare count matrix ===
genes <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
rownames(count_matrix) <- genes
colnames(count_matrix) <- colnames(counts)[-1]

# === Filter metadata to valid barcodes ===
metadata <- metadata[grepl("Batch[3-6]", Batch)]
valid_barcodes <- intersect(colnames(count_matrix), metadata$Barcode_full)
count_matrix <- count_matrix[, valid_barcodes, drop = FALSE]
metadata <- metadata[match(valid_barcodes, metadata$Barcode_full), ]
rownames(metadata) <- metadata$Barcode_full

# === Create base DESeq2 object ===
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ VTR
)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

# === Output directory ===
dir.create("../results/deseq2", showWarnings = FALSE, recursive = TRUE)

# === Combined comparisons ===
message("Running combined comparisons across all batches...")

# (1) Stimulated: each VTR vs vector
stim_meta <- metadata[Plate == "stimulated"]
vector_cells <- stim_meta[VTR == "vector"]
unique_vtrs <- setdiff(unique(stim_meta$VTR), "vector")

for (vtr in unique_vtrs) {
  vtr_cells <- stim_meta[VTR == vtr]
  if (nrow(vtr_cells) > 0 & nrow(vector_cells) > 0) {
    sel_meta <- rbind(vtr_cells, vector_cells)
    rownames(sel_meta) <- sel_meta$Barcode_full
    dds_sub <- dds[, sel_meta$Barcode_full]
    dds_sub$VTR <- droplevels(factor(sel_meta$VTR))
    design(dds_sub) <- ~ VTR
    dds_sub <- estimateSizeFactors(dds_sub, type = "poscounts")
    dds_sub <- DESeq(dds_sub)
    res <- results(dds_sub, contrast = c("VTR", vtr, "vector"))
    res_dt <- as.data.table(res, keep.rownames = "gene")
    fwrite(res_dt, paste0("../results/deseq2/DE_all_stim_", vtr, "_vs_vector.tsv"), sep = "\t")
  }
}

# (2) Stim vector vs Unstim vector across all batches
vector_all <- metadata[VTR == "vector"]
if (length(unique(vector_all$Plate)) == 2) {
  rownames(vector_all) <- vector_all$Barcode_full
  dds_sub <- dds[, vector_all$Barcode_full]
  dds_sub$Plate <- droplevels(factor(vector_all$Plate))
  design(dds_sub) <- ~ Plate
  dds_sub <- estimateSizeFactors(dds_sub, type = "poscounts")
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("Plate", "stimulated", "unstimulated"))
  res_dt <- as.data.table(res, keep.rownames = "gene")
  fwrite(res_dt, "../results/deseq2/DE_all_stim_vector_vs_unstim_vector.tsv", sep = "\t")
}

# === Per-batch comparisons ===
message("Running per-batch comparisons...")

for (batch in unique(metadata$Batch)) {
  message("Processing ", batch)
  batch_meta <- metadata[Batch == batch]

  # (3) Per-batch: Stim VTR vs vector
  stim_meta <- batch_meta[Plate == "stimulated"]
  vector_cells <- stim_meta[VTR == "vector"]
  unique_vtrs <- setdiff(unique(stim_meta$VTR), "vector")

  for (vtr in unique_vtrs) {
    vtr_cells <- stim_meta[VTR == vtr]
    if (nrow(vtr_cells) > 0 & nrow(vector_cells) > 0) {
      sel_meta <- rbind(vtr_cells, vector_cells)
      rownames(sel_meta) <- sel_meta$Barcode_full
      sel_counts <- count_matrix[, sel_meta$Barcode_full, drop = FALSE]
      dds_sub <- DESeqDataSetFromMatrix(
        countData = sel_counts,
        colData = sel_meta,
        design = ~ VTR
      )
      dds_sub <- estimateSizeFactors(dds_sub, type = "poscounts")
      dds_sub <- DESeq(dds_sub)
      res <- results(dds_sub, contrast = c("VTR", vtr, "vector"))
      res_dt <- as.data.table(res, keep.rownames = "gene")
      out_file <- paste0("../results/deseq2/DE_", batch, "_stim_", vtr, "_vs_vector.tsv")
      fwrite(res_dt, out_file, sep = "\t")
    }
  }

  # (4) Per-batch: stim vector vs unstim vector
  stim_vector <- batch_meta[Plate == "stimulated" & VTR == "vector"]
  unstim_vector <- batch_meta[Plate == "unstimulated" & VTR == "vector"]
  if (nrow(stim_vector) > 0 & nrow(unstim_vector) > 0) {
    sel_meta <- rbind(stim_vector, unstim_vector)
    rownames(sel_meta) <- sel_meta$Barcode_full
    sel_counts <- count_matrix[, sel_meta$Barcode_full, drop = FALSE]
    dds_sub <- DESeqDataSetFromMatrix(
      countData = sel_counts,
      colData = sel_meta,
      design = ~ Plate
    )
    dds_sub <- estimateSizeFactors(dds_sub, type = "poscounts")
    dds_sub <- DESeq(dds_sub)
    res <- results(dds_sub, contrast = c("Plate", "stimulated", "unstimulated"))
    res_dt <- as.data.table(res, keep.rownames = "gene")
    out_file <- paste0("../results/deseq2/DE_", batch, "_stim_vector_vs_unstim_vector.tsv")
    fwrite(res_dt, out_file, sep = "\t")
  }
}

message("All DE comparisons completed successfully.")
