library(data.table)
library(Matrix)

# === Define input/output ===
batches <- c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5", "Batch6")
conditions <- c("stimulated", "unstimulated")
base_dir <- "../data"
output_file <- "../results/matrix/all.umi.counts.txt"

# === Initialize containers ===
combined_mat <- NULL
combined_barcodes <- c()
features_ref <- NULL

# === Process each batch-condition pair ===
for (batch in batches) {
  for (stim in conditions) {
    label <- paste0(batch, ".", stim)
    matrix_dir <- file.path(base_dir, label, "out.Solo.out", "Gene", "raw")

    # Skip if missing
    if (!file.exists(file.path(matrix_dir, "matrix.mtx"))) {
      warning("Missing matrix for ", label)
      next
    }

    # Read matrix, features, and barcodes
    mat <- as.matrix(readMM(file.path(matrix_dir, "matrix.mtx")))
    features <- fread(file.path(matrix_dir, "features.tsv"), header = FALSE)
    barcodes <- fread(file.path(matrix_dir, "barcodes.tsv"), header = FALSE)[[1]]
    barcodes <- paste0(label, "_", barcodes)  # Prefix to ensure uniqueness

    # Set rownames/colnames
    colnames(mat) <- barcodes

    # Set reference features on first iteration
    if (is.null(features_ref)) {
      features_ref <- features
    } else {
      stopifnot(all(features_ref$V1 == features$V1))  # Check feature consistency
    }

    # Combine matrices
    if (is.null(combined_mat)) {
      combined_mat <- mat
    } else {
      combined_mat <- cbind(combined_mat, mat)
    }

    combined_barcodes <- c(combined_barcodes, barcodes)
  }
}

# Finalize matrix
rownames(combined_mat) <- features_ref$V1
combined_df <- as.data.frame(combined_mat)

# Write combined count matrix
fwrite(
  data.table(gene = rownames(combined_df), combined_df),
  file = output_file,
  sep = "\t",
  quote = FALSE
)