# ============================================================
# Add Taxonomy to DESeq2 Significant Taxa Results
# Uses the phyloseq tax_table to assign taxonomy to ASVs
# ============================================================

library(phyloseq)

# Load phyloseq object
ps <- readRDS("~/Desktop/WB2501_phyloseq_collapsed.rds")

# Extract taxonomy table as a data frame
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

cat("Taxonomy table dimensions:", dim(tax), "\n")
cat("Taxonomy ranks:", colnames(tax), "\n\n")

# Find all significant taxa files
sig_files <- list.files("~/Desktop/DESeq2_Sediment_Results",
                        pattern = "_significant_taxa\\.csv$",
                        full.names = TRUE)

cat("Found", length(sig_files), "significant taxa files:\n")
cat(paste(" ", basename(sig_files), collapse = "\n"), "\n\n")

# Process each file
for (f in sig_files) {
  df <- read.csv(f, row.names = 1)

  # Match ASVs to taxonomy
  matched <- merge(df, tax, by.x = "row.names", by.y = "ASV", all.x = TRUE)
  colnames(matched)[1] <- "ASV"

  # Sort by padj
  matched <- matched[order(matched$padj), ]

  # Count matches
  n_matched <- sum(!is.na(matched$Kingdom) | !is.na(matched[, ncol(df) + 2]))
  cat(basename(f), ":", nrow(matched), "taxa,", n_matched, "matched to taxonomy\n")

  # Save with taxonomy
  out_name <- sub("_significant_taxa\\.csv$", "_significant_taxa_with_taxonomy.csv", f)
  write.csv(matched, file = out_name, row.names = FALSE)
  cat("  -> Saved:", basename(out_name), "\n")
}

cat("\nDone!\n")
