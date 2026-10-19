# ============================================================
# DESeq2 Pairwise Comparison of Sediment Types
# Sediment Classification:
#   - Nor'easter Overwash : 25–35 cm depth
#   - Hurricane Overwash  : > 45 cm depth
#   - Peat                : everything else
# ============================================================

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. LOAD PHYLOSEQ OBJECT
# ============================================================
# Update this path to match your file location and format
# Supported loaders: readRDS() for .rds, load() for .RData
ps <- readRDS("~/Desktop/WB2501_phyloseq_collapsed.rds")
# If .RData: load("~/Desktop/your_phyloseq_file.RData")  # then call the object by name

# ============================================================
# 2. CLASSIFY SEDIMENT TYPES BY DEPTH
# ============================================================
# Extract sample data
samp_df <- as.data.frame(sample_data(ps))

# Check your depth column name — update "depth" below if yours differs
# print(colnames(samp_df))  # Uncomment to inspect column names

depth_col <- "depth"  # <-- UPDATE if your depth column is named differently

samp_df$sediment_type <- case_when(
  samp_df[[depth_col]] >= 25 & samp_df[[depth_col]] <= 35 ~ "Noreaster_Overwash",
  samp_df[[depth_col]] > 45                               ~ "Hurricane_Overwash",
  TRUE                                                    ~ "Peat"
)

samp_df$sediment_type <- factor(samp_df$sediment_type,
                                 levels = c("Peat",
                                            "Noreaster_Overwash",
                                            "Hurricane_Overwash"))

# Write classification back to phyloseq object
sample_data(ps) <- sample_data(samp_df)

cat("Sample counts per sediment type:\n")
print(table(samp_df$sediment_type))

# ============================================================
# 3. HELPER FUNCTION: RUN DESEQ2 PAIRWISE COMPARISON
# ============================================================
run_deseq2_comparison <- function(ps_obj, group1, group2,
                                  group_var = "sediment_type",
                                  alpha = 0.05) {

  cat("\n", strrep("=", 60), "\n")
  cat("Comparing:", group1, "vs", group2, "\n")
  cat(strrep("=", 60), "\n")

  # Subset to the two groups
  ps_sub <- subset_samples(ps_obj, get(group_var) %in% c(group1, group2))
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  # Drop unused factor levels
  sample_data(ps_sub)[[group_var]] <- droplevels(
    sample_data(ps_sub)[[group_var]]
  )

  # Convert to DESeq2 dataset
  dds <- phyloseq_to_deseq2(ps_sub, as.formula(paste("~", group_var)))

  # Estimate size factors (use geometric mean method robust to sparse data)
  geoMeans <- apply(counts(dds), 1, function(x) {
    exp(mean(log(x[x > 0]), na.rm = TRUE))
  })
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds <- DESeq(dds, fitType = "parametric", test = "Wald")

  # Extract results (group2 vs group1 — positive LFC = higher in group2)
  res <- results(dds,
                 contrast = c(group_var, group2, group1),
                 alpha = alpha)
  res <- res[order(res$padj), ]

  cat("Summary of results:\n")
  summary(res)

  # Significant taxa
  sig <- as.data.frame(res) %>%
    filter(!is.na(padj), padj < alpha) %>%
    arrange(padj)

  cat("\nSignificant taxa (padj <", alpha, "):", nrow(sig), "\n")

  return(list(dds = dds, results = res, significant = sig,
              comparison = paste0(group2, "_vs_", group1)))
}

# ============================================================
# 4. RUN ALL THREE PAIRWISE COMPARISONS
# ============================================================
comp1 <- run_deseq2_comparison(ps, "Peat",               "Noreaster_Overwash")
comp2 <- run_deseq2_comparison(ps, "Peat",               "Hurricane_Overwash")
comp3 <- run_deseq2_comparison(ps, "Noreaster_Overwash", "Hurricane_Overwash")

comparisons <- list(comp1, comp2, comp3)

# ============================================================
# 5. VISUALIZE: VOLCANO PLOTS
# ============================================================
plot_volcano <- function(comp_obj, alpha = 0.05, lfc_thresh = 1) {
  df <- as.data.frame(comp_obj$results) %>%
    mutate(
      significant = !is.na(padj) & padj < alpha & abs(log2FoldChange) >= lfc_thresh,
      direction   = case_when(
        significant & log2FoldChange > 0 ~ "Up",
        significant & log2FoldChange < 0 ~ "Down",
        TRUE                             ~ "NS"
      )
    )

  ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = direction)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey60")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey40") +
    labs(
      title = paste("Volcano Plot:", gsub("_", " ", comp_obj$comparison)),
      x = "Log2 Fold Change",
      y = "-Log10(p-value)",
      color = "Direction"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))
}

# ============================================================
# 6. VISUALIZE: MA PLOTS
# ============================================================
plot_ma <- function(comp_obj, alpha = 0.05) {
  plotMA(comp_obj$results,
         main = paste("MA Plot:", gsub("_", " ", comp_obj$comparison)),
         alpha = alpha,
         colSig = "#E64B35",
         colNonSig = "grey60",
         ylim = c(-5, 5))
}

# ============================================================
# 7. SAVE OUTPUTS
# ============================================================
output_dir <- "~/Desktop/DESeq2_Sediment_Results"
dir.create(output_dir, showWarnings = FALSE)

# Save significant results tables as CSVs
for (comp in comparisons) {
  write.csv(
    as.data.frame(comp$results),
    file = file.path(output_dir, paste0(comp$comparison, "_full_results.csv")),
    row.names = TRUE
  )
  write.csv(
    comp$significant,
    file = file.path(output_dir, paste0(comp$comparison, "_significant_taxa.csv")),
    row.names = TRUE
  )
}

# Save volcano plots as PDFs
pdf(file.path(output_dir, "volcano_plots.pdf"), width = 8, height = 6)
for (comp in comparisons) print(plot_volcano(comp))
dev.off()

# Save MA plots as PDFs
pdf(file.path(output_dir, "MA_plots.pdf"), width = 8, height = 6)
for (comp in comparisons) plot_ma(comp)
dev.off()

cat("\n", strrep("=", 60), "\n")
cat("All results saved to:", output_dir, "\n")
cat(strrep("=", 60), "\n")