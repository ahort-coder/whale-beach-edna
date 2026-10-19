# ============================================================
# DESeq2 Pairwise Comparison of Sediment Types — WB2501
# Sediment Classification:
#   - Nor'easter Overwash : 20–30 cm depth (midpoint)
#   - Hurricane Overwash  : > 45 cm depth (midpoint)
#   - Peat                : everything else
# ============================================================

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. LOAD PHYLOSEQ OBJECT
# ============================================================
ps <- readRDS("~/Desktop/WB2501_phyloseq_collapsed.rds")

cat("Loaded WB2501 phyloseq object\n")
cat("  Samples:", nsamples(ps), "\n")
cat("  Taxa:   ", ntaxa(ps), "\n\n")

# ============================================================
# 2. CLASSIFY SEDIMENT TYPES BY DEPTH
# ============================================================
samp_df <- as.data.frame(sample_data(ps))
samp_df$depth_mid <- (samp_df$depth_start + samp_df$depth_end) / 2

samp_df$sediment_type <- dplyr::case_when(
  samp_df$depth_mid >= 25 & samp_df$depth_mid <= 35 ~ "Noreaster_Overwash",
  samp_df$depth_mid > 45                            ~ "Hurricane_Overwash",
  TRUE                                              ~ "Peat"
)

samp_df$sediment_type <- factor(samp_df$sediment_type,
                                 levels = c("Peat",
                                            "Noreaster_Overwash",
                                            "Hurricane_Overwash"))

sample_data(ps) <- sample_data(samp_df)

cat("Sample classification:\n")
print(samp_df[order(samp_df$depth_mid), c("depth_start", "depth_end", "depth_mid", "sediment_type")])
cat("\nSample counts per sediment type:\n")
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

  keep_samples <- sample_data(ps_obj)[[group_var]] %in% c(group1, group2)
  ps_sub <- prune_samples(keep_samples, ps_obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  sample_data(ps_sub)[[group_var]] <- droplevels(
    sample_data(ps_sub)[[group_var]]
  )

  dds <- phyloseq_to_deseq2(ps_sub, as.formula(paste("~", group_var)))

  geoMeans <- apply(counts(dds), 1, function(x) {
    exp(mean(log(x[x > 0]), na.rm = TRUE))
  })
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds <- DESeq(dds, fitType = "parametric", test = "Wald")

  # Positive LFC = higher in group2
  res <- results(dds,
                 contrast = c(group_var, group2, group1),
                 alpha = alpha)
  res <- res[order(res$padj), ]

  cat("Summary of results:\n")
  summary(res)

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
# 5. ADD TAXONOMY TO RESULTS
# ============================================================
tax <- as.data.frame(tax_table(ps))
tax$ASV <- rownames(tax)

add_taxonomy <- function(df) {
  merged <- merge(df, tax, by.x = "row.names", by.y = "ASV", all.x = TRUE)
  colnames(merged)[1] <- "ASV"
  merged[order(merged$padj), ]
}

# ============================================================
# 6. VOLCANO PLOTS
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
# 7. SAVE OUTPUTS
# ============================================================
output_dir <- "~/Desktop/DESeq2_WB2501_Results"

for (comp in comparisons) {
  # Full results with taxonomy
  full_tax <- add_taxonomy(as.data.frame(comp$results))
  write.csv(full_tax,
            file = file.path(output_dir, paste0(comp$comparison, "_full_results.csv")),
            row.names = FALSE)

  # Significant taxa with taxonomy
  sig_tax <- add_taxonomy(comp$significant)
  write.csv(sig_tax,
            file = file.path(output_dir, paste0(comp$comparison, "_significant_taxa.csv")),
            row.names = FALSE)
}

# Volcano plots
pdf(file.path(output_dir, "volcano_plots.pdf"), width = 8, height = 6)
for (comp in comparisons) print(plot_volcano(comp))
dev.off()

# MA plots
pdf(file.path(output_dir, "MA_plots.pdf"), width = 8, height = 6)
for (comp in comparisons) {
  plotMA(comp$results,
         main = paste("MA Plot:", gsub("_", " ", comp$comparison)),
         alpha = 0.05,
         colSig = "#E64B35",
         colNonSig = "grey60",
         ylim = c(-10, 10))
}
dev.off()

cat("\n", strrep("=", 60), "\n")
cat("All results saved to:", output_dir, "\n")
cat(strrep("=", 60), "\n")
cat("\nFiles saved:\n")
for (comp in comparisons) {
  cat(" -", comp$comparison, "_full_results.csv\n")
  cat(" -", comp$comparison, "_significant_taxa.csv\n")
}
cat(" - volcano_plots.pdf\n")
cat(" - MA_plots.pdf\n")
