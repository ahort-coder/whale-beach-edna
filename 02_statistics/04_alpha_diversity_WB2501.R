# ============================================================
# Alpha Diversity Analysis — WB2501
# Shannon (H') and Simpson's (D) Diversity Indices
# ============================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

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
  samp_df$depth_mid >= 20 & samp_df$depth_mid <= 30 ~ "Noreaster_Overwash",
  samp_df$depth_mid > 40                            ~ "Hurricane_Overwash",
  TRUE                                              ~ "Peat"
)

samp_df$sediment_type <- factor(samp_df$sediment_type,
                                 levels = c("Peat",
                                            "Noreaster_Overwash",
                                            "Hurricane_Overwash"))
sample_data(ps) <- sample_data(samp_df)

cat("Sample classification:\n")
print(samp_df[, c("depth_start", "depth_end", "depth_mid", "sediment_type")])
cat("\nSample counts per sediment type:\n")
print(table(samp_df$sediment_type))

# ============================================================
# 3. CALCULATE DIVERSITY INDICES
# ============================================================
# Get OTU table (samples as rows)
otu_mat <- as.data.frame(otu_table(ps))
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
otu_mat <- as.data.frame(otu_mat)

# Shannon H'
shannon <- diversity(otu_mat, index = "shannon")

# Simpson's D (1 - dominance)
simpson <- diversity(otu_mat, index = "simpson")

# Observed richness
observed <- rowSums(otu_mat > 0)

# Pielou's evenness (J' = H' / ln(S))
evenness <- shannon / log(observed)

# ============================================================
# 4. BUILD SUMMARY TABLE
# ============================================================
div_df <- data.frame(
  Sample         = rownames(otu_mat),
  Depth_Start    = samp_df$depth_start[match(rownames(otu_mat), rownames(samp_df))],
  Depth_End      = samp_df$depth_end[match(rownames(otu_mat), rownames(samp_df))],
  Depth_Mid      = samp_df$depth_mid[match(rownames(otu_mat), rownames(samp_df))],
  Sediment_Type  = as.character(samp_df$sediment_type[match(rownames(otu_mat), rownames(samp_df))]),
  Observed_ASVs  = observed,
  Shannon_H      = round(shannon, 4),
  Simpson_D      = round(simpson, 6),
  Pielou_J       = round(evenness, 4),
  stringsAsFactors = FALSE
)

div_df <- div_df[order(div_df$Depth_Mid), ]

output_dir <- "~/Desktop/Diversity_WB2501_Results"

write.csv(div_df, file.path(output_dir, "diversity_indices_WB2501.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("DIVERSITY SUMMARY — WB2501\n")
cat("========================================\n\n")
print(div_df, row.names = FALSE)

cat("\n\nSummary statistics by sediment type:\n")
cat("------------------------------------\n")
summary_stats <- div_df %>%
  group_by(Sediment_Type) %>%
  summarise(
    n = n(),
    Shannon_mean = round(mean(Shannon_H), 4),
    Shannon_sd   = round(sd(Shannon_H), 4),
    Shannon_min  = round(min(Shannon_H), 4),
    Shannon_max  = round(max(Shannon_H), 4),
    Simpson_mean = round(mean(Simpson_D), 6),
    Simpson_sd   = round(sd(Simpson_D), 6),
    Simpson_min  = round(min(Simpson_D), 6),
    Simpson_max  = round(max(Simpson_D), 6),
    Evenness_mean = round(mean(Pielou_J), 4),
    Evenness_sd   = round(sd(Pielou_J), 4),
    .groups = "drop"
  )
print(as.data.frame(summary_stats))

write.csv(summary_stats, file.path(output_dir, "diversity_summary_by_sediment_WB2501.csv"), row.names = FALSE)

# ============================================================
# 5. PRETTY LABELS & COLORS
# ============================================================
pretty_names <- c("Peat"               = "Peat",
                  "Noreaster_Overwash"  = "Nor'easter Overwash",
                  "Hurricane_Overwash"  = "Hurricane Overwash")

sediment_colors <- c("Peat"                = "#8B6914",
                     "Nor'easter Overwash" = "#2166AC",
                     "Hurricane Overwash"  = "#B2182B")

div_df$Sediment_Pretty <- pretty_names[as.character(div_df$Sediment_Type)]
div_df$Sediment_Pretty <- factor(div_df$Sediment_Pretty,
                                  levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash"))

# ============================================================
# 6. BOX PLOT — SHANNON H'
# ============================================================
p_shannon <- ggplot(div_df, aes(x = Sediment_Pretty, y = Shannon_H, fill = Sediment_Pretty)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, color = "grey20") +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = NULL,
    y = "Shannon Diversity Index (H')",
    title = "Shannon Diversity (H') by Sediment Type - WB2501",
    subtitle = "Higher values indicate greater species richness and evenness"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    axis.text.x      = element_text(size = 11, face = "bold"),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 11),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin      = margin(15, 15, 10, 10)
  )

ggsave(file.path(output_dir, "Shannon_boxplot_WB2501.png"),
       p_shannon, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Shannon_boxplot_WB2501.pdf"),
       p_shannon, width = 8, height = 6, bg = "white")
cat("\nSaved Shannon boxplot\n")

# ============================================================
# 7. BOX PLOT — SIMPSON'S D
# ============================================================
p_simpson <- ggplot(div_df, aes(x = Sediment_Pretty, y = Simpson_D, fill = Sediment_Pretty)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, color = "grey20") +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    x = NULL,
    y = "Simpson Diversity Index (D)",
    title = "Simpson Diversity (D) by Sediment Type - WB2501",
    subtitle = "Values closer to 1 indicate higher diversity and lower dominance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    axis.text.x      = element_text(size = 11, face = "bold"),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 11),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin      = margin(15, 15, 10, 10)
  )

ggsave(file.path(output_dir, "Simpson_boxplot_WB2501.png"),
       p_simpson, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Simpson_boxplot_WB2501.pdf"),
       p_simpson, width = 8, height = 6, bg = "white")
cat("Saved Simpson boxplot\n")

# ============================================================
# 8. COMBINED PANEL PLOT
# ============================================================
library(patchwork)

combined <- p_shannon + p_simpson +
  plot_annotation(
    title = "Alpha Diversity by Sediment Type - WB2501",
    subtitle = "Shannon (H') and Simpson's (D) diversity indices",
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0.5,
                                    margin = margin(b = 15))
    )
  )

ggsave(file.path(output_dir, "Combined_diversity_boxplots_WB2501.png"),
       combined, width = 14, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Combined_diversity_boxplots_WB2501.pdf"),
       combined, width = 14, height = 6, bg = "white")
cat("Saved combined diversity boxplot\n")

# ============================================================
# 9. DEPTH PROFILE SCATTER PLOTS
# ============================================================
p_depth_shannon <- ggplot(div_df, aes(x = Depth_Mid, y = Shannon_H,
                                       color = Sediment_Pretty)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_line(aes(group = 1), color = "grey60", linewidth = 0.4, alpha = 0.5) +
  scale_color_manual(values = sediment_colors, name = "Sediment Type") +
  scale_x_continuous(breaks = seq(0, 70, by = 5)) +
  labs(
    x = "Depth (cm, midpoint)",
    y = "Shannon H'",
    title = "Shannon Diversity vs Depth - WB2501"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin     = margin(10, 10, 10, 10)
  )

p_depth_simpson <- ggplot(div_df, aes(x = Depth_Mid, y = Simpson_D,
                                       color = Sediment_Pretty)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_line(aes(group = 1), color = "grey60", linewidth = 0.4, alpha = 0.5) +
  scale_color_manual(values = sediment_colors, name = "Sediment Type") +
  scale_x_continuous(breaks = seq(0, 70, by = 5)) +
  labs(
    x = "Depth (cm, midpoint)",
    y = "Simpson D",
    title = "Simpson Diversity vs Depth - WB2501"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin     = margin(10, 10, 10, 10)
  )

depth_combined <- p_depth_shannon + p_depth_simpson +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Diversity_vs_depth_WB2501.png"),
       depth_combined, width = 14, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Diversity_vs_depth_WB2501.pdf"),
       depth_combined, width = 14, height = 6, bg = "white")
cat("Saved diversity vs depth plots\n")

cat("\n========================================\n")
cat("All outputs saved to:", output_dir, "\n")
cat("========================================\n")
