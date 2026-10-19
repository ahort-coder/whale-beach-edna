# ============================================================
# Indicator Species Analysis & Publication-Quality Figures
# Sediment Types: Peat, Nor'easter Overwash, Hurricane Overwash
# ============================================================

library(phyloseq)
library(indicspecies)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(viridis)
library(gridExtra)
library(patchwork)
library(ggrepel)
library(reshape2)

# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================
ps <- readRDS("~/Desktop/WB2501_phyloseq_collapsed.rds")

samp_df <- as.data.frame(sample_data(ps))
samp_df$depth_mid <- (samp_df$depth_start + samp_df$depth_end) / 2
samp_df$sediment_type <- dplyr::case_when(
  samp_df$depth_mid >= 25 & samp_df$depth_mid <= 35 ~ "Noreaster_Overwash",
  samp_df$depth_mid > 45                            ~ "Hurricane_Overwash",
  TRUE                                              ~ "Peat"
)
samp_df$sediment_type <- factor(samp_df$sediment_type,
                                 levels = c("Peat", "Noreaster_Overwash", "Hurricane_Overwash"))
sample_data(ps) <- sample_data(samp_df)

# Get OTU table and taxonomy
otu <- as.data.frame(otu_table(ps))
if (!taxa_are_rows(ps)) otu <- t(otu)
otu <- as.data.frame(otu)
tax <- as.data.frame(tax_table(ps))

# Create clean taxonomy labels
tax$label <- ifelse(!is.na(tax$Genus), tax$Genus,
             ifelse(!is.na(tax$Family), paste0("f__", tax$Family),
             ifelse(!is.na(tax$Order), paste0("o__", tax$Order),
             ifelse(!is.na(tax$Class), paste0("c__", tax$Class),
             ifelse(!is.na(tax$Phylum), paste0("p__", tax$Phylum),
                    "Unclassified")))))

output_dir <- "~/Desktop/DESeq2_Sediment_Results"

# ============================================================
# 2. INDICATOR SPECIES ANALYSIS (multipatt)
# ============================================================
cat("\n========== INDICATOR SPECIES ANALYSIS ==========\n\n")

# Transpose OTU table: samples as rows, taxa as columns
otu_t <- t(otu)
groups <- as.character(sample_data(ps)$sediment_type)

# Run indicator species analysis
set.seed(42)
indval <- multipatt(otu_t, groups,
                    func = "IndVal.g",
                    control = how(nperm = 999))

cat("Indicator species summary:\n")
summary(indval, indvalcomp = TRUE)

# Extract results
indval_df <- indval$sign
indval_df$ASV <- rownames(indval_df)
indval_df <- indval_df[!is.na(indval_df$p.value) & indval_df$p.value < 0.05, ]

# Add taxonomy
indval_df <- merge(indval_df, tax, by.x = "ASV", by.y = "row.names", all.x = TRUE)

# Determine which group(s) each taxon is indicator for
# Column names s.Peat, s.Noreaster_Overwash, s.Hurricane_Overwash indicate group membership
group_cols <- grep("^s\\.", colnames(indval_df), value = TRUE)
indval_df$indicator_group <- apply(indval_df[, group_cols], 1, function(x) {
  groups_assigned <- gsub("^s\\.", "", group_cols[x == 1])
  paste(groups_assigned, collapse = " + ")
})

# Add IndVal components from the full indval object
# A = specificity, B = fidelity
A_mat <- indval$A
B_mat <- indval$B
indval_df$specificity <- sapply(1:nrow(indval_df), function(i) {
  asv <- indval_df$ASV[i]
  grp_idx <- which(indval_df[i, group_cols] == 1)
  if (length(grp_idx) == 1 && asv %in% rownames(A_mat)) {
    A_mat[asv, grp_idx]
  } else { NA }
})
indval_df$fidelity <- sapply(1:nrow(indval_df), function(i) {
  asv <- indval_df$ASV[i]
  grp_idx <- which(indval_df[i, group_cols] == 1)
  if (length(grp_idx) == 1 && asv %in% rownames(B_mat)) {
    B_mat[asv, grp_idx]
  } else { NA }
})

# Calculate IndVal = sqrt(A * B)
indval_df$IndVal <- sqrt(indval_df$specificity * indval_df$fidelity)

# Save full results
write.csv(indval_df, file.path(output_dir, "indicator_species_results.csv"), row.names = FALSE)
cat("\nSaved indicator species results CSV.\n")

# ============================================================
# 3. FIGURE 1: INDICATOR SPECIES DOT PLOT (publication style)
# ============================================================

# Filter to single-group indicators only and take top 15 per group
single_group <- indval_df[!grepl("\\+", indval_df$indicator_group) &
                           !is.na(indval_df$IndVal), ]

# Pretty group names
pretty_names <- c("Peat" = "Peat",
                  "Noreaster_Overwash" = "Nor'easter Overwash",
                  "Hurricane_Overwash" = "Hurricane Overwash")
single_group$group_pretty <- pretty_names[single_group$indicator_group]
single_group$group_pretty <- factor(single_group$group_pretty,
                                     levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash"))

top_indicators <- single_group %>%
  group_by(group_pretty) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 15) %>%
  ungroup()

# Make labels unique
top_indicators$label <- make.unique(top_indicators$label, sep = " ")

# Phylum for color
top_indicators$phylum_clean <- ifelse(is.na(top_indicators$Phylum),
                                       "Unclassified", top_indicators$Phylum)

# Order within facets by IndVal
top_indicators <- top_indicators %>%
  arrange(group_pretty, IndVal) %>%
  mutate(label = factor(label, levels = unique(label)))

# Color palette for phyla
phyla_present <- sort(unique(top_indicators$phylum_clean))
n_phyla <- length(phyla_present)
phyla_colors <- setNames(
  c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC9A7D", "#7E6148", "#B09C85",
    "#E8C51D", "#A65628", "#984EA3", "#FF7F00", "#377EB8",
    "#4DAF4A", "#F781BF", "#999999")[1:n_phyla],
  phyla_present
)

fig1 <- ggplot(top_indicators, aes(x = IndVal, y = label)) +
  geom_segment(aes(x = 0, xend = IndVal, yend = label),
               color = "grey70", linewidth = 0.4) +
  geom_point(aes(size = -log10(p.value), fill = phylum_clean),
             shape = 21, color = "grey30", stroke = 0.4) +
  facet_grid(group_pretty ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_size_continuous(range = c(2, 7),
                        name = expression(-log[10](italic(p)))) +
  scale_fill_manual(values = phyla_colors, name = "Phylum") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, NA)) +
  labs(
    x = "Indicator Value (IndVal)",
    y = NULL,
    title = "Indicator Taxa by Sediment Type",
    subtitle = expression(italic("multipatt")~"analysis, IndVal.g, 999 permutations, "~italic(p)~"< 0.05")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "grey40", hjust = 0,
                                     margin = margin(b = 10)),
    strip.text.y.left = element_text(size = 10, face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.placement  = "outside",
    axis.text.y      = element_text(size = 9, face = "italic"),
    axis.text.x      = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1))

ggsave(file.path(output_dir, "Fig1_indicator_species_dotplot.png"),
       fig1, width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig1_indicator_species_dotplot.pdf"),
       fig1, width = 10, height = 10, bg = "white")
cat("Saved Figure 1: Indicator species dot plot\n")

# ============================================================
# 4. FIGURE 2: HEATMAP OF INDICATOR TAXA (publication style)
# ============================================================

# Get top 10 indicator taxa per group for heatmap
heat_taxa <- single_group %>%
  group_by(indicator_group) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Extract relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_rel <- as.data.frame(otu_table(ps_rel))
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)
otu_rel <- as.data.frame(otu_rel)

# Subset to indicator taxa
heat_mat <- otu_rel[heat_taxa$ASV, , drop = FALSE]
rownames(heat_mat) <- make.unique(heat_taxa$label, sep = " ")

# Order samples by sediment type
sample_order <- samp_df[order(samp_df$sediment_type, samp_df$depth_mid), ]
heat_mat <- heat_mat[, rownames(sample_order), drop = FALSE]

# Log-transform for better visualization (add pseudocount)
heat_mat_log <- log10(heat_mat + 1e-5)

# Annotation for columns (samples)
col_annot <- data.frame(
  Sediment = pretty_names[as.character(sample_order$sediment_type)],
  row.names = rownames(sample_order)
)

# Annotation for rows (taxa)
row_annot <- data.frame(
  Indicator = pretty_names[heat_taxa$indicator_group],
  Phylum = ifelse(is.na(heat_taxa$Phylum), "Unclassified", heat_taxa$Phylum),
  row.names = make.unique(heat_taxa$label, sep = " ")
)

# Color schemes
sediment_colors <- c("Peat" = "#8B6914",
                     "Nor'easter Overwash" = "#2166AC",
                     "Hurricane Overwash" = "#B2182B")

ann_colors <- list(
  Sediment  = sediment_colors,
  Indicator = sediment_colors
)

# Generate heatmap
png(file.path(output_dir, "Fig2_indicator_heatmap.png"),
    width = 12, height = 9, units = "in", res = 300)
pheatmap(
  as.matrix(heat_mat_log),
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  annotation_col = col_annot,
  annotation_row = row_annot[, "Indicator", drop = FALSE],
  annotation_colors = ann_colors,
  gaps_col = cumsum(table(col_annot$Sediment)),
  fontsize_row = 9,
  fontsize_col = 7,
  fontsize = 10,
  main = "Indicator Taxa Relative Abundance Across Sediment Types",
  angle_col = 45,
  border_color = NA,
  cellwidth = 18,
  cellheight = 16
)
dev.off()

pdf(file.path(output_dir, "Fig2_indicator_heatmap.pdf"),
    width = 12, height = 9)
pheatmap(
  as.matrix(heat_mat_log),
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  annotation_col = col_annot,
  annotation_row = row_annot[, "Indicator", drop = FALSE],
  annotation_colors = ann_colors,
  gaps_col = cumsum(table(col_annot$Sediment)),
  fontsize_row = 9,
  fontsize_col = 7,
  fontsize = 10,
  main = "Indicator Taxa Relative Abundance Across Sediment Types",
  angle_col = 45,
  border_color = NA,
  cellwidth = 18,
  cellheight = 16
)
dev.off()
cat("Saved Figure 2: Indicator taxa heatmap\n")

# ============================================================
# 5. FIGURE 3: BUBBLE PLOT – Specificity vs Fidelity
# ============================================================

bubble_df <- single_group %>%
  filter(!is.na(specificity) & !is.na(fidelity)) %>%
  group_by(group_pretty) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 15) %>%
  ungroup()

bubble_df$label <- make.unique(bubble_df$label, sep = " ")

fig3 <- ggplot(bubble_df, aes(x = specificity, y = fidelity)) +
  geom_point(aes(size = IndVal, fill = group_pretty),
             shape = 21, color = "grey30", stroke = 0.4, alpha = 0.8) +
  geom_text_repel(aes(label = label),
                  size = 2.8, fontface = "italic",
                  max.overlaps = 20,
                  segment.color = "grey60",
                  segment.size = 0.3,
                  box.padding = 0.4) +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_size_continuous(range = c(2, 9), name = "IndVal",
                        breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  labs(
    x = "Specificity (A)",
    y = "Fidelity (B)",
    title = "Indicator Taxa: Specificity vs Fidelity",
    subtitle = "Bubble size proportional to IndVal; top 15 per sediment type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave(file.path(output_dir, "Fig3_specificity_fidelity_bubble.png"),
       fig3, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig3_specificity_fidelity_bubble.pdf"),
       fig3, width = 10, height = 8, bg = "white")
cat("Saved Figure 3: Specificity vs Fidelity bubble plot\n")

# ============================================================
# 6. FIGURE 4: DESeq2 COMBINED LOG2FC BAR CHART (publication)
# ============================================================

# Read all three significant taxa files with taxonomy
files <- list(
  list(file = "Hurricane_Overwash_vs_Peat_significant_taxa_with_taxonomy.csv",
       comparison = "Hurricane vs Peat",
       group_pos = "Hurricane Overwash", group_neg = "Peat"),
  list(file = "Noreaster_Overwash_vs_Peat_significant_taxa_with_taxonomy.csv",
       comparison = "Nor'easter vs Peat",
       group_pos = "Nor'easter Overwash", group_neg = "Peat"),
  list(file = "Hurricane_Overwash_vs_Noreaster_Overwash_significant_taxa_with_taxonomy.csv",
       comparison = "Hurricane vs Nor'easter",
       group_pos = "Hurricane Overwash", group_neg = "Nor'easter Overwash")
)

all_deseq <- do.call(rbind, lapply(files, function(f) {
  df <- read.csv(file.path(output_dir, f$file), stringsAsFactors = FALSE)
  df$comparison <- f$comparison
  df$group_pos <- f$group_pos
  df$group_neg <- f$group_neg
  df
}))

# Create labels
all_deseq$label <- ifelse(!is.na(all_deseq$Genus), all_deseq$Genus,
                   ifelse(!is.na(all_deseq$Family), paste0("f__", all_deseq$Family),
                   ifelse(!is.na(all_deseq$Order), paste0("o__", all_deseq$Order),
                   ifelse(!is.na(all_deseq$Class), paste0("c__", all_deseq$Class),
                   ifelse(!is.na(all_deseq$Phylum), paste0("p__", all_deseq$Phylum),
                          "Unclassified")))))

all_deseq$phylum_clean <- ifelse(is.na(all_deseq$Phylum), "Unclassified", all_deseq$Phylum)

# Top 10 per comparison by |LFC|
top_deseq <- all_deseq %>%
  group_by(comparison) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 10) %>%
  ungroup()

top_deseq$label <- make.unique(top_deseq$label, sep = " ")
top_deseq$direction <- ifelse(top_deseq$log2FoldChange > 0, "Enriched", "Depleted")
top_deseq$comparison <- factor(top_deseq$comparison,
                                levels = c("Hurricane vs Peat",
                                           "Nor'easter vs Peat",
                                           "Hurricane vs Nor'easter"))

# Order labels within facets
top_deseq <- top_deseq %>%
  arrange(comparison, log2FoldChange) %>%
  mutate(label = factor(label, levels = unique(label)))

fig4 <- ggplot(top_deseq, aes(x = log2FoldChange, y = label, fill = phylum_clean)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
  geom_col(color = "white", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1e", padj),
                x = ifelse(log2FoldChange > 0,
                           log2FoldChange + 0.5,
                           log2FoldChange - 0.5)),
            hjust = ifelse(top_deseq$log2FoldChange > 0, 0, 1),
            size = 2.5, color = "grey40") +
  facet_grid(comparison ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = phyla_colors, name = "Phylum") +
  scale_x_continuous(
    breaks = seq(-30, 30, by = 5),
    labels = abs(seq(-30, 30, by = 5))
  ) +
  labs(
    x = expression("|"~log[2]~"Fold Change|"),
    y = NULL,
    title = "Differentially Abundant Taxa Across Sediment Comparisons",
    subtitle = expression("DESeq2 Wald test,"~italic(p)[adj]~"< 0.05 (BH correction); top 10 by |LFC| per comparison")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.y      = element_text(size = 9, face = "italic"),
    axis.text.x      = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 8, face = "italic"),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(output_dir, "Fig4_DESeq2_combined_barplot.png"),
       fig4, width = 11, height = 10, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig4_DESeq2_combined_barplot.pdf"),
       fig4, width = 11, height = 10, bg = "white")
cat("Saved Figure 4: Combined DESeq2 bar plot\n")

# ============================================================
# 7. FIGURE 5: RELATIVE ABUNDANCE BOX PLOTS OF KEY INDICATORS
# ============================================================

# Pick top 3 indicator taxa per sediment type
key_taxa <- single_group %>%
  group_by(indicator_group) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 3) %>%
  ungroup()

# Get relative abundances for these taxa across all samples
ps_rel_df <- as.data.frame(otu_table(ps_rel))
if (taxa_are_rows(ps_rel)) ps_rel_df <- t(ps_rel_df)
ps_rel_df <- as.data.frame(ps_rel_df)
ps_rel_df$sample <- rownames(ps_rel_df)

# Map sediment types using the updated sample data
sd_df <- as.data.frame(sample_data(ps))
sed_map <- setNames(as.character(sd_df$sediment_type), rownames(sd_df))
ps_rel_df$sediment_type <- pretty_names[sed_map[ps_rel_df$sample]]
ps_rel_df$sediment_type <- factor(ps_rel_df$sediment_type,
                                   levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash"))

# Melt to long format for key taxa
key_long <- ps_rel_df %>%
  select(sample, sediment_type, all_of(key_taxa$ASV)) %>%
  pivot_longer(cols = -c(sample, sediment_type),
               names_to = "ASV", values_to = "rel_abund")

# Add labels
key_labels <- setNames(make.unique(key_taxa$label, sep = " "), key_taxa$ASV)
key_long$taxon <- key_labels[key_long$ASV]

# Add which group it indicates
key_indicator <- setNames(pretty_names[key_taxa$indicator_group], key_taxa$ASV)
key_long$indicator_for <- key_indicator[key_long$ASV]

fig5 <- ggplot(key_long, aes(x = sediment_type, y = rel_abund, fill = sediment_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey20") +
  facet_wrap(~ taxon, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_y_continuous(labels = label_percent(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL,
    y = "Relative Abundance",
    title = "Top Indicator Taxa: Relative Abundance by Sediment Type",
    subtitle = "Top 3 indicator taxa per sediment type (IndVal analysis)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text       = element_text(size = 9, face = "bold.italic"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.x      = element_text(size = 8, angle = 35, hjust = 1),
    axis.text.y      = element_text(size = 8),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.6, "lines"),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(file.path(output_dir, "Fig5_indicator_boxplots.png"),
       fig5, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig5_indicator_boxplots.pdf"),
       fig5, width = 12, height = 8, bg = "white")
cat("Saved Figure 5: Indicator taxa box plots\n")

# ============================================================
# SUMMARY
# ============================================================
cat("\n========================================\n")
cat("All outputs saved to:", output_dir, "\n")
cat("========================================\n")
cat("\nFiles generated:\n")
cat("  - indicator_species_results.csv\n")
cat("  - Fig1_indicator_species_dotplot (.png/.pdf)\n")
cat("  - Fig2_indicator_heatmap (.png/.pdf)\n")
cat("  - Fig3_specificity_fidelity_bubble (.png/.pdf)\n")
cat("  - Fig4_DESeq2_combined_barplot (.png/.pdf)\n")
cat("  - Fig5_indicator_boxplots (.png/.pdf)\n")
