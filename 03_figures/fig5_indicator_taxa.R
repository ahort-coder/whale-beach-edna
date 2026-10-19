# ============================================================
# Regenerated Publication Figures with Improved Labels
# Includes Pseudoalteromonas and descriptive annotations
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

tax <- as.data.frame(tax_table(ps))

pretty_names <- c("Peat" = "Peat",
                  "Noreaster_Overwash" = "Nor'easter Overwash",
                  "Hurricane_Overwash" = "Hurricane Overwash")

sediment_colors <- c("Peat" = "#8B6914",
                     "Nor'easter Overwash" = "#2166AC",
                     "Hurricane Overwash" = "#B2182B")

output_dir <- "~/Desktop/DESeq2_Sediment_Results"

# Helper: create descriptive labels for uncultured lineages
make_descriptive_label <- function(kingdom, phylum, class, order, family, genus) {
  # If genus is available, use it
  if (!is.na(genus)) return(genus)
  # Otherwise build informative label with taxonomic context
  if (!is.na(family)) {
    ctx <- ifelse(!is.na(phylum), paste0(" (", phylum, ")"), "")
    return(paste0(family, ctx))
  }
  if (!is.na(order)) {
    ctx <- ifelse(!is.na(phylum), paste0(" (", phylum, ")"), "")
    return(paste0(order, ctx))
  }
  if (!is.na(class)) {
    ctx <- ifelse(!is.na(phylum), paste0(" (", phylum, ")"), "")
    return(paste0(class, ctx))
  }
  if (!is.na(phylum)) return(phylum)
  return("Unclassified")
}

tax$descriptive_label <- mapply(make_descriptive_label,
                                 tax$Kingdom, tax$Phylum, tax$Class,
                                 tax$Order, tax$Family, tax$Genus)

# ============================================================
# 2. LOAD INDICATOR RESULTS
# ============================================================
indval_df <- read.csv(file.path(output_dir, "indicator_species_results.csv"),
                       stringsAsFactors = FALSE)

# Add descriptive labels from tax table
tax_lookup <- setNames(tax$descriptive_label, rownames(tax))
indval_df$desc_label <- tax_lookup[indval_df$ASV]
indval_df$desc_label[is.na(indval_df$desc_label)] <- "Unclassified"

# Single-group indicators
single_group <- indval_df[!grepl("\\+", indval_df$indicator_group) &
                           !is.na(indval_df$IndVal), ]
single_group$group_pretty <- pretty_names[single_group$indicator_group]
single_group$group_pretty <- factor(single_group$group_pretty,
                                     levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash"))

# ============================================================
# 3. FIGURE 1: INDICATOR SPECIES DOT PLOT (top 15 per group)
# ============================================================
top_indicators <- single_group %>%
  group_by(group_pretty) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 15) %>%
  ungroup()

top_indicators$desc_label <- make.unique(top_indicators$desc_label, sep = " ")
top_indicators$phylum_clean <- ifelse(is.na(top_indicators$Phylum),
                                       "Unclassified", top_indicators$Phylum)

# Order within facets by IndVal
top_indicators <- top_indicators %>%
  arrange(group_pretty, IndVal) %>%
  mutate(desc_label = factor(desc_label, levels = unique(desc_label)))

# Phylum color palette
phyla_present <- sort(unique(top_indicators$phylum_clean))
n_phyla <- length(phyla_present)
phyla_pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
               "#8491B4", "#91D1C2", "#DC9A7D", "#7E6148", "#B09C85",
               "#E8C51D", "#A65628", "#984EA3", "#FF7F00", "#377EB8",
               "#4DAF4A", "#F781BF", "#999999")
phyla_colors <- setNames(phyla_pal[1:n_phyla], phyla_present)

fig1 <- ggplot(top_indicators, aes(x = IndVal, y = desc_label)) +
  geom_segment(aes(x = 0, xend = IndVal, yend = desc_label),
               color = "grey70", linewidth = 0.4) +
  geom_point(aes(size = -log10(p.value), fill = phylum_clean),
             shape = 21, color = "grey30", stroke = 0.4) +
  facet_grid(group_pretty ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_size_continuous(range = c(2.5, 8),
                        name = expression(-log[10](italic(p)))) +
  scale_fill_manual(values = phyla_colors, name = "Phylum") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, NA)) +
  labs(
    x = "Indicator Value (IndVal)",
    y = NULL,
    title = "Indicator Taxa by Sediment Type",
    subtitle = expression(italic("multipatt")~"IndVal.g | 999 permutations |"~italic(p)~"< 0.05 | Top 15 per group")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "grey40", hjust = 0,
                                     margin = margin(b = 10)),
    strip.text.y.left = element_text(size = 11, face = "bold", angle = 0),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.placement  = "outside",
    axis.text.y      = element_text(size = 9),
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
       fig1, width = 11, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig1_indicator_species_dotplot.pdf"),
       fig1, width = 11, height = 11, bg = "white")
cat("Saved Figure 1\n")

# ============================================================
# 4. FIGURE 2: HEATMAP with descriptive labels
# ============================================================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_rel <- as.data.frame(otu_table(ps_rel))
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)
otu_rel <- as.data.frame(otu_rel)

heat_taxa <- single_group %>%
  group_by(indicator_group) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 10) %>%
  ungroup()

heat_mat <- otu_rel[heat_taxa$ASV, , drop = FALSE]
heat_labels <- make.unique(heat_taxa$desc_label, sep = " ")
rownames(heat_mat) <- heat_labels

sample_order <- samp_df[order(samp_df$sediment_type, samp_df$depth_mid), ]
heat_mat <- heat_mat[, rownames(sample_order), drop = FALSE]
heat_mat_log <- log10(heat_mat + 1e-5)

col_annot <- data.frame(
  Sediment = pretty_names[as.character(sample_order$sediment_type)],
  row.names = rownames(sample_order)
)

row_annot <- data.frame(
  Indicator = pretty_names[heat_taxa$indicator_group],
  row.names = heat_labels
)

ann_colors <- list(
  Sediment  = sediment_colors,
  Indicator = sediment_colors
)

png(file.path(output_dir, "Fig2_indicator_heatmap.png"),
    width = 12, height = 9, units = "in", res = 300)
pheatmap(
  as.matrix(heat_mat_log),
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  annotation_col = col_annot,
  annotation_row = row_annot,
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
  annotation_row = row_annot,
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
cat("Saved Figure 2\n")

# ============================================================
# 5. FIGURE 3: SPECIFICITY vs FIDELITY BUBBLE PLOT
# ============================================================
bubble_df <- single_group %>%
  filter(!is.na(specificity) & !is.na(fidelity)) %>%
  group_by(group_pretty) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 15) %>%
  ungroup()

bubble_df$desc_label <- make.unique(bubble_df$desc_label, sep = " ")

fig3 <- ggplot(bubble_df, aes(x = specificity, y = fidelity)) +
  geom_point(aes(size = IndVal, fill = group_pretty),
             shape = 21, color = "grey30", stroke = 0.4, alpha = 0.8) +
  geom_text_repel(aes(label = desc_label),
                  size = 2.6, fontface = "italic",
                  max.overlaps = 25,
                  segment.color = "grey60",
                  segment.size = 0.3,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  force = 2) +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_size_continuous(range = c(2, 9), name = "IndVal",
                        breaks = c(0.3, 0.5, 0.7, 0.9)) +
  scale_x_continuous(limits = c(0.55, 1.05), breaks = seq(0.6, 1, 0.1),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0.4, 1.08), breaks = seq(0.4, 1, 0.1),
                     labels = label_number(accuracy = 0.1)) +
  labs(
    x = "Specificity (A) — How exclusive to the sediment type",
    y = "Fidelity (B) — How consistently present in the sediment type",
    title = "Indicator Taxa: Specificity vs Fidelity",
    subtitle = "Top 15 per sediment type | Bubble size = IndVal | Upper-right = strongest indicators"
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
       fig3, width = 11, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig3_specificity_fidelity_bubble.pdf"),
       fig3, width = 11, height = 8, bg = "white")
cat("Saved Figure 3\n")

# ============================================================
# 6. FIGURE 4: COMBINED DESeq2 BAR CHART with descriptive labels
# ============================================================
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

all_deseq$desc_label <- mapply(make_descriptive_label,
                                all_deseq$Kingdom, all_deseq$Phylum, all_deseq$Class,
                                all_deseq$Order, all_deseq$Family, all_deseq$Genus)
all_deseq$phylum_clean <- ifelse(is.na(all_deseq$Phylum), "Unclassified", all_deseq$Phylum)

top_deseq <- all_deseq %>%
  group_by(comparison) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 10) %>%
  ungroup()

top_deseq$desc_label <- make.unique(top_deseq$desc_label, sep = " ")
top_deseq$comparison <- factor(top_deseq$comparison,
                                levels = c("Hurricane vs Peat",
                                           "Nor'easter vs Peat",
                                           "Hurricane vs Nor'easter"))

top_deseq <- top_deseq %>%
  arrange(comparison, log2FoldChange) %>%
  mutate(desc_label = factor(desc_label, levels = unique(desc_label)))

# Phylum colors for DESeq2 plot
deseq_phyla <- sort(unique(top_deseq$phylum_clean))
n_dp <- length(deseq_phyla)
deseq_phyla_colors <- setNames(phyla_pal[1:n_dp], deseq_phyla)

fig4 <- ggplot(top_deseq, aes(x = log2FoldChange, y = desc_label, fill = phylum_clean)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
  geom_col(color = "white", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1e", padj),
                x = ifelse(log2FoldChange > 0,
                           log2FoldChange + 0.5,
                           log2FoldChange - 0.5)),
            hjust = ifelse(top_deseq$log2FoldChange > 0, 0, 1),
            size = 2.5, color = "grey40") +
  facet_grid(comparison ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = deseq_phyla_colors, name = "Phylum") +
  scale_x_continuous(
    breaks = seq(-30, 30, by = 5),
    labels = abs(seq(-30, 30, by = 5))
  ) +
  labs(
    x = expression("|"~log[2]~"Fold Change|"),
    y = NULL,
    title = "Differentially Abundant Taxa Across Sediment Comparisons",
    subtitle = expression("DESeq2 Wald test |"~italic(p)[adj]~"< 0.05 (BH) | Top 10 by |LFC| per comparison")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.y      = element_text(size = 9),
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
       fig4, width = 12, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig4_DESeq2_combined_barplot.pdf"),
       fig4, width = 12, height = 11, bg = "white")
cat("Saved Figure 4\n")

# ============================================================
# 7. FIGURE 5: BOX PLOTS - Top 4 per group (includes Pseudoalteromonas)
# ============================================================
key_taxa <- single_group %>%
  group_by(indicator_group) %>%
  arrange(desc(IndVal)) %>%
  slice_head(n = 4) %>%
  ungroup()

cat("\nKey indicator taxa for boxplots:\n")
for (i in 1:nrow(key_taxa)) {
  cat("  ", key_taxa$desc_label[i], "->", key_taxa$indicator_group[i],
      "(IndVal:", round(key_taxa$IndVal[i], 3), ")\n")
}

otu_rel_samples <- as.data.frame(otu_table(ps_rel))
if (taxa_are_rows(ps_rel)) otu_rel_samples <- t(otu_rel_samples)
otu_rel_samples <- as.data.frame(otu_rel_samples)
otu_rel_samples$sample <- rownames(otu_rel_samples)

sd_df <- as.data.frame(sample_data(ps))
sed_map <- setNames(as.character(sd_df$sediment_type), rownames(sd_df))
otu_rel_samples$sediment_type <- pretty_names[sed_map[otu_rel_samples$sample]]
otu_rel_samples$sediment_type <- factor(otu_rel_samples$sediment_type,
                                         levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash"))

avail_asvs <- key_taxa$ASV[key_taxa$ASV %in% colnames(otu_rel_samples)]

key_long <- otu_rel_samples %>%
  select(sample, sediment_type, all_of(avail_asvs)) %>%
  pivot_longer(cols = -c(sample, sediment_type),
               names_to = "ASV", values_to = "rel_abund")

key_labels <- setNames(make.unique(
  key_taxa$desc_label[key_taxa$ASV %in% avail_asvs], sep = " "),
  avail_asvs)
key_long$taxon <- key_labels[key_long$ASV]

# Add indicator group for ordering
key_indicator <- setNames(
  pretty_names[key_taxa$indicator_group[key_taxa$ASV %in% avail_asvs]],
  avail_asvs)
key_long$indicator_for <- key_indicator[key_long$ASV]

# Order facets by indicator group
key_long <- key_long %>%
  mutate(indicator_for = factor(indicator_for,
                                 levels = c("Peat", "Nor'easter Overwash", "Hurricane Overwash")))
key_long$taxon <- factor(key_long$taxon,
                          levels = unique(key_long$taxon[order(key_long$indicator_for)]))

fig5 <- ggplot(key_long, aes(x = sediment_type, y = rel_abund, fill = sediment_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5, color = "grey20") +
  facet_wrap(~ taxon, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = sediment_colors, name = "Sediment Type") +
  scale_y_continuous(labels = label_percent(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL,
    y = "Relative Abundance (%)",
    title = "Top Indicator Taxa: Relative Abundance Distribution by Sediment Type",
    subtitle = "Top 4 indicator taxa per sediment type | Points = individual samples"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text       = element_text(size = 8.5, face = "bold.italic"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    axis.text.x      = element_text(size = 8, angle = 35, hjust = 1),
    axis.text.y      = element_text(size = 8),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.6, "lines"),
    plot.margin      = margin(10, 10, 10, 10)
  )

ggsave(file.path(output_dir, "Fig5_indicator_boxplots.png"),
       fig5, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Fig5_indicator_boxplots.pdf"),
       fig5, width = 14, height = 9, bg = "white")
cat("Saved Figure 5\n")

cat("\n========================================\n")
cat("All figures regenerated successfully!\n")
cat("========================================\n")
