# =============================================================
# PCoA figure rebuilt in the style of the original Figure 3 (WB25PCoA), with
#   - the sediment-type groupings used in the dbRDA figure (see `groups` below)
#   - new colours (nor'easter blue, hurricane red, peat brown = same as the dbRDA figure)
#   - a depth label on every sample
# Bray-Curtis on square-root transformed, median-normalised counts; PCoA via phyloseq::ordinate
# =============================================================
suppressMessages({library(phyloseq); library(vegan); library(ggplot2); library(ggrepel); library(patchwork)})
# square-root transformed, median-normalised counts -> Bray-Curtis (median step is a no-op for Bray-Curtis; kept to match the stated method)
sqrt_bray <- function(ps) {
  m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) m <- t(m)
  vegan::vegdist(sqrt(m * median(rowSums(m)) / rowSums(m)), "bray")
}

IN  <- path.expand("~/Desktop/Phyloseq Files")
OUT <- path.expand("~/Desktop/WB Data Analysis /PCoA_dbRDA_Figures")
FONT <- "Arial"; SERIF <- "Times New Roman"

groups <- list(   # same sediment-type definitions as the dbRDA figure (depth_start, cm); hurricane = intervals from 40 cm down
  WB2501 = list(ne = function(d) d >= 20 & d <= 26, hu = function(d) d >= 40),
  WB2506 = list(ne = function(d) d >= 24 & d <= 27, hu = function(d) d >= 40),
  WB2508 = list(ne = function(d) d >= 20 & d <= 30, hu = function(d) d >= 40),
  WB2509 = list(ne = function(d) d >= 20 & d <= 24, hu = function(d) d >= 40))
# (earlier Figure 3 grouping, kept for reference: WB2501 ne 20-28 | WB2506 ne 20,25,26,30 | WB2508 ne 15,20 | WB2509 ne 20-22; hurricane > 44 in all)
LV  <- c("Hurricane", "Noreaster", "Peat")                       # legend order of the original figure
LAB <- c(Hurricane = "Hurricane", Noreaster = "Nor’easter", Peat = "Salt Marsh Peat")
COL <- c(Noreaster = "#5B8DB8", Hurricane = "#C46B6B", Peat = "#A0855B")
tint <- function(col, f = 0.45) rgb(t((1 - f) * col2rgb(col) + f * 255), maxColorValue = 255)
ELL <- sapply(COL, tint)                                          # lighter ellipse outlines, as in the original

# ellipse that encloses every point of a group (original ellipses were drawn around the whole group)
enclose <- function(x, y, pad = 1.12, min_ratio = 0.35, n = 120) {
  d <- cbind(x, y); ctr <- colMeans(d); ang <- seq(0, 2 * pi, length.out = n); circ <- rbind(cos(ang), sin(ang))
  if (nrow(d) == 2) {
    dv <- d[2, ] - d[1, ]; len <- sqrt(sum(dv^2)); u <- dv / len; w <- c(-u[2], u[1]); a <- 0.75 * len; b <- 0.40 * a
    e <- ctr + a * outer(u, circ[1, ]) + b * outer(w, circ[2, ])
  } else {
    ev <- eigen(cov(d)); l <- ev$values; l[2] <- max(l[2], min_ratio^2 * l[1]); V <- ev$vectors
    z <- sweep(d, 2, ctr) %*% V; r <- max(sqrt(z[, 1]^2 / l[1] + z[, 2]^2 / l[2])) * pad
    e <- ctr + V %*% (r * sqrt(l) * circ)
  }
  data.frame(x = e[1, ], y = e[2, ])
}

th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.3), axis.ticks.length = unit(1.2, "mm"),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6, colour = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 2)),
        plot.tag = element_text(size = 11, family = SERIF), plot.margin = margin(2, 5, 2, 2), legend.position = "none")

plots <- list(); cls <- list()
for (core in names(groups)) {
  ps <- readRDS(file.path(IN, sprintf("%s_phyloseq_collapsed.rds", core)))
  sd <- as.data.frame(sample_data(ps)); ds <- as.numeric(as.character(sd$depth_start)); de <- as.numeric(as.character(sd$depth_end))
  type <- rep("Peat", length(ds)); type[groups[[core]]$hu(ds)] <- "Hurricane"; type[groups[[core]]$ne(ds)] <- "Noreaster"
  sd <- data.frame(sample = rownames(sd), type = factor(type, levels = LV), label = paste0(ds, "–", de), depth_start = ds, depth_end = de, row.names = rownames(sd))
  cls[[core]] <- cbind(core = core, sd[, c("sample", "depth_start", "depth_end", "type")]); cat(core, ":", paste(names(table(sd$type)), table(sd$type), collapse = ", "), "\n")
  pc <- ordinate(ps, "PCoA", distance = sqrt_bray(ps)); eig <- pc$values$Relative_eig * 100
  pt <- cbind(as.data.frame(pc$vectors[, 1:2]), sd[rownames(pc$vectors), c("type", "label")]); colnames(pt)[1:2] <- c("x", "y")
  pt <- pt[order(pt$type != "Peat"), ]   # peat drawn first so overwash points are never hidden
  ell <- do.call(rbind, lapply(c("Noreaster", "Hurricane"), function(t) { s <- pt[pt$type == t, ]; e <- enclose(s$x, s$y); e$type <- t; e }))
  g <- ggplot(pt, aes(x, y)) +
    geom_blank() +
    geom_point(aes(colour = type), size = 2, alpha = 0.95) +
    geom_text_repel(aes(label = label), size = 1.9, family = FONT, colour = "grey15", segment.size = 0.2, segment.colour = "grey45",
                    box.padding = 0.32, point.padding = 0.25, min.segment.length = 0, max.overlaps = Inf, force = 4, force_pull = 0.5,
                    max.iter = 30000, seed = 7, bg.color = "white", bg.r = 0.12) +
    scale_colour_manual(values = COL, breaks = LV, labels = LAB, name = "Sediment Type") +
    scale_x_continuous(expand = expansion(mult = 0.10)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste(core, "Bray–Curtis PCoA"), x = sprintf("PCoA1 (%.1f%%)", eig[1]), y = sprintf("PCoA2 (%.1f%%)", eig[2])) + th
  for (t in c("Noreaster", "Hurricane")) g <- g + geom_path(data = ell[ell$type == t, ], aes(x, y), colour = ELL[[t]], linewidth = 0.55, inherit.aes = FALSE)
  plots[[core]] <- g
}
# legend inside panel d (bottom right, boxed), as in the original
plots$WB2509 <- plots$WB2509 + scale_x_continuous(expand = expansion(mult = c(0.10, 0.34))) + guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "inside", legend.position.inside = c(0.985, 0.985), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5), legend.key.size = unit(3.6, "mm"),
        legend.margin = margin(3, 5, 3, 4), legend.key = element_blank())
fig <- wrap_plots(plots[c("WB2501", "WB2506", "WB2508", "WB2509")], ncol = 2) + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(file.path(OUT, "Figure_PCoA.pdf"), fig, width = 180, height = 150, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_PCoA.png"), fig, width = 180, height = 150, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
write.csv(do.call(rbind, cls), file.path(OUT, "PCoA_sample_classification.csv"), row.names = FALSE)
