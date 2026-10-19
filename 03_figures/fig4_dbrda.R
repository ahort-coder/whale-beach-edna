# =============================================================
# Matching PCoA and dbRDA figures for the four Whale Beach cores
#  - same sediment-type definitions, colours and labels in both figures
#  - every sample labelled with its depth interval (cm); labels repelled so they stay readable
#  PCoA : Bray-Curtis on square-root transformed, median-normalised counts
#  dbRDA: same distances, constrained by sediment type
# =============================================================
suppressMessages({library(phyloseq); library(vegan); library(ggplot2); library(ggrepel); library(patchwork); library(MASS)})
# square-root transformed, median-normalised counts -> Bray-Curtis (median step is a no-op for Bray-Curtis; kept to match the stated method)
sqrt_bray <- function(ps) {
  m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) m <- t(m)
  vegan::vegdist(sqrt(m * median(rowSums(m)) / rowSums(m)), "bray")
}

set.seed(42)
IN  <- path.expand("~/Desktop/Phyloseq Files")
OUT <- path.expand("~/Desktop/WB Data Analysis /PCoA_dbRDA_Figures")
FONT <- "Arial"

# ---- sediment types (depth_start, cm) ---------------------------------------
# nor'easter: depth_start inside [lo, hi] (inclusive); hurricane: interval above 40 cm (mid-depth > 40); peat: everything else
rules <- list(WB2501 = c(20, 26), WB2506 = c(24, 27), WB2508 = c(20, 30), WB2509 = c(20, 24))
LV   <- c("Noreaster", "Hurricane", "Peat")
LAB  <- c(Noreaster = "Nor’easter overwash", Hurricane = "Hurricane overwash", Peat = "Salt marsh peat")
COL  <- c(Noreaster = "#5B8DB8", Hurricane = "#C46B6B", Peat = "#A0855B")   # blue / red / brown (same in both figures)

classify <- function(core, ps) {
  sd <- as.data.frame(sample_data(ps))
  for (cl in colnames(sd)) if (is.list(sd[[cl]])) sd[[cl]] <- unlist(sd[[cl]])
  sd$depth_start <- as.numeric(as.character(sd$depth_start)); sd$depth_end <- as.numeric(as.character(sd$depth_end))
  sd$mid <- (sd$depth_start + sd$depth_end) / 2
  r <- rules[[core]]
  sd$type <- "Peat"
  sd$type[sd$mid > 40] <- "Hurricane"
  sd$type[sd$depth_start >= r[1] & sd$depth_start <= r[2]] <- "Noreaster"
  sd$type <- factor(sd$type, levels = LV)
  sd$label <- paste0(sd$depth_start, "–", sd$depth_end)
  sd$sample <- rownames(sd)
  sd
}

# ---- ellipses: 80% t-ellipse (as ggplot2::stat_ellipse); two-point groups get an ellipse enclosing both points ----
ellipse_df <- function(x, y, level = 0.80, segments = 60) {
  d <- cbind(x, y); n <- nrow(d); ang <- seq(0, 2 * pi, length.out = segments + 1); circ <- cbind(cos(ang), sin(ang))
  if (n >= 3) {
    v <- tryCatch(MASS::cov.trob(d), error = function(e) NULL)
    if (is.null(v)) v <- list(cov = cov(d), center = colMeans(d))
    shape <- v$cov + diag(1e-10, 2); ctr <- v$center; rad <- sqrt(2 * qf(level, 2, n - 1))
    e <- t(ctr + rad * t(circ %*% chol(shape)))
  } else {
    ctr <- colMeans(d); dv <- d[2, ] - d[1, ]; len <- sqrt(sum(dv^2)); u <- dv / len; w <- c(-u[2], u[1])
    a <- 0.75 * len; b <- 0.40 * a
    e <- t(ctr + a * outer(u, circ[, 1]) + b * outer(w, circ[, 2]))
  }
  data.frame(x = e[, 1], y = e[, 2])
}
ell_all <- function(df, xv, yv, types = levels(df$type)) do.call(rbind, lapply(types, function(t) {
  s <- df[df$type == t, ]; if (nrow(s) < 2) return(NULL); e <- ellipse_df(s[[xv]], s[[yv]]); e$type <- t; e }))

# ---- common look ------------------------------------------------------------
base_theme <- function(grid) theme_bw(base_size = 7, base_family = FONT) +
  theme(panel.border = element_rect(colour = "grey25", linewidth = 0.4),
        panel.grid.major = if (grid) element_line(colour = "grey90", linewidth = 0.25) else element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 7), axis.text = element_text(size = 6, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey25"), axis.ticks.length = unit(1.2, "mm"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 2)),
        plot.tag = element_text(size = 9, face = "bold", family = FONT),
        legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5),
        legend.key.size = unit(4.5, "mm"), legend.margin = margin(0, 0, 0, 2),
        plot.margin = margin(2, 4, 2, 2))
fills <- scale_fill_manual(values = COL, labels = LAB, name = "Sediment type", drop = FALSE)
cols  <- scale_colour_manual(values = COL, labels = LAB, name = "Sediment type", drop = FALSE)
depth_labels <- function(df, xv, yv) geom_text_repel(data = df, aes(x = .data[[xv]], y = .data[[yv]], label = label),
  size = 1.9, family = FONT, colour = "grey15", segment.size = 0.2, segment.colour = "grey45",
  box.padding = 0.28, point.padding = 0.15, min.segment.length = 0, max.overlaps = Inf, force = 3, force_pull = 0.6,
  max.iter = 20000, seed = 7, show.legend = FALSE)

# ---- build per-core plots -----------------------------------------------------
pcoa_plots <- list(); db_plots <- list(); stats <- list(); classes <- list()
for (core in names(rules)) {
  ps <- readRDS(file.path(IN, sprintf("%s_phyloseq_collapsed.rds", core)))
  sd <- classify(core, ps); classes[[core]] <- cbind(core = core, sd[, c("sample", "depth_start", "depth_end", "type")])
  cat("\n==", core, "==\n"); print(table(sd$type))

  ## PCoA
  bray <- sqrt_bray(ps); pc <- ordinate(ps, "PCoA", distance = bray)
  pt <- as.data.frame(pc$vectors[, 1:2]); colnames(pt) <- c("x", "y"); pt <- cbind(pt, sd[rownames(pt), c("type", "label")]); pt <- pt[order(-as.integer(pt$type)), ]
  eig <- pc$values$Relative_eig * 100
  e <- ell_all(pt, "x", "y", types = c("Noreaster", "Hurricane"))   # PCoA: ellipses around the two overwash types only (as in the original figure)
  pcoa_plots[[core]] <- ggplot(pt, aes(x, y)) +
    geom_path(data = e, aes(x, y, colour = type, group = type), linewidth = 0.45, show.legend = FALSE) +
    geom_point(aes(fill = type), shape = 21, size = 2, colour = "grey25", stroke = 0.25) +
    depth_labels(pt, "x", "y") + fills + cols +
    scale_x_continuous(expand = expansion(mult = 0.10)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste(core, "Bray–Curtis PCoA"), x = sprintf("PCoA1 (%.1f%%)", eig[1]), y = sprintf("PCoA2 (%.1f%%)", eig[2])) +
    base_theme(FALSE)

  ## dbRDA (Bray-Curtis on relative abundance, constrained by sediment type)
  d <- sqrt_bray(ps); otu <- data.frame(row.names = sample_names(ps)); env_df <- data.frame(sediment_type = sd$type[match(rownames(otu), rownames(sd))]); rownames(env_df) <- rownames(otu)
  fit <- dbrda(d ~ sediment_type, data = env_df)
  set.seed(42); an <- anova(fit, permutations = 999)
  sc <- as.data.frame(scores(fit, display = "sites", choices = 1:2)); colnames(sc) <- c("x", "y")
  sc <- cbind(sc, sd[rownames(sc), c("type", "label")]); sc <- sc[order(-as.integer(sc$type)), ]
  ev <- eigenvals(fit); pct <- round(ev[1:2] / sum(ev) * 100, 1)
  e2 <- ell_all(sc, "x", "y")
  db_plots[[core]] <- ggplot(sc, aes(x, y)) +
    geom_polygon(data = e2, aes(x, y, fill = type, colour = type, group = type), alpha = 0.12, linewidth = 0.4, show.legend = FALSE) +
    geom_point(aes(fill = type), shape = 21, size = 2, colour = "grey25", stroke = 0.25) +
    depth_labels(sc, "x", "y") + fills + cols +
    annotate("text", x = Inf, y = Inf, label = core, hjust = 1.08, vjust = 1.5, fontface = "bold", size = 3, family = FONT) +
    scale_x_continuous(expand = expansion(mult = 0.10)) + scale_y_continuous(expand = expansion(mult = 0.14)) +
    labs(x = sprintf("CAP1 (%.1f%% of total variance)", pct[1]), y = sprintf("CAP2 (%.1f%% of total variance)", pct[2])) +
    base_theme(FALSE)
  stats[[core]] <- data.frame(core = core, n_noreaster = sum(sd$type == "Noreaster"), n_hurricane = sum(sd$type == "Hurricane"), n_peat = sum(sd$type == "Peat"),
                              CAP1_pct = pct[1], CAP2_pct = pct[2], F = an$F[1], p = an$`Pr(>F)`[1])
}

assemble <- function(plots) wrap_plots(plots[c("WB2501", "WB2506", "WB2508", "WB2509")], ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "right")

W <- 180; H <- 150   # mm, Nature two-column width
for (nm in c("dbRDA")) {   # PCoA figure is built by make_PCoA_original_style.R
  g <- assemble(db_plots)
  ggsave(file.path(OUT, sprintf("Figure_%s.pdf", nm)), g, width = W, height = H, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT, sprintf("Figure_%s.png", nm)), g, width = W, height = H, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
}
write.csv(do.call(rbind, classes), file.path(OUT, "sample_classification.csv"), row.names = FALSE)
st <- do.call(rbind, stats); write.csv(st, file.path(OUT, "dbRDA_summary_all_cores.csv"), row.names = FALSE); print(st)
