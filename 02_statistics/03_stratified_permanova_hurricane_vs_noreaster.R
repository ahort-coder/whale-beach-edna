# =============================================================
# Hurricane vs Nor'easter: are the overwash deposits' microbial
# eDNA communities distinguishable? (supports Figure 3 / PCoAV2.R)
#
# Tests (Bray-Curtis, same distance as the PCoA in PCoAV2.R):
#   1. PERMANOVA (adonis2) per core, hurricane vs nor'easter samples only
#   2. PERMDISP (betadisper) per core -- checks PERMANOVA isn't just a
#      dispersion artefact
#   3. PERMANOVA on all samples of each core (3 groups) + pairwise
#   4. Combined all-core test: stratified permutation (labels shuffled
#      within core only), statistic = sum of per-core pseudo-F
# =============================================================
suppressMessages({library(phyloseq); library(vegan)})
set.seed(42)
NPERM <- 99999

# Group definitions reverse-engineered from Figure 3 (depth_start, cm).
# WB2501 follows PCoAV2.R exactly; the other cores were matched to the
# plotted point positions in Figure3.png.
groups <- list(  # definitions supplied for the revised Figures (depth_start, cm; hurricane = interval above 40 cm)
  WB2501 = list(noreaster = function(d) d >= 20 & d <= 26, hurricane = function(d) d >= 40),
  WB2506 = list(noreaster = function(d) d >= 24 & d <= 27, hurricane = function(d) d >= 40),
  WB2508 = list(noreaster = function(d) d >= 20 & d <= 30, hurricane = function(d) d >= 40),
  WB2509 = list(noreaster = function(d) d >= 20 & d <= 24, hurricane = function(d) d >= 40)
)

pseudoF <- function(D, g) {            # D: dist matrix (full), g: factor
  n <- length(g); a <- nlevels(g)
  D2 <- as.matrix(D)^2
  sst <- sum(D2[upper.tri(D2)]) / n
  ssw <- sum(sapply(levels(g), function(l) { i <- g == l; m <- D2[i, i]; sum(m[upper.tri(m)]) / sum(i) }))
  ((sst - ssw) / (a - 1)) / (ssw / (n - a))
}

res <- list(); store <- list()
for (core in names(groups)) {
  ps <- readRDS(sprintf("~/Desktop/Phyloseq Files/%s_phyloseq_collapsed.rds", core))
  sd <- as.data.frame(sample_data(ps))
  sd$depth_start <- as.numeric(as.character(sd$depth_start))
  sd$type <- "Salt Marsh Peat"
  sd$type[groups[[core]]$hurricane(sd$depth_start)] <- "Hurricane"
  sd$type[groups[[core]]$noreaster(sd$depth_start)] <- "Nor'easter"
  sd$type <- factor(sd$type, levels = c("Hurricane", "Nor'easter", "Salt Marsh Peat"))
  cat("\n=====", core, "=====\n"); print(table(sd$type))

  D_all <- phyloseq::distance(ps, method = "bray")          # as in PCoAV2.R
  keep  <- sd$type %in% c("Hurricane", "Nor'easter")
  D <- as.dist(as.matrix(D_all)[keep, keep])
  g <- droplevels(sd$type[keep])
  nH <- sum(g == "Hurricane"); nN <- sum(g == "Nor'easter")

  # 1. PERMANOVA: F and R2 from adonis2; p-value by EXACT enumeration of every
  #    possible relabelling (group sizes are tiny, so this is exact, not Monte Carlo)
  ad <- adonis2(D ~ g, permutations = 99)
  nall <- choose(length(g), nH)
  Fobs <- pseudoF(D, g)
  stopifnot(abs(Fobs - ad$F[1]) < 1e-8)
  Fperm <- apply(combn(length(g), nH), 2, function(i) { gg <- factor(rep("Nor'easter", length(g)), levels = levels(g)); gg[i] <- "Hurricane"; pseudoF(D, gg) })
  p_exact <- mean(Fperm >= Fobs - 1e-12)
  # 2. PERMDISP
  bd <- betadisper(D, g); bdt <- permutest(bd, permutations = 9999)
  # ANOSIM for good measure
  an <- anosim(D, g, permutations = 9999)
  # 3. all-sample, 3 group
  ad3 <- adonis2(D_all ~ sd$type, permutations = NPERM)

  res[[core]] <- data.frame(Core = core, n_Hurricane = nH, n_Noreaster = nN,
    Pseudo_F = ad$F[1], R2 = ad$R2[1], PERMANOVA_p = p_exact,
    Min_possible_p = 1 / nall,
    PERMDISP_p = bdt$tab$`Pr(>F)`[1], ANOSIM_R = an$statistic, ANOSIM_p = an$signif,
    ThreeGroup_R2 = ad3$R2[1], ThreeGroup_p = ad3$`Pr(>F)`[1])
  store[[core]] <- list(D = D, g = g)
  print(ad); print(bdt$tab)
}
tab <- do.call(rbind, res)

# 4. Combined all-core test (stratified permutation, labels shuffled within each core)
obsF <- sapply(store, function(s) pseudoF(s$D, s$g))
obs  <- sum(obsF)
perm <- replicate(NPERM, sum(sapply(store, function(s) pseudoF(s$D, sample(s$g)))))
p_comb <- (sum(perm >= obs) + 1) / (NPERM + 1)

# Fisher's combined p (secondary; p-values are discrete so slightly conservative)
fisher_X2 <- -2 * sum(log(tab$PERMANOVA_p)); fisher_p <- pchisq(fisher_X2, 2 * nrow(tab), lower.tail = FALSE)

cat("\n\n############ SUMMARY ############\n")
print(format(tab, digits = 3), row.names = FALSE)
cat(sprintf("\nCombined stratified permutation test (sum of per-core pseudo-F = %.2f): p = %.5f\n", obs, p_comb))
cat(sprintf("Fisher's combined p across cores: chi2 = %.2f, df = %d, p = %.5f\n", fisher_X2, 2 * nrow(tab), fisher_p))
write.csv(tab, "~/Desktop/WB Data Analysis /PCoA_dbRDA_Figures/PERMANOVA_Hurricane_vs_Noreaster_newgroups_results.csv", row.names = FALSE)
