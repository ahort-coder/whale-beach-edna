# =============================================================
# Do the four cores differ from each other in community composition
# (PERMANOVA) and in dispersion (PERMDISP)? -> justifies analysing cores
# individually.
#
# Cores were collapsed separately, so they must be put on a common
# feature set before pooling. Done at three levels:
#   Phylum / Class : ASVs summed within taxon; unresolved ASVs ("NA")
#                    kept as one "Unassigned" category
#   ASV            : exact-sequence match across cores (union of ASVs,
#                    zero where an ASV is absent from a core)
# Bray-Curtis on relative abundance (cores differ ~3x in reads/sample).
# Sensitivity: rarefy every sample to the same depth before collapsing.
# =============================================================
suppressMessages({library(phyloseq); library(vegan)})
set.seed(42)
NPERM <- 9999
cores <- c("WB2501", "WB2506", "WB2508", "WB2509")

groups <- list(  # final deposit-type definitions (same as the dbRDA figure); depth_start in cm
  WB2501 = list(ne = function(d) d >= 20 & d <= 26, hu = function(d) d >= 40),
  WB2506 = list(ne = function(d) d >= 24 & d <= 27, hu = function(d) d >= 40),
  WB2508 = list(ne = function(d) d >= 20 & d <= 30, hu = function(d) d >= 40),
  WB2509 = list(ne = function(d) d >= 20 & d <= 24, hu = function(d) d >= 40))

ps_list <- lapply(setNames(cores, cores), function(c) readRDS(sprintf("~/Desktop/Phyloseq Files/%s_phyloseq_collapsed.rds", c)))

meta <- do.call(rbind, lapply(cores, function(c) {
  sd <- as.data.frame(sample_data(ps_list[[c]])); d <- as.numeric(as.character(sd$depth_start))
  type <- rep("Peat", length(d)); type[groups[[c]]$hu(d)] <- "Hurricane"; type[groups[[c]]$ne(d)] <- "Noreaster"
  data.frame(sample = paste(c, rownames(sd), sep = "_"), core = c, depth_start = d, type = type,
             reads = sample_sums(ps_list[[c]]), stringsAsFactors = FALSE)
}))
rownames(meta) <- meta$sample
cat("Samples per core x sediment type:\n"); print(table(meta$core, meta$type))

# ---- build pooled sample x feature count matrices ----------------------
get_counts <- function(ps) { m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) m <- t(m); m }  # samples x ASV
pool <- function(mats) {                      # union of features, zero-filled
  feats <- unique(unlist(lapply(mats, colnames)))
  do.call(rbind, lapply(mats, function(m) { out <- matrix(0, nrow(m), length(feats), dimnames = list(rownames(m), feats)); out[, colnames(m)] <- m; out }))
}
asv_mats <- lapply(cores, function(c) { m <- get_counts(ps_list[[c]]); rownames(m) <- paste(c, rownames(m), sep = "_"); m })
names(asv_mats) <- cores

collapse_rank <- function(mats, rank) {
  lapply(cores, function(c) {
    m <- mats[[c]]; tt <- as.data.frame(tax_table(ps_list[[c]]))[colnames(m), rank]
    tt[is.na(tt) | tt == ""] <- "Unassigned"
    out <- rowsum(t(m), group = tt); t(out)
  }) -> L; names(L) <- cores; L
}
rarefy_mat <- function(m, depth) { out <- suppressWarnings(vegan::rrarefy(m, depth)); dimnames(out) <- dimnames(m); out }

pw_store <- list()
run_level <- function(label, mats) {
  X <- pool(mats); X <- X[rownames(meta), ]
  rel <- sqrt(X * median(rowSums(X)) / rowSums(X)); D <- vegdist(rel, "bray")   # median-normalised, sqrt-transformed
  cat("\n\n########", label, "-", ncol(X), "features,", nrow(X), "samples ########\n")
  # (1) core effect, all samples
  a1 <- adonis2(D ~ core, data = meta, permutations = NPERM)
  # (2) PERMDISP by core
  bd <- betadisper(D, meta$core); pd <- permutest(bd, permutations = NPERM)
  # (3) storm deposits only
  k <- meta$type != "Peat"; Dk <- as.dist(as.matrix(D)[k, k]); mk <- droplevels(meta[k, ])
  a2 <- adonis2(Dk ~ core, data = mk, permutations = NPERM)
  bd2 <- betadisper(Dk, mk$core); pd2 <- permutest(bd2, permutations = NPERM)
  # (4) core, sediment type and their interaction (does the storm signal differ among cores?)
  a3 <- adonis2(D ~ core * type, data = meta, permutations = NPERM, by = "terms")
  # (5) is dispersion tied to sequencing depth?
  rho <- suppressWarnings(cor.test(bd$distances, meta$reads, method = "spearman"))
  # (6) pairwise core PERMANOVA, Holm-adjusted
  pr <- combn(cores, 2, simplify = FALSE); pw <- do.call(rbind, lapply(pr, function(p) {
    i <- meta$core %in% p; Dp <- as.dist(as.matrix(D)[i, i]); mp <- meta[i, ]
    a <- adonis2(Dp ~ core, data = mp, permutations = NPERM); b <- permutest(betadisper(Dp, mp$core), permutations = NPERM)
    data.frame(pair = paste(p, collapse = " vs "), F = a$F[1], R2 = a$R2[1], p = a$`Pr(>F)`[1], PERMDISP_p = b$tab$`Pr(>F)`[1])
  })); pw$p_holm <- p.adjust(pw$p, "holm")
  pw_store[[label]] <<- cbind(Level = label, pw)
  cat("\n-- PERMANOVA ~ core (all samples)\n"); print(a1)
  cat("-- PERMDISP by core (all samples)\n"); print(pd$tab); print(round(tapply(bd$distances, meta$core, mean), 3))
  cat("-- Storm deposits only: PERMANOVA ~ core\n"); print(a2); cat("   PERMDISP\n"); print(pd2$tab)
  cat("-- PERMANOVA ~ core * sediment type\n"); print(a3)
  cat(sprintf("-- Spearman, distance-to-centroid vs reads/sample: rho = %.2f, p = %.3f\n", rho$estimate, rho$p.value))
  cat("-- Pairwise cores\n"); print(pw, digits = 3)
  data.frame(Level = label, Features = ncol(X),
    Core_F = a1$F[1], Core_R2 = a1$R2[1], Core_p = a1$`Pr(>F)`[1],
    PERMDISP_F = pd$tab$F[1], PERMDISP_p = pd$tab$`Pr(>F)`[1],
    StormOnly_F = a2$F[1], StormOnly_R2 = a2$R2[1], StormOnly_p = a2$`Pr(>F)`[1], StormOnly_PERMDISP_p = pd2$tab$`Pr(>F)`[1],
    Storm_R2 = a3["type", "R2"], Storm_p = a3["type", "Pr(>F)"], CoreXStorm_R2 = a3["core:type", "R2"], CoreXStorm_p = a3["core:type", "Pr(>F)"],
    Rho_dist_reads = rho$estimate, Rho_p = rho$p.value)
}

out <- list()
for (r in c("Phylum", "Class")) { out[[r]] <- run_level(r, collapse_rank(asv_mats, r)); }
out[["ASV"]] <- run_level("ASV (exact sequence)", asv_mats)

# sensitivity: equal sequencing effort (rarefy to the smallest library, 6040 reads)
dmin <- min(meta$reads); cat("\n\nRarefying to", dmin, "reads per sample\n")
asv_r <- lapply(asv_mats, rarefy_mat, depth = dmin)
for (r in c("Phylum", "Class")) {
  ps_r <- asv_r; mats <- lapply(cores, function(c) { m <- ps_r[[c]]; tt <- as.data.frame(tax_table(ps_list[[c]]))[colnames(m), r]; tt[is.na(tt) | tt == ""] <- "Unassigned"; t(rowsum(t(m), group = tt)) }); names(mats) <- cores
  out[[paste0(r, "_rarefied")]] <- run_level(paste0(r, " (rarefied to ", dmin, ")"), mats)
}
out[["ASV_rarefied"]] <- run_level(paste0("ASV (rarefied to ", dmin, ")"), asv_r)

tab <- do.call(rbind, out); rownames(tab) <- NULL
cat("\n\n############ SUMMARY ############\n"); print(format(tab, digits = 3), row.names = FALSE)

# assignment coverage per core (reads assigned at rank)
cov <- do.call(rbind, lapply(cores, function(c) { tt <- as.data.frame(tax_table(ps_list[[c]])); m <- get_counts(ps_list[[c]])
  data.frame(core = c, sapply(c("Phylum", "Class", "Order", "Family", "Genus"), function(r) round(100 * sum(colSums(m)[!(is.na(tt[colnames(m), r]) | tt[colnames(m), r] == "")]) / sum(m), 1))) }))
cat("\n% of reads assigned at each rank:\n"); print(cov)
outdir <- "~/Desktop/WB Data Analysis /Supplementary_Tables_Update/sqrt_mednorm"
write.csv(tab, file.path(outdir, "S_across_cores_finalgroups.csv"), row.names = FALSE)
write.csv(do.call(rbind, pw_store), file.path(outdir, "S_pairwise_cores.csv"), row.names = FALSE)
