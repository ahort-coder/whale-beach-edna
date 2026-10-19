# SQRT / MEDIAN-NORMALISED VERSION. Statistics for the Master Supplementary Document, using the final deposit-type definitions (same as the dbRDA figure).
#   nor'easter: depth_start inside the range below (inclusive); hurricane: intervals from 40 cm down; peat: everything else
#   Bray-Curtis on relative abundance of the collapsed phyloseq objects.
suppressMessages({library(phyloseq); library(vegan)})
set.seed(42)
IN  <- path.expand("~/Desktop/Phyloseq Files"); OUT <- path.expand("~/Desktop/WB Data Analysis /Supplementary_Tables_Update/sqrt_mednorm")
ne_rng <- list(WB2501 = c(20, 26), WB2506 = c(24, 27), WB2508 = c(20, 30), WB2509 = c(20, 24))
TYPES <- c("Hurricane", "Noreaster", "Peat")

load_core <- function(core) {
  ps <- readRDS(file.path(IN, sprintf("%s_phyloseq_collapsed.rds", core)))
  sd <- as.data.frame(sample_data(ps)); ds <- as.numeric(as.character(sd$depth_start)); de <- as.numeric(as.character(sd$depth_end))
  type <- rep("Peat", length(ds)); type[ds >= 40] <- "Hurricane"; r <- ne_rng[[core]]; type[ds >= r[1] & ds <= r[2]] <- "Noreaster"
  m <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) m <- t(m)
  m <- sqrt(m * median(rowSums(m)) / rowSums(m))   # median-normalised, then square-root transformed
  list(D = vegdist(m, "bray"), g = factor(type, levels = TYPES), label = paste0(ds, "-", de), ds = ds, raw = ps)
}
cores <- lapply(setNames(names(ne_rng), names(ne_rng)), load_core)

pseudoF <- function(D2, g) { n <- length(g); a <- nlevels(g); sst <- sum(D2[upper.tri(D2)]) / n
  ssw <- sum(sapply(levels(g), function(l) { i <- g == l; m <- D2[i, i, drop = FALSE]; sum(m[upper.tri(m)]) / sum(i) })); c(F = ((sst - ssw) / (a - 1)) / (ssw / (n - a)), SST = sst, SSW = ssw) }

## 1. definitions table
defs <- do.call(rbind, lapply(names(cores), function(k) { x <- cores[[k]]; do.call(rbind, lapply(TYPES, function(t) { i <- x$g == t
  data.frame(core = k, type = t, n = sum(i), intervals = paste(x$label[i][order(x$ds[i])], collapse = ", ")) })) }))
write.csv(defs, file.path(OUT, "S_deposit_definitions.csv"), row.names = FALSE)

## 2. average between-group Bray-Curtis dissimilarity
pairs <- list(c("Hurricane", "Peat"), c("Noreaster", "Peat"), c("Hurricane", "Noreaster"))
bc <- do.call(rbind, lapply(names(cores), function(k) { x <- cores[[k]]; Dm <- as.matrix(x$D)
  do.call(rbind, lapply(pairs, function(p) data.frame(core = k, comparison = paste(p, collapse = " vs "), dissimilarity = 100 * mean(Dm[x$g == p[1], x$g == p[2]])))) }))
bc$similarity <- 100 - bc$dissimilarity; write.csv(bc, file.path(OUT, "S_bc_between_groups.csv"), row.names = FALSE)

## 3. PERMANOVA (3 deposit types) per core, 9,999 permutations
pm <- do.call(rbind, lapply(names(cores), function(k) { x <- cores[[k]]; a <- adonis2(x$D ~ x$g, permutations = 9999)
  data.frame(core = k, df_model = a$Df[1], df_resid = a$Df[2], df_total = a$Df[3], SS_model = a$SumOfSqs[1], SS_resid = a$SumOfSqs[2], SS_total = a$SumOfSqs[3],
             MS_model = a$SumOfSqs[1] / a$Df[1], MS_resid = a$SumOfSqs[2] / a$Df[2], F = a$F[1], R2 = a$R2[1], p = a$`Pr(>F)`[1], perms = 9999) }))
write.csv(pm, file.path(OUT, "S_permanova_3groups.csv"), row.names = FALSE)

## 4. pairwise PERMANOVA between deposit types (exact enumeration of all relabellings), Holm within core, PERMDISP
pw <- do.call(rbind, lapply(names(cores), function(k) { x <- cores[[k]]
  res <- do.call(rbind, lapply(pairs, function(p) { i <- x$g %in% p; D <- as.matrix(x$D)[i, i]; g <- droplevels(x$g[i]); g <- factor(g, levels = p); D2 <- D^2
    n <- length(g); n1 <- sum(g == p[1]); obs <- pseudoF(D2, g)
    perm <- apply(combn(n, n1), 2, function(ix) { gg <- factor(rep(p[2], n), levels = p); gg[ix] <- p[1]; pseudoF(D2, gg)["F"] })
    bd <- tryCatch(permutest(betadisper(as.dist(D), g), permutations = 9999)$tab$`Pr(>F)`[1], error = function(e) NA)
    data.frame(core = k, comparison = paste(p, collapse = " vs "), n1 = n1, n2 = n - n1, F = obs["F"], R2 = (obs["SST"] - obs["SSW"]) / obs["SST"],
               p = mean(perm >= obs["F"] - 1e-12), min_p = 1 / choose(n, n1), PERMDISP_p = bd) }))
  res$p_holm <- p.adjust(res$p, "holm"); res }))
rownames(pw) <- NULL; write.csv(pw, file.path(OUT, "S_permanova_pairwise.csv"), row.names = FALSE)

## 5. combined hurricane vs nor'easter test (labels shuffled within core; statistic = sum of per-core pseudo-F)
hn <- lapply(cores, function(x) { i <- x$g %in% c("Hurricane", "Noreaster"); list(D2 = as.matrix(x$D)[i, i]^2, g = factor(as.character(x$g[i]), levels = c("Hurricane", "Noreaster"))) })
obs <- sum(sapply(hn, function(s) pseudoF(s$D2, s$g)["F"])); NP <- 99999
perm <- replicate(NP, sum(sapply(hn, function(s) pseudoF(s$D2, sample(s$g))["F"])))
pc <- (sum(perm >= obs) + 1) / (NP + 1); pf <- pw[pw$comparison == "Hurricane vs Noreaster", ]
fish <- pchisq(-2 * sum(log(pf$p)), 2 * nrow(pf), lower.tail = FALSE)
write.csv(data.frame(sum_F = obs, p_combined = pc, perms = NP, fisher_p = fish), file.path(OUT, "S_hurricane_vs_noreaster_combined.csv"), row.names = FALSE)

## 6. dbRDA (constrained by deposit type)
db <- do.call(rbind, lapply(names(cores), function(k) { x <- cores[[k]]; env_df <- data.frame(sediment_type = x$g)
  fit <- dbrda(x$D ~ sediment_type, data = env_df); an <- anova(fit, permutations = 9999); ev <- eigenvals(fit)
  data.frame(core = k, n = length(x$g), constrained_pct = 100 * fit$CCA$tot.chi / fit$tot.chi, CAP1_pct = 100 * ev[1] / sum(ev), CAP2_pct = 100 * ev[2] / sum(ev), F = an$F[1], p = an$`Pr(>F)`[1], perms = 9999) }))
write.csv(db, file.path(OUT, "S_dbRDA.csv"), row.names = FALSE)
print(bc, digits = 3); print(pm, digits = 4); print(pw, digits = 3); print(read.csv(file.path(OUT, "S_hurricane_vs_noreaster_combined.csv"))); print(db, digits = 3)
