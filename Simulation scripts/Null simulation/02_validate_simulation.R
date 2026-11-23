#!/usr/bin/env Rscript
# 02_validate_simulation.R
#
# Checks that the simulated phenotypes have the properties each scenario is
# supposed to have, on a sample of replicates. Run this before committing the
# server to 128 x 257 x 4 GWAS runs.
#
# Checks per scenario:
#   * realised Pearson correlation vs the target R      (S1, S2, S4)
#   * realised Spearman correlation vs the S1 Spearman  (all, incl. S3)
#   * marginal skewness / excess kurtosis               (all)
#   * missingness rate and pairwise complete overlap    (S4)
#
# Usage:
#   Rscript 02_validate_simulation.R --reps=1,2,3

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  corr      = "/HGCNT95FS/ADDLIE/work/Test/new_GWAS/residuals/T1_pooled_correlation_matrix.txt",
  outdir    = "/HGCNT95FS/ADDLIE/work/Test/JAGWAS_null",
  k         = "128",
  reps      = "1,2,3",
  scenarios = "all",
  miss_rate = "0.01",
  header    = "false"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k         <- as_int(opt$k, "k")
miss_rate <- as_num(opt$miss_rate, "miss_rate")
has_hdr   <- as_bool(opt$header, "header")
reps      <- as.integer(trimws(strsplit(opt$reps, ",", fixed = TRUE)[[1]]))
scen      <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])

log_open(file.path(opt$outdir, "logs", "02_validate_simulation.log"))
on.exit(log_close(), add = TRUE)

R <- read_corr_matrix(opt$corr, k = k)
# Spearman correlation implied by a Gaussian copula with correlation R.
S_target <- (6 / pi) * asin(R / 2)

read_pheno <- function(s, r) {
  f <- file.path(opt$outdir, "pheno", s, sprintf("rep%03d.pheno", r))
  if (!file.exists(f)) { msg("  MISSING: ", f); return(NULL) }
  dt <- fread(f, header = has_hdr, na.strings = c("NA", "-9"), data.table = FALSE)
  X <- as.matrix(dt[, -(1:2), drop = FALSE])
  storage.mode(X) <- "double"
  if (ncol(X) != k) { msg(sprintf("  WRONG WIDTH: %s has %d trait columns", f, ncol(X))); return(NULL) }
  X
}

lt <- lower.tri(R)
res <- vector("list", 0L)

for (s in scen) {
  msg("")
  msg("=== ", s, " ==========================================================")
  for (r in reps) {
    X <- read_pheno(s, r)
    if (is.null(X)) next

    n_finite <- sum(is.finite(X)); n_na <- sum(is.na(X))
    sk <- apply(X, 2, skewness); ku <- apply(X, 2, ex_kurtosis)

    P <- suppressWarnings(cor(X, use = "pairwise.complete.obs", method = "pearson"))
    S <- suppressWarnings(cor(X, use = "pairwise.complete.obs", method = "spearman"))

    dP <- max(abs(P[lt] - R[lt]))
    dS <- max(abs(S[lt] - S_target[lt]))
    rmseP <- sqrt(mean((P[lt] - R[lt])^2))

    msg(sprintf("rep %03d  n=%d  NA=%d (%.3f%%)", r, nrow(X), n_na,
                100 * n_na / length(X)))
    msg(sprintf("   skewness      median %8.3f   range [%.3f, %.3f]",
                median(sk, na.rm = TRUE), min(sk, na.rm = TRUE), max(sk, na.rm = TRUE)))
    msg(sprintf("   ex. kurtosis  median %8.3f   max %.3g",
                median(ku, na.rm = TRUE), max(ku, na.rm = TRUE)))

    if (s == "S3_cauchy") {
      # Pearson correlation is undefined for Cauchy (no finite variance); the
      # sample value is noise. Only the rank correlation is meaningful.
      msg(sprintf("   Spearman vs target: max|dS| = %.4f   (Pearson not applicable)", dS))
      msg(sprintf("   |X| max = %.4g   (extreme values are expected; check fastGWA copes)",
                  max(abs(X[is.finite(X)]))))
    } else if (grepl("^S5_t", s)) {
      # t marginals are rank-matched, not Pearson-matched: the copula transform
      # attenuates Pearson below R. Report the attenuation, do not treat it as
      # an error. Slope through the origin is the cleanest summary.
      slope <- sum(P[lt] * R[lt]) / sum(R[lt]^2)
      msg(sprintf("   Spearman vs target: max|dS| = %.4f   (this is the matched quantity)", dS))
      msg(sprintf("   Pearson attenuation: slope P ~ R = %.4f   (expected < 1 for heavy tails)",
                  slope))
    } else {
      msg(sprintf("   Pearson  vs R     : max|dP| = %.4f, RMSE = %.4f", dP, rmseP))
      msg(sprintf("   Spearman vs target: max|dS| = %.4f", dS))
    }

    if (s == "S4_normal_miss1pct") {
      obs  <- !is.na(X)
      rate <- 1 - mean(obs)
      n_i  <- colSums(obs)                    # per-trait GWAS sample sizes
      n_ij <- crossprod(obs)                  # pairwise complete counts
      ov   <- n_ij / nrow(X)
      # The phenotype correlation is still R (pairwise-complete Pearson is
      # unbiased). What shrinks is the correlation OF THE Z-SCORES, by the
      # sample-overlap factor n_ij / sqrt(n_i n_j).
      shrink <- n_ij / sqrt(outer(n_i, n_i))
      msg(sprintf("   missing rate        = %.5f   (target %.5f)", rate, miss_rate))
      msg(sprintf("   per-trait n         = %.1f    (predicted %.1f)",
                  mean(n_i), (1 - miss_rate) * nrow(X)))
      msg(sprintf("   pairwise overlap    = %.5f   (predicted %.5f)",
                  mean(ov[lt]), (1 - miss_rate)^2))
      msg(sprintf("   complete cases      = %.5f   (predicted %.5f, not the driver)",
                  mean(rowSums(obs) == k), (1 - miss_rate)^k))
      msg(sprintf("   z-corr shrink factor= %.5f   (predicted %.5f)",
                  mean(shrink[lt]), 1 - miss_rate))
    }

    res[[length(res) + 1L]] <- data.table(
      scenario = s, replicate = r, n = nrow(X), n_na = n_na,
      pct_na = 100 * n_na / length(X),
      skew_median = median(sk, na.rm = TRUE), kurt_median = median(ku, na.rm = TRUE),
      max_abs_dPearson = if (s == "S3_cauchy" || grepl("^S5_t", s)) NA_real_ else dP,
      rmse_Pearson = if (s == "S3_cauchy" || grepl("^S5_t", s)) NA_real_ else rmseP,
      max_abs_dSpearman = dS,
      max_abs_value = max(abs(X[is.finite(X)]))
    )
  }
}

if (length(res)) {
  out <- rbindlist(res)
  f <- file.path(opt$outdir, "qc", "validation_summary.txt")
  fwrite(out, f, sep = "\t")
  msg("")
  msg("Wrote ", f)
}

msg("")
msg("Expected values, for reference:")
msg("  S1  skewness ~ 0,     ex. kurtosis ~ 0,     max|dPearson| < ~0.03 at n = 23000")
msg("  S2  skewness ~ 1.75,  ex. kurtosis ~ 5.9,   max|dPearson| < ~0.05")
msg("  S3  skewness and kurtosis are meaningless (no finite moments);")
msg("      only max|dSpearman| should be small")
msg("  S4  missing rate ~ 0.01, pairwise overlap ~ 0.9801, complete cases ~ 0.276")
msg("  S5_t3  ex. kurtosis undefined (infinite 4th moment); S5_t5 ~ 6.")
msg("         Only max|dSpearman| should be small; Pearson is attenuated by design.")
msg("Done.")
