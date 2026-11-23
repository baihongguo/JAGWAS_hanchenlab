#!/usr/bin/env Rscript
# 03_check_scenario_misspecification.R
#
# Answers one question per scenario: does the REALISED phenotype correlation
# differ from the R that JAGWAS tests with, badly enough to miscalibrate the
# test on its own?
#
# This matters because R is ill-conditioned. A correlation mismatch that looks
# negligible entry-wise can be severe if it aligns with R's small eigenvectors
# -- the same mechanism that makes S4 blow up. Entry-wise RMSE cannot tell you
# this; the eigenvalues of R^-1 P can.
#
#   z ~ MVN(0, P)  tested as  T = z' R^-1 z   =>   T ~ sum_k lambda_k chi2_1
#   lambda = eigenvalues of R^(-1/2) P R^(-1/2)
#   Perfectly matched  =>  all lambda = 1, E[T] = k, tail = nominal.
#
# Use it two ways:
#   * S1 / S2 -- diagnostic. If E[T] is far from k, the scenario is confounded
#     by a correlation mismatch and is no longer a clean test of the marginal
#     distribution. Regenerate with a smaller --sigma_log.
#   * S4 -- confirmation. This should empirically reproduce the analytic
#     prediction from 00_check_correlation_matrix.R.
#
# Not meaningful for S3_cauchy (Pearson correlation does not exist) or the
# opt-in S5_t* scenarios (rank-matched by design, Pearson attenuated on purpose).
#
# Usage:
#   Rscript 03_check_scenario_misspecification.R --scenarios=S1_normal,S2_lognormal --reps=1

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
  reps      = "1",
  scenarios = "S1_normal,S2_lognormal,S4_normal_miss1pct",
  alpha     = "5e-8",
  header    = "false"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k       <- as_int(opt$k, "k")
alpha   <- as_num(opt$alpha, "alpha")
has_hdr <- as_bool(opt$header, "header")
reps    <- as.integer(trimws(strsplit(opt$reps, ",", fixed = TRUE)[[1]]))
scen    <- trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])

log_open(file.path(opt$outdir, "logs", "03_check_scenario_misspecification.log"))
on.exit(log_close(), add = TRUE)

R  <- read_corr_matrix(opt$corr, k = k)
eR <- eigen(R, symmetric = TRUE)
if (min(eR$values) <= 0) stop("R is not positive definite; cannot form R^(-1/2).")
Rmhalf <- eR$vectors %*% (diag(1 / sqrt(eR$values)) %*% t(eR$vectors))
qcrit  <- qchisq(alpha, df = k, lower.tail = FALSE)

msg(sprintf("Testing matrix R: lambda_min = %.6g, condition = %.4g, tr(R^-1) = %.1f",
            min(eR$values), max(eR$values) / min(eR$values), sum(1 / eR$values)))
msg(sprintf("Nominal threshold alpha = %.3g, chi2_%d critical value = %.3f", alpha, k, qcrit))

res <- vector("list", 0L)

for (s in scen) {
  if (s == "S3_cauchy" || grepl("^S5_t", s)) {
    msg("")
    msg("Skipping ", s, ": Pearson correlation is not the matched (or defined) quantity.")
    next
  }
  msg("")
  msg("=== ", s, " ==========================================================")

  for (r in reps) {
    f <- file.path(opt$outdir, "pheno", s, sprintf("rep%03d.pheno", r))
    if (!file.exists(f)) { msg("  MISSING: ", f); next }
    dt <- fread(f, header = has_hdr, na.strings = c("NA", "-9"), data.table = FALSE)
    X  <- as.matrix(dt[, -(1:2), drop = FALSE]); storage.mode(X) <- "double"

    P <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))
    if (anyNA(P)) { msg("  Realised correlation has NA; skipping."); next }

    # For S4 the z-score correlation is the phenotype correlation scaled by the
    # sample-overlap factor n_ij / sqrt(n_i n_j); for complete data that factor
    # is 1 and P is used as-is.
    obs  <- !is.na(X)
    n_i  <- colSums(obs)
    ovl  <- crossprod(obs) / sqrt(outer(n_i, n_i))
    Pz   <- P * ovl
    diag(Pz) <- 1

    M   <- Rmhalf %*% Pz %*% Rmhalf
    lam <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values

    ET  <- sum(lam)
    c2  <- 2 * sum(lam^2)
    a_s <- c2 / (2 * ET); df_e <- 2 * ET^2 / c2
    p_sw <- pchisq(qcrit / a_s, df = df_e, lower.tail = FALSE)
    p_dav <- if (requireNamespace("CompQuadForm", quietly = TRUE)) {
      CompQuadForm::davies(qcrit, lambda = lam, h = rep(1, k), acc = 1e-12)$Qq
    } else NA_real_

    msg(sprintf("rep %03d", r))
    msg(sprintf("   max|P - R|           = %.4f    RMSE = %.4f",
                max(abs(P - R)), sqrt(mean((P - R)^2))))
    msg(sprintf("   lambda range         = [%.4f, %.4f]   (1.0 = perfectly matched)",
                min(lam), max(lam)))
    msg(sprintf("   E[T]                 = %.3f    (null expectation %d, ratio %.4f)",
                ET, k, ET / k))
    msg(sprintf("   type I error at %.1g  = %.4g   (ratio to nominal %.3f)",
                alpha, p_sw, p_sw / alpha))
    if (!is.na(p_dav))
      msg(sprintf("   Davies (exact)       = %.4g   (ratio %.3f)", p_dav, p_dav / alpha))

    ratio <- if (is.na(p_dav)) p_sw / alpha else p_dav / alpha
    verdict <- if (ratio < 2) {
      "CLEAN -- the correlation mismatch does not confound this scenario."
    } else if (ratio < 20) {
      "MARGINAL -- report the mismatch alongside the result."
    } else {
      "CONFOUNDED -- inflation here is driven by the correlation mismatch, not the marginals."
    }
    msg("   => ", verdict)

    res[[length(res) + 1L]] <- data.table(
      scenario = s, replicate = r, max_abs_dP = max(abs(P - R)),
      rmse_dP = sqrt(mean((P - R)^2)), lambda_min = min(lam), lambda_max = max(lam),
      E_T = ET, mean_ratio = ET / k, type1_satterthwaite = p_sw,
      type1_davies = p_dav, tail_ratio = ratio
    )
  }
}

if (length(res)) {
  out <- rbindlist(res)
  f <- file.path(opt$outdir, "qc", "misspecification_summary.txt")
  fwrite(out, f, sep = "\t")
  msg("")
  msg("Wrote ", f)
}

msg("")
msg("Interpretation. For S1 and S2 a tail ratio near 1 means the scenario is a clean")
msg("test of the MARGINAL distribution. A large ratio means you are measuring a")
msg("correlation mismatch instead, and should regenerate with a smaller --sigma_log.")
msg("For S4 a large ratio is the expected result and should match 00's prediction.")
msg("Done.")
