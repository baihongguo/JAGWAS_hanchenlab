#!/usr/bin/env Rscript
# 00_check_correlation_matrix.R
#
# Inspects the 128 x 128 pooled residual correlation matrix BEFORE any
# phenotypes are generated, and reports the quantities that determine whether
# each of the four null scenarios will show anything:
#
#   * eigenvalue spectrum and condition number      -> numerical stability of R^-1
#   * tr(R^-1)                                      -> analytic prediction for S4
#   * Li & Ji (2005) effective number of tests      -> the "77 vs 128" number
#   * feasibility of the log-normal Pearson match   -> whether sigma_log is usable
#
# Run this first. If it fails or warns, fix the matrix before simulating.
#
# # Usage:
#   Rscript 00_check_correlation_matrix.R \
#     --corr=/HGCNT95FS/ADDLIE/work/Test/new_GWAS/residuals/T1_pooled_correlation_matrix.txt \
#     --outdir=/HGCNT95FS/ADDLIE/work/Test/JAGWAS_null

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
  sigma_log = "0.5",
  miss_rate = "0.01"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k         <- as_int(opt$k, "k")
sigma_log <- as_num(opt$sigma_log, "sigma_log")
miss_rate <- as_num(opt$miss_rate, "miss_rate")

qcdir <- file.path(opt$outdir, "qc")
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$outdir, "logs", "00_check_correlation_matrix.log"))
on.exit(log_close(), add = TRUE)

msg("Reading correlation matrix: ", opt$corr)
R <- read_corr_matrix(opt$corr, k = k)
msg(sprintf("Parsed %d x %d. max|R - R'| = %.3g, max|diag - 1| = %.3g",
            k, k, attr(R, "max_asymmetry"), attr(R, "max_diag_dev")))

# ------------------------------------------------------------------ spectrum --

ev   <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
lmin <- min(ev); lmax <- max(ev)
cond <- lmax / lmin

msg("")
msg("--- Eigenvalue spectrum -------------------------------------------------")
msg(sprintf("  lambda_min      = %.6g", lmin))
msg(sprintf("  lambda_max      = %.6g", lmax))
msg(sprintf("  condition number= %.6g", cond))
msg(sprintf("  n eigenvalues < 1e-6 : %d", sum(ev < 1e-6)))
msg(sprintf("  n eigenvalues < 1e-3 : %d", sum(ev < 1e-3)))

if (lmin <= 0) {
  msg("  ** R IS NOT POSITIVE DEFINITE. JAGWAS needs R^-1; use the generalized")
  msg("     inverse with df = rank(R), and note this in the robustness section.")
} else if (cond > 1e4) {
  msg("  ** Condition number is large. R^-1 is numerically delicate and the S4")
  msg("     missingness scenario is likely to show real inflation (see below).")
}

# ------------------------------------------- tr(R^-1) and the S4 prediction --

# Under phenotype-specific missingness at rate m, independent across traits,
# the per-trait GWAS sample sizes are n_i ~ (1-m)N and the pairwise overlap is
# n_ij ~ (1-m)^2 N, so
#       corr(z_i, z_j) = r_ij * n_ij / sqrt(n_i n_j) = (1-m) * r_ij,
# i.e. the z-score covariance is  Sigma = (1-m) R + m I  while JAGWAS still
# tests with R. Then T = z' R^-1 z ~ sum_k lambda_k chi2_1 with
#       lambda_k = (1-m) + m / mu_k,     mu_k = eigenvalues of R,
# and  E[T] = (1-m) k + m * tr(R^-1).
# If R is well conditioned, tr(R^-1) ~ k and E[T] = k exactly: no inflation.

Rinv    <- tryCatch(solve(R), error = function(e) NULL)
trRinv  <- if (is.null(Rinv)) sum(1 / ev[ev > 1e-12]) else sum(diag(Rinv))
m       <- miss_rate
lambda  <- (1 - m) + m / ev
ET      <- (1 - m) * k + m * trRinv

msg("")
msg("--- Scenario S4: analytic prediction ------------------------------------")
msg(sprintf("  missingness rate m           = %.4f", m))
msg(sprintf("  per-trait n / N              = %.4f", 1 - m))
msg(sprintf("  pairwise overlap n_ij / N    = %.4f", (1 - m)^2))
msg(sprintf("  complete-case fraction       = %.4f  (0.99^%d; NOT the driver)", (1 - m)^k, k))
msg(sprintf("  implied z-correlation shrink = %.4f x r_ij", 1 - m))
msg(sprintf("  tr(R^-1)                     = %.4f   (equals %d if R = I)", trRinv, k))
msg(sprintf("  E[T] under S4                = %.4f   (null expectation is %d)", ET, k))
msg(sprintf("  inflation in mean            = %.4f", ET / k))
msg(sprintf("  max mixture weight lambda_max= %.4f", max(lambda)))
if (ET / k < 1.01) {
  msg("  => S4 is predicted to be essentially calibrated. Report the analytic")
  msg("     result alongside the empirical curve; a null result here is a result.")
} else {
  msg("  => S4 is predicted to be INFLATED. This is driven by the small")
  msg("     eigenvalues of R, and ties directly to the robustness analysis.")
}

# Exact tail probability of the mixture at genome-wide significance, by the
# Satterthwaite moment match (cheap sanity check; Davies' method is exact if the
# CompQuadForm package is available).
c1 <- sum(lambda); c2 <- 2 * sum(lambda^2)
a_scale <- c2 / (2 * c1); df_eff <- 2 * c1^2 / c2
qcrit <- qchisq(5e-8, df = k, lower.tail = FALSE)
p_sw  <- pchisq(qcrit / a_scale, df = df_eff, lower.tail = FALSE)
msg(sprintf("  P(T > chi2_%d critical value at 5e-8) under S4:", k))
msg(sprintf("    Satterthwaite approx = %.4g  (nominal 5e-8, ratio %.3f)", p_sw, p_sw / 5e-8))
if (requireNamespace("CompQuadForm", quietly = TRUE)) {
  p_dav <- CompQuadForm::davies(qcrit, lambda = lambda, h = rep(1, k), acc = 1e-12)$Qq
  msg(sprintf("    Davies (exact)       = %.4g  (ratio %.3f)", p_dav, p_dav / 5e-8))
} else {
  msg("    (install CompQuadForm for the exact Davies tail probability)")
}

# ------------------------------------------------------ effective test count --

meff <- li_ji_meff(R)
msg("")
msg("--- Effective number of independent tests -------------------------------")
msg(sprintf("  Li & Ji (2005) M_eff = %.2f  of %d traits", meff, k))
msg(sprintf("  Bonferroni thresholds: 5e-8/%d = %.4g   vs   5e-8/%.0f = %.4g",
            k, 5e-8 / k, meff, 5e-8 / meff))

# ------------------------------------------- log-normal Pearson feasibility --

rmin <- min(R[lower.tri(R)]); rmax <- max(R[lower.tri(R)])
msg("")
msg("--- Scenario S2: log-normal Pearson matching ----------------------------")
msg(sprintf("  off-diagonal range of R  = [%.4f, %.4f]", rmin, rmax))
msg(sprintf("  requested sigma_log      = %.3f  (marginal skewness %.3f)",
            sigma_log, (exp(sigma_log^2) + 2) * sqrt(exp(sigma_log^2) - 1)))

smax <- lognormal_max_sigma(rmin)
msg(sprintf("  max sigma_log attainable = %.4f  (set by r_min = %.4f)", smax, rmin))

rho_ln <- lognormal_latent_rho(R, sigma_log)
if (is.null(rho_ln) || any(!is.finite(rho_ln)) || min(rho_ln) < -1) {
  msg("  ** sigma_log is TOO LARGE: the target Pearson correlations cannot be")
  msg(sprintf("     reproduced with log-normal marginals. Use --sigma_log <= %.3f", smax))
} else {
  diag(rho_ln) <- 1
  ev_ln <- eigen(rho_ln, symmetric = TRUE, only.values = TRUE)$values
  msg(sprintf("  latent rho range         = [%.4f, %.4f]",
              min(rho_ln[lower.tri(rho_ln)]), max(rho_ln[lower.tri(rho_ln)])))
  msg(sprintf("  latent lambda_min        = %.6g %s", min(ev_ln),
              if (min(ev_ln) <= 0) "(will be projected to nearest PD matrix)" else ""))
  msg(sprintf("  max |rho_latent - r|     = %.4f", max(abs(rho_ln - R))))
}

# ---------------------------------------------------------------- artefacts --

fwrite(data.table(index = seq_len(k), eigenvalue = ev),
       file.path(qcdir, "R_eigenvalues.txt"), sep = "\t")

summary_dt <- data.table(
  quantity = c("k", "lambda_min", "lambda_max", "condition_number", "tr_Rinv",
               "li_ji_meff", "miss_rate", "E_T_under_S4", "mean_inflation_S4",
               "tail_ratio_S4_satterthwaite", "r_offdiag_min", "r_offdiag_max",
               "sigma_log_requested", "sigma_log_max_feasible"),
  value    = c(k, lmin, lmax, cond, trRinv, meff, m, ET, ET / k,
               p_sw / 5e-8, rmin, rmax, sigma_log, smax)
)
fwrite(summary_dt, file.path(qcdir, "R_summary.txt"), sep = "\t")

msg("")
msg("Wrote ", file.path(qcdir, "R_eigenvalues.txt"))
msg("Wrote ", file.path(qcdir, "R_summary.txt"))
msg("Done.")
