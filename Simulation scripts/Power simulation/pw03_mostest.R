#!/usr/bin/env Rscript
# pw03_mostest.R -- MOSTest (van der Meer et al. 2020) on the power simulation.
#
# MOSTest is not a summary-statistic method: it runs its own GWAS on
# individual-level data, including a second GWAS with the genotypes permuted
# across individuals, and uses the permuted run to set up the test. Its three
# defining steps, as in mostest.m (github.com/precimed/mostest), are
# reproduced here from the simulated individual-level data:
#
# 1. CORRELATION. C0 = correlation of the permuted-genotype z-scores across
#    variants. With the phenotypes fixed, those null z-scores have covariance
#    equal to the sample correlation of the phenotypes, so C0 is the data
#    set's own phenotypic correlation -- the matrix used here.
#
# 2. REGULARIZATION. Eigenvalues of C0 below its 10th largest are raised to the
#    10th largest (mostest.m: max_lambda = s(min(10, length(s))),
#    C0_reg = U * diag(max(max_lambda, s)) * U'). The statistic is
#    z' C0_reg^-1 z. This down-weights exactly the low-variance directions that
#    JAGWAS's full inverse up-weights, which is the main reason the two can
#    differ in power.
#
# 3. NULL CALIBRATION. A gamma distribution is fitted by maximum likelihood to
#    the statistic of null variants (MOSTest: all permuted variants), and each
#    variant's p-value is its upper tail. Here the null variants are the K per
#    data set simulated in pw02, tested against pure-noise phenotypes -- the
#    same null MOSTest's permutation produces.
#
# MOSTest also rank-normalizes phenotypes by default; the simulated phenotypes
# are normal, so that step changes nothing and is omitted. Check step 2 against
# the mostest.m version you cite; --keep sets the eigenvalue index.
#
# Writes <root>/pvalues/mostest/<scenario>.txt and
#        <root>/pvalues/mostest/gamma_fit.txt.
#
# Usage:
#   Rscript pw03_mostest.R            (the null variants are always done first)

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

M <- "mostest"
S <- method_setup(M, list(keep = "10"))
keep <- as_int(S$opt$keep, "keep")

# Regularized inverse, as in mostest.m.
reg_inverse <- function(R) {
  e <- eigen(R, symmetric = TRUE)
  s <- pmax(e$values, e$values[min(keep, length(e$values))])
  e$vectors %*% (t(e$vectors) / s)
}
mostest_stat <- function(z, R) {
  A <- reg_inverse(R)
  list(stat = rowSums((z %*% A) * z), p = NA_real_)
}

# Maximum-likelihood gamma fit (MATLAB's fitdist(x, 'gamma')): the shape solves
# log(a) - digamma(a) = log(mean(x)) - mean(log(x)); scale = mean(x) / a.
fit_gamma <- function(x) {
  x <- x[is.finite(x) & x > 0]
  s <- log(mean(x)) - mean(log(x))
  a <- uniroot(function(a) log(a) - digamma(a) - s, c(1e-4, 1e7), tol = 1e-12)$root
  c(shape = a, scale = mean(x) / a)
}

fgam <- file.path(S$opt$root, "pvalues", M, "gamma_fit.txt")

# Null variants first: they define the gamma distribution.
if (S$overwrite || !file.exists(fgam) || !file.exists(pval_file(S$opt$root, M, "null"))) {
  inp <- test_inputs(S$opt$root, S$design, S$D, "null")
  r <- timed(by_dataset(inp, mostest_stat))
  gfit <- fit_gamma(r$value$stat)
  r$value[, p := pgamma(stat, shape = gfit[["shape"]], scale = gfit[["scale"]], lower.tail = FALSE)]
  dir.create(dirname(fgam), recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(shape = gfit[["shape"]], scale = gfit[["scale"]], keep = keep,
                    n_null = nrow(r$value)), fgam, sep = "\t")
  save_pvals(r$value, S$opt$root, M, "null", r$secs, r$mem_mb)
  msg(sprintf("Gamma fit to %d null statistics: shape %.4f, scale %.4f (mean %.2f; a chi2_%d would have mean %d)",
              nrow(r$value), gfit[["shape"]], gfit[["scale"]], prod(gfit), S$design$k, S$design$k))
  report(M, "null", r$value, S$design$alpha, r$secs)
}
gfit <- unlist(fread(fgam)[, .(shape, scale)])

for (scen in setdiff(S$scen, "null")) {
  if (!S$overwrite && file.exists(pval_file(S$opt$root, M, scen))) next
  inp <- test_inputs(S$opt$root, S$design, S$D, scen)
  r <- timed(by_dataset(inp, mostest_stat))
  r$value[, p := pgamma(stat, shape = gfit[["shape"]], scale = gfit[["scale"]], lower.tail = FALSE)]
  save_pvals(r$value, S$opt$root, M, scen, r$secs, r$mem_mb)
  report(M, scen, r$value, S$design$alpha, r$secs)
}
msg("Done.")
