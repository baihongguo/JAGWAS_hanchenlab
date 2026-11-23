#!/usr/bin/env Rscript
# pw03_multiabel.R -- MultiABEL on the power-simulation summary statistics.
#
# MultiABEL's summary-statistic test reconstructs the classical MANOVA for one
# variant from the marginal statistics and the phenotypic correlation matrix.
# With a single predictor the MANOVA statistics (Pillai, Wilks, Hotelling-
# Lawley, Roy) coincide, and Pillai's trace is the squared multiple
# correlation of the genotype on the k traits:
#
#   r_t = t_t / sqrt(t_t^2 + n - 2)          marginal genotype-trait correlation
#   V   = r' Rhat^-1 r
#   F   = (V / k) / ((1 - V) / (n - k - 1))  ~ F(k, n - k - 1) under the null
#
# This is the exact finite-sample version of JAGWAS's large-sample chi2_k test,
# so the two should agree closely; the F reference makes MultiABEL very
# slightly more conservative.
#
# Rhat is the phenotypic correlation of the variant's own data set, as for
# JAGWAS. Writes <root>/pvalues/multiabel/<scenario>.txt.
#
# To confirm against the MultiABEL package on a few variants, pass the same
# summary statistics and correlation matrix to MultiABEL::load.summary() and
# MultiABEL::MultiSummary() and compare p-values.
#
# Usage:
#   Rscript pw03_multiabel.R

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

M <- "multiabel"
S <- method_setup(M)
k <- S$design$k; n <- S$design$n

multiabel_test <- function(z, R) {
  r <- z / sqrt(z^2 + n - 2)
  V <- colSums(forwardsolve(t(chol(R)), t(r))^2)
  F <- (V / k) / ((1 - V) / (n - k - 1))
  list(stat = F, p = pf(F, k, n - k - 1, lower.tail = FALSE))
}

for (scen in S$scen) {
  if (!S$overwrite && file.exists(pval_file(S$opt$root, M, scen))) next
  inp <- test_inputs(S$opt$root, S$design, S$D, scen)
  r <- timed(by_dataset(inp, multiabel_test))
  save_pvals(r$value, S$opt$root, M, scen, r$secs, r$mem_mb)
  report(M, scen, r$value, S$design$alpha, r$secs)
}
msg("Done.")
