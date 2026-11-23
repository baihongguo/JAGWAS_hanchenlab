#!/usr/bin/env Rscript
# pw03_minp.R -- minP on the power-simulation summary statistics, two ways:
#
#   minp_bonf   the single-phenotype workflow of the paper: significant when the
#               smallest per-trait p-value is below alpha / k (5e-8 / 128)
#   minp_eff    the correlation-adjusted version the referees asked for:
#               alpha / M_eff, with M_eff from Li & Ji (2005) on R (77 for T1)
#
# Each is written as an adjusted p-value, min(1, k * p_min) and
# min(1, M_eff * p_min), so every method is judged at the same alpha.
# Per-trait p-values use the t distribution with n - 2 df.
#
# Writes <root>/pvalues/minp_bonf/<scenario>.txt and
#        <root>/pvalues/minp_eff/<scenario>.txt.
#
# Usage:
#   Rscript pw03_minp.R

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

S <- method_setup("minp")
k <- S$design$k; n <- S$design$n; meff <- S$design$meff
msg(sprintf("Bonferroni over k = %d; effective number of tests M_eff = %.2f", k, meff))

for (scen in S$scen) {
  if (!S$overwrite && file.exists(pval_file(S$opt$root, "minp_bonf", scen)) &&
      file.exists(pval_file(S$opt$root, "minp_eff", scen))) next
  inp <- test_inputs(S$opt$root, S$design, S$D, scen)
  r <- timed({
    pmin_ <- apply(2 * pt(-abs(inp$z), df = n - 2), 1, min)
    data.table(ds = inp$ds, snp = inp$snp, stat = pmin_)
  })
  for (v in c("bonf", "eff")) {
    mult <- if (v == "bonf") k else meff
    dt <- copy(r$value)[, p := pmin(1, mult * stat)]
    save_pvals(dt, S$opt$root, paste0("minp_", v), scen, r$secs, r$mem_mb)
    report(paste0("minp_", v), scen, dt, S$design$alpha, r$secs)
  }
}
msg("Done.")
