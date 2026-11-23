#!/usr/bin/env Rscript
# pw03_jagwas.R -- JAGWAS on the power-simulation summary statistics.
#
#   T = z' Rhat^-1 z,   p = P(chi2_k > T)
#
# with Rhat the phenotypic correlation of the variant's own data set, the
# matrix JAGWAS uses in practice. Computed via a Cholesky solve; this is the
# same computation checked against the JAGWAS binary in the null study (06:
# maximum relative difference 5e-6, the binary's output rounding).
#
# --binary_check=N reruns the first N data sets of the first scenario through
# the JAGWAS binary itself (causal plus null variants, the data set's own
# Rhat) and reports the largest relative difference in -log10 p.
#
# Writes <root>/pvalues/jagwas/<scenario>.txt (ds, snp, stat, p); snp 0 is the
# causal variant, 1..K the null variants.
#
# Usage:
#   Rscript pw03_jagwas.R [--scenarios=all] [--binary_check=3]

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

M <- "jagwas"
S <- method_setup(M, list(binary_check = "0",
                          bin = "/data5/bguo1/JAGWAS_hanchenlab/JAGWAS"))
k <- S$design$k

jagwas_test <- function(z, R) {
  L <- t(chol(R))
  T <- colSums(forwardsolve(L, t(z))^2)
  list(stat = T, p = pchisq(T, df = k, lower.tail = FALSE))
}

for (scen in S$scen) {
  if (!S$overwrite && file.exists(pval_file(S$opt$root, M, scen))) next
  inp <- test_inputs(S$opt$root, S$design, S$D, scen)
  r <- timed(by_dataset(inp, jagwas_test))
  save_pvals(r$value, S$opt$root, M, scen, r$secs, r$mem_mb)
  report(M, scen, r$value, S$design$alpha, r$secs)
}

# ------------------------------------------------ check against the binary --

n_check <- as_int(S$opt$binary_check, "binary_check")
if (n_check > 0) {
  bin <- S$opt$bin
  if (!file.exists(bin)) stop("JAGWAS binary not found: ", bin)
  scen <- S$design$scen$name[1]
  cz   <- test_inputs(S$opt$root, S$design, S$D, scen)
  nz   <- test_inputs(S$opt$root, S$design, S$D, "null")
  tmp  <- tempfile("jagwas_check_"); dir.create(tmp)
  worst <- 0
  for (j in head(sort(unique(cz$ds)), n_check)) {
    z  <- rbind(cz$z[cz$ds == j, , drop = FALSE], nz$z[nz$ds == j, , drop = FALSE])
    R  <- cz$rhat(j)
    ids <- c(sprintf("ds%04d_causal", j), sprintf("ds%04d_null%03d", j, nz$snp[nz$ds == j]))
    files <- file.path(tmp, sprintf("QT%d.txt", seq_len(k) - 1L))
    for (tt in seq_len(k))
      fwrite(data.table(CHR = 1L, SNP = ids, POS = seq_along(ids), A1 = "A", A2 = "G",
                        N = S$design$n, AF1 = 0.3, BETA = z[, tt], SE = 1),
             files[tt], sep = "\t")
    cf <- file.path(tmp, "R.txt")
    writeLines(apply(R, 1, function(x) paste(sprintf("%.17g", x), collapse = " ")), cf)
    of <- file.path(tmp, "out.txt")
    # Arguments are read by position; this order is required.
    status <- system2(bin, c("--outputFilePath", of, "--cor_matrix", cf, "--nrow", "50000",
                             "--MAF", "0", "--score_test", "0", "--beta_se", "1",
                             "--logP", "1", "--delim", "t", "--fileNames", files),
                      stdout = file.path(tmp, "log.txt"), stderr = file.path(tmp, "log.txt"))
    if (status != 0 || !file.exists(of)) stop("JAGWAS binary failed; see ", file.path(tmp, "log.txt"))
    x <- -fread(of)$log10_P
    y <- -pchisq(colSums(forwardsolve(t(chol(R)), t(z))^2), df = k,
                 lower.tail = FALSE, log.p = TRUE) / log(10)
    d <- max(abs(x - y) / pmax(1, y))
    worst <- max(worst, d)
    msg(sprintf("  binary check, data set %d (%d variants): max relative difference in -log10 p %.2e",
                j, length(x), d))
  }
  msg(sprintf("Binary check: worst relative difference %.2e (%s)", worst,
              if (worst < 1e-4) "PASS" else "FAIL"))
  unlink(tmp, recursive = TRUE)
}
msg("Done.")
