#!/usr/bin/env Rscript
# 08_check_layout.R
#
# Decides how tstat.f32 is laid out, by reading a sample of variants under both
# hypotheses and comparing the resulting cross-trait covariance with Rhat.
#
#   trait-major   value (variant v, trait t) at index (t-1)*n_variant + v
#                 -- what 05 assumed
#   variant-major value (variant v, trait t) at index (v-1)*n_trait + t
#
# The right layout gives z-scores whose cross-trait covariance is Rhat, so
# E[T]/k = tr(Rhat^-1 C)/k is ~1. The wrong one shuffles unrelated variants and
# traits into each "variant", so C is ~I and E[T]/k is ~tr(Rhat^-1)/k -- 64 for
# this matrix, which puts every variant past any significance threshold.
#
# Reads about 100 MB per hypothesis; takes seconds.
#
# Usage:
#   Rscript 08_check_layout.R --scenario=S1_normal --rep=1

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  tstat_root = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current",
  inbase     = "/data5/bguo1/JAGWAS_null_input",
  scenario   = "S1_normal",
  rep        = "1",
  k          = "128",
  n_variant  = "8931083",
  n_check    = "200000",     # variants sampled
  first      = "1"           # first variant of the sample
)
opt     <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k       <- as_int(opt$k, "k")
nv      <- as_num(opt$n_variant, "n_variant")
n_check <- as_int(opt$n_check, "n_check")
first   <- as_num(opt$first, "first")
r       <- as_int(opt$rep, "rep")

tag  <- sprintf("rep%03d", r)
fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
frh  <- file.path(opt$inbase, opt$scenario, sprintf("Rhat_rep%03d.txt", r))

cat(sprintf("binary : %s\n", fbin))
cat(sprintf("Rhat   : %s\n", frh))
if (!file.exists(fbin)) stop("Missing ", fbin)
if (!file.exists(frh))  stop("Missing ", frh)
cat(sprintf("size   : %.0f bytes = %.0f values = %.0f variants x %d traits\n",
            file.size(fbin), file.size(fbin) / 4, nv, k))

Rhat <- as.matrix(fread(frh, header = FALSE, showProgress = FALSE))
dimnames(Rhat) <- NULL
eR   <- eigen(Rhat, symmetric = TRUE)
Rinv <- eR$vectors %*% (t(eR$vectors) / eR$values)
lt   <- lower.tri(Rhat)
cat(sprintf("Rhat   : lambda_min %.6g, tr(Rhat^-1) %.1f, mean |off-diagonal| %.4f\n\n",
            min(eR$values), sum(1 / eR$values), mean(abs(Rhat[lt]))))

report <- function(name, Z) {
  Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
  n <- nrow(Z)
  C <- crossprod(Z) / n
  ET <- sum(diag(Rinv %*% C)) / k
  cr <- cov2cor(C)
  cat(sprintf("--- %s (%d variants) ---\n", name, n))
  cat(sprintf("   mean z^2 per trait : %.4f  [%.4f, %.4f]\n",
              mean(diag(C)), min(diag(C)), max(diag(C))))
  cat(sprintf("   mean |off-diagonal correlation| : %.4f   (Rhat: %.4f)\n",
              mean(abs(cr[lt])), mean(abs(Rhat[lt]))))
  cat(sprintf("   max |cor(Z) - Rhat| : %.4f\n", max(abs(cr - Rhat)[lt])))
  cat(sprintf("   E[T]/k = tr(Rhat^-1 C)/k : %.4f   (1 = correct layout)\n\n", ET))
  ET
}

con <- file(fbin, "rb")
on.exit(close(con))

# trait-major: trait t is one contiguous block of n_variant values
Zt <- matrix(0, n_check, k)
for (t in seq_len(k)) {
  seek(con, where = ((t - 1) * nv + (first - 1)) * 4, origin = "start")
  Zt[, t] <- readBin(con, "numeric", n = n_check, size = 4, endian = "little")
}
et_trait <- report("trait-major  (t-1)*n_variant + v", Zt)

# variant-major: each variant's 128 traits are consecutive
seek(con, where = (first - 1) * k * 4, origin = "start")
Zv <- matrix(readBin(con, "numeric", n = n_check * k, size = 4, endian = "little"),
             nrow = n_check, ncol = k, byrow = TRUE)
et_var <- report("variant-major  (v-1)*n_trait + t", Zv)

cat("VERDICT: ")
if (abs(et_trait - 1) < 0.1 && abs(et_var - 1) >= 0.1) {
  cat("trait-major. 05 is reading the file correctly; look elsewhere.\n")
} else if (abs(et_var - 1) < 0.1 && abs(et_trait - 1) >= 0.1) {
  cat("VARIANT-major. 05 is reading the file wrongly and must be changed.\n")
} else if (abs(et_trait - 1) < 0.1 && abs(et_var - 1) < 0.1) {
  cat("ambiguous -- both look calibrated, which should not happen. Re-run with --first=4000000.\n")
} else {
  cat("NEITHER layout matches Rhat.\n")
  cat("  If the traits are strongly correlated under one of them (mean |off-diagonal|\n")
  cat("  close to Rhat's) but max |cor(Z) - Rhat| is large, the traits are in a\n")
  cat("  different ORDER from the phenotype columns.\n")
  cat("  If both give mean |off-diagonal| ~0, neither layout is right.\n")
}
