#!/usr/bin/env Rscript
# pw02_simulate_sumstats.R
#
# Simulates the individual-level data sets defined by pw01 and runs the
# per-trait GWAS, writing summary statistics:
#
#   <root>/sumstats/<scenario>/ss_<range>.rds   causal variant: t, beta, se per trait
#   <root>/sumstats/<scenario>/ss_<range>.tsv.gz  the same in long format (--tsv=true)
#   <root>/data/ds_<range>.rds                  per data set, for the methods:
#                                               E'E, E'g, g'g, and the null
#                                               variants' t statistics
#
# DATA SET j (seeded by j alone, so it is identical in every scenario):
#   E  : n x 128 noise, rows ~ MVN(0, R)
#   g  : causal genotype ~ Binomial(2, MAF_j), standardized
#   G0 : --K null genotypes, independent of everything
#   Y  = g beta' + E, with beta the scenario's effect vector for data set j
#
# THE GWAS is ordinary least squares of each trait on each variant, reported
# as t statistics. With simulated unrelated individuals there is no kinship
# matrix, and GMMAT without one fits exactly this linear model, so its score
# test gives the same statistics; running GMMAT 128 x 1000 x 25 times would
# only be slower. Per trait t:
#     b_t = g'y_t / g'g,   se_t = sqrt((y_t'y_t - b_t^2 g'g) / (n-2) / g'g),
#     t_t = b_t / se_t
# which is exactly the slope t statistic of lm(y_t ~ g). Only the causal term
# changes between scenarios, so each scenario reuses E'E, E'g and g'g and costs
# almost nothing; --check=true refits lm() on the first data set and compares.
#
# The null variants are tested against the pure-noise phenotypes E, once. They
# give every method its calibration check and MOSTest its null distribution.
#
# Usage (four parallel jobs):
#   for r in "1 250" "251 500" "501 750" "751 1000"; do set -- $r
#     nohup Rscript pw02_simulate_sumstats.R --ds_start=$1 --ds_end=$2 > pw02_$1.log 2>&1 &
#   done

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})
RNGkind("Mersenne-Twister", "Inversion")

defaults <- list(
  root      = POWER_ROOT,
  ds_start  = "1",
  ds_end    = "0",          # 0: last data set
  tsv       = "false",
  check     = "false",
  check_only = "false",     # simulate data set ds_start, run the check, write nothing
  overwrite = "false",
  threads   = "1"
)
opt       <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
design    <- load_design(opt$root)
k <- design$k; n <- design$n; K <- design$K
ds_start  <- as_int(opt$ds_start, "ds_start")
ds_end    <- as_int(opt$ds_end, "ds_end"); if (ds_end == 0L) ds_end <- design$n_ds
write_tsv <- as_bool(opt$tsv, "tsv")
do_check  <- as_bool(opt$check, "check")
check_only <- as_bool(opt$check_only, "check_only")
if (check_only) { ds_end <- ds_start; do_check <- TRUE }
overwrite <- as_bool(opt$overwrite, "overwrite")
thr <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
if (ds_start < 1 || ds_end > design$n_ds || ds_start > ds_end) stop("Invalid data-set range.")

tag   <- sprintf("%04d-%04d", ds_start, ds_end)
fdata <- file.path(opt$root, "data", sprintf("ds_%s.rds", tag))
dir.create(dirname(fdata), recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$root, "logs", sprintf("pw02_%s.log", tag)))
if (file.exists(fdata) && !overwrite && !check_only) {
  msg(fdata, " exists; use --overwrite=true to redo.")
  quit(save = "no")
}

scen_names <- design$scen$name
U  <- chol(design$R)                      # rows of W %*% U are MVN(0, R)
nr <- ds_end - ds_start + 1L
msg(sprintf("Data sets %d..%d: n = %d, %d traits, %d null variants each, %d scenarios",
            ds_start, ds_end, n, k, K, length(scen_names)))

SE     <- array(0, c(k, k, nr))
cvec   <- matrix(0, nr, k)
gg     <- numeric(nr); sdg <- numeric(nr); af <- numeric(nr)
z_null <- matrix(NA_real_, nr * K, k)
out    <- lapply(scen_names, function(s) list(t = matrix(0, nr, k), b = matrix(0, nr, k),
                                               se = matrix(0, nr, k)))
names(out) <- scen_names
t0 <- Sys.time()

for (ii in seq_len(nr)) {
  j <- ds_start + ii - 1L
  set.seed(design$seed_data + j)
  E  <- matrix(rnorm(n * k), n, k) %*% U
  E  <- sweep(E, 2, colMeans(E))                                  # centred traits
  g  <- rbinom(n, 2, design$maf_causal[j])
  G0 <- matrix(rbinom(n * K, 2, rep(design$maf_null[j, ], each = n)), n, K)

  gc_ <- g - mean(g)
  sdg[ii] <- sqrt(sum(gc_^2) / (n - 1))
  af[ii]  <- mean(g) / 2
  gs <- gc_ / sdg[ii]                                             # standardized, g'g = n - 1
  gg[ii]  <- sum(gs^2)
  SEj <- crossprod(E)
  cj  <- drop(crossprod(E, gs))
  SE[, , ii] <- SEj
  cvec[ii, ] <- cj
  yy0 <- diag(SEj)

  # Null variants against the pure-noise phenotypes.
  G0 <- sweep(G0, 2, colMeans(G0))
  s0 <- sqrt(colSums(G0^2) / (n - 1))
  ok <- s0 > 0
  G0 <- sweep(G0[, ok, drop = FALSE], 2, s0[ok], "/")
  b0 <- crossprod(G0, E) / (n - 1)
  t0n <- b0 / sqrt(sweep(-(b0^2) * (n - 1), 2, yy0, "+") / (n - 2) / (n - 1))
  z_null[(ii - 1L) * K + which(ok), ] <- t0n

  # The causal variant, in every scenario.
  for (s in scen_names) {
    beta <- design$beta[[s]][j, ]
    gy <- cj + gg[ii] * beta                                      # g'Y
    yy <- yy0 + 2 * beta * cj + gg[ii] * beta^2                   # diag(Y'Y)
    b  <- gy / gg[ii]
    se <- sqrt((yy - b^2 * gg[ii]) / (n - 2) / gg[ii])
    out[[s]]$t[ii, ]  <- b / se
    out[[s]]$b[ii, ]  <- b
    out[[s]]$se[ii, ] <- se
  }

  if (do_check && ii == 1L) {
    s    <- scen_names[length(scen_names) %/% 2 + 1L]
    beta <- design$beta[[s]][j, ]
    Y    <- E + outer(gs, beta)
    tr   <- c(1L, which.max(abs(beta)), k)
    t_lm <- vapply(tr, function(tt) coef(summary(lm(Y[, tt] ~ g)))[2, 3], 0)
    d_t  <- max(abs(t_lm - out[[s]]$t[ii, tr]))
    d_R  <- max(abs(cor(Y) - rhat_update(SEj, cj, gg[ii], beta)))
    msg(sprintf("  check on data set %d, %s: max |t(lm) - t| = %.2e, max |cor(Y) - Rhat| = %.2e",
                j, s, d_t, d_R))
  }

  if (ii %% 25L == 0L || ii == nr) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    msg(sprintf("  data set %d (%d/%d), %.1f min elapsed, ~%.1f min left",
                j, ii, nr, el, el / ii * (nr - ii)))
  }
}

# ------------------------------------------------------------------- output --

if (check_only) {
  msg("Check only: nothing written.")
  quit(save = "no")
}
ds <- ds_start:ds_end
for (s in scen_names) {
  sdir <- file.path(opt$root, "sumstats", s)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(scenario = s, ds = ds, t = out[[s]]$t, b = out[[s]]$b, se = out[[s]]$se,
               af = af, sdg = sdg), file.path(sdir, sprintf("ss_%s.rds", tag)))
  if (write_tsv) {
    # Long format; BETA and SE per allele in trait-SD units.
    long <- data.table(
      SNP   = rep(sprintf("ds%04d_causal", ds), times = k),
      TRAIT = rep(sprintf("QT%d", seq_len(k) - 1L), each = nr),
      N     = n,
      AF    = rep(af, times = k),
      BETA  = as.vector(out[[s]]$b / sdg),
      SE    = as.vector(out[[s]]$se / sdg),
      T     = as.vector(out[[s]]$t))
    long[, P := 2 * pt(-abs(T), df = n - 2)]
    fwrite(long, file.path(sdir, sprintf("ss_%s.tsv.gz", tag)), sep = "\t", compress = "gzip")
  }
}
saveRDS(list(ds = ds, SE = SE, cvec = cvec, gg = gg, sdg = sdg, af = af, z_null = z_null,
             null_ds = rep(ds, each = K), null_snp = rep(seq_len(K), times = nr)), fdata)

msg(sprintf("Wrote %s and summary statistics for %d scenarios (%.1f min).", fdata,
            length(scen_names), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
