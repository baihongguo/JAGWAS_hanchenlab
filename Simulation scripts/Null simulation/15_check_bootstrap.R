#!/usr/bin/env Rscript
# 15_check_bootstrap.R
#
# Checks the replicate-bootstrap intervals of 07 against quantities computed
# directly from the per-replicate counts, without resampling:
#
#   1. Spread. For each scenario and alpha: the standard deviation, across the
#      replicates, of the number of tests with P < alpha, next to the binomial
#      standard deviation sqrt(n p (1 - p)) that independent tests would give
#      (n tests per replicate, p the pooled rate). The interval
#      p +/- 1.96 SD(rate) / sqrt(R), the usual Monte Carlo standard error over
#      independent replicates, is printed next to the bootstrap interval
#      (recomputed exactly as in 07: same seed, same resamples) and the exact
#      Poisson interval, which treats every test as independent.
#   2. Clustering at 5e-8. From direct_hits.txt: how many replicates have hits,
#      and how many regions they form (hits of one replicate fewer than --gap
#      variants apart are one region; variant indices follow genome order).
#
# Reads only direct_counts.txt and direct_hits.txt in each replicate folder
# and writes nothing but its log.
#
# Usage:
#   Rscript 15_check_bootstrap.R
#   Rscript 15_check_bootstrap.R --scenarios=S1_normal,S4_adjfill --gap=3000

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  inbase    = "/data5/bguo1/JAGWAS_null_input",
  scenarios = "S1_normal,S1_normal_k64,S1_normal_k16",
  alphas    = "0.05,1e-3,1e-5,5e-8",
  B         = "10000",      # as in 07
  seed      = "20260919",   # as in 07
  gap       = "3000"        # variants; about 1 Mb at 8.93 million variants
)
opt    <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
B      <- as_int(opt$B, "B")
seed   <- as_int(opt$seed, "seed")
gap    <- as_num(opt$gap, "gap")
alphas <- as.numeric(strsplit(opt$alphas, ",", fixed = TRUE)[[1]])
scen   <- trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])

log_open(file.path(opt$inbase, "logs", "15_check_bootstrap.log"))
on.exit(log_close(), add = TRUE)

g4  <- function(x) formatC(x, digits = 4, format = "g")   # 4 digits: at 0.05 the intervals are narrow
ci  <- function(lo, hi) sprintf("%s to %s", g4(lo), g4(hi))

for (s in scen) {
  sdir <- file.path(opt$inbase, s)
  reps <- sort(list.files(sdir, pattern = "^rep[0-9]{3}$"))
  # the replicates 07 uses: finished (qc_summary.txt) and with a histogram
  reps <- reps[file.exists(file.path(sdir, reps, "qc_summary.txt")) &
               file.exists(file.path(sdir, reps, "direct_hist.txt.gz"))]
  R <- length(reps)
  if (!R) { msg(s, ": no completed replicates, skipping."); next }

  cnt <- rbindlist(lapply(reps, function(rp) {
    d <- fread(file.path(sdir, rp, "direct_counts.txt")); d[, replicate := rp]; d }))
  cnt[, `:=`(observed = as.numeric(observed), expected = as.numeric(expected))]
  setorder(cnt, threshold, replicate)                 # the order 07 resamples in
  set.seed(seed)
  idx <- matrix(sample.int(R, R * B, replace = TRUE), nrow = B)

  msg("")
  msg(sprintf("%s: %d replicates", s, R))
  msg(sprintf("   %-6s %11s %10s %10s %7s   %-24s %-24s %-24s", "alpha", "mean count", "SD", "binom SD",
              "SD/bin", "replicate SE, +/-1.96", "bootstrap (as 07)", "exact Poisson"))
  for (a in alphas) {
    d <- cnt[abs(threshold - a) < a * 1e-6]
    if (nrow(d) != R) { msg("   alpha ", a, ": not in every replicate, skipped"); next }
    o <- d$observed
    n <- d$expected / a                               # tests in each replicate
    p <- sum(o) / sum(n)
    sd_bin <- sqrt(mean(n) * p * (1 - p))
    se_rep <- sd(o / n) / sqrt(R)
    boot <- rowSums(matrix(o[idx], nrow = B)) / rowSums(matrix(n[idx], nrow = B))
    bq <- quantile(boot, c(0.025, 0.975), names = FALSE)
    O <- sum(o); N <- sum(n)
    pl <- if (O == 0) 0 else qchisq(0.025, 2 * O) / 2 / N
    ph <- qchisq(0.975, 2 * (O + 1)) / 2 / N
    msg(sprintf("   %-6s %11.1f %10.1f %10.1f %7.2f   %-24s %-24s %-24s", g4(a), mean(o), sd(o), sd_bin,
                sd(o) / sd_bin, ci(p - 1.96 * se_rep, p + 1.96 * se_rep), ci(bq[1], bq[2]), ci(pl, ph)))
  }

  # clustering of the hits at 5e-8
  o58 <- cnt[abs(threshold - 5e-8) < 5e-14, observed]
  hits <- rbindlist(lapply(seq_along(reps), function(i) {
    f <- file.path(sdir, reps[i], "direct_hits.txt")
    if (!file.exists(f)) return(NULL)
    h <- fread(f, select = "index"); h[, rep := i]; h }))
  if (sum(o58) == 0) { msg("   P < 5e-8: no hits"); next }
  if (nrow(hits) != sum(o58))
    msg(sprintf("   ** direct_hits.txt has %d hits, direct_counts.txt %.0f: check --hit_p of 05/12 **",
                nrow(hits), sum(o58)))
  setorder(hits, rep, index)
  hits[, region := cumsum(c(TRUE, diff(index) > gap | diff(rep) != 0))]
  reg <- hits[, .N, by = region]
  msg(sprintf("   P < 5e-8: %d hits in %d of %d replicates, forming %d regions (hits < %s variants apart);",
              nrow(hits), uniqueN(hits$rep), R, nrow(reg), format(gap, big.mark = ",")))
  msg(sprintf("             %d regions with one hit, the largest with %d; hits per replicate: %s",
              sum(reg$N == 1), max(reg$N),
              paste(sprintf("%s x %d", names(table(o58)), as.integer(table(o58))), collapse = ", ")))
}
msg("")
msg("SD/bin = 1 when tests are independent; the bootstrap interval should match the replicate-SE one.")
