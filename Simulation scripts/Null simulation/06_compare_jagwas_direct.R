#!/usr/bin/env Rscript
# 06_compare_jagwas_direct.R
#
# For the replicates run through the real JAGWAS binary, checks that its output
# matches the direct z' Rhat^-1 z computed by 05. Passing this is what justifies
# using the direct statistic for every other replicate.
#
# JAGWAS writes log10_P (log10 p, so negative) to 6 significant digits, so the
# two agree to within ~5e-6 relative, not exactly. JAGWAS's output rows are the rows of
# its input files -- the variants 05 kept, in order -- which is also the order
# of direct_neglog10p.f32.
#
# Needs 05 run with --save_direct_p=true; run_jagwas_null.sh does that by default.
#
# Usage:
#   Rscript 06_compare_jagwas_direct.R --scenarios=all --reps=1-5

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  inbase    = "/data5/bguo1/JAGWAS_null_input",
  outbase   = "/data5/bguo1/JAGWAS_null_output",
  scenarios = "all",
  reps      = "1-5",       # "1-5", "1,2,7", or "all"
  tol       = "1e-4"       # max relative difference allowed in -log10 p
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
tol <- as_num(opt$tol, "tol")

scen <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
reps <- if (identical(opt$reps, "all")) 1:257 else
  unlist(lapply(strsplit(opt$reps, ",", fixed = TRUE)[[1]], function(x) {
    ab <- as.integer(strsplit(trimws(x), "-", fixed = TRUE)[[1]])
    if (length(ab) == 2) ab[1]:ab[2] else ab
  }))

log_open(file.path(opt$outbase, "logs", "06_compare_jagwas_direct.log"))
on.exit(log_close(), add = TRUE)

cut58 <- -log10(5e-8)
res <- vector("list", 0L)

for (s in scen) for (r in reps) {
  tag <- sprintf("rep%03d", r)
  fj  <- file.path(opt$outbase, s, paste0(tag, ".jagwas.txt.gz"))
  fd  <- file.path(opt$inbase, s, tag, "direct_neglog10p.f32")
  if (!file.exists(fj)) next
  if (!file.exists(fd)) {
    msg(sprintf("%s %s: JAGWAS output present but no direct_neglog10p.f32 (05 run without --save_direct_p)", s, tag))
    next
  }

  j <- fread(cmd = paste("zcat", shQuote(fj)), select = "log10_P", showProgress = FALSE)
  x <- -j$log10_P
  n <- length(x)
  d <- readBin(fd, what = "numeric", n = n + 1L, size = 4)
  if (length(d) != n) {
    msg(sprintf("%s %s: FAIL -- JAGWAS has %d rows, direct has %d", s, tag, n, length(d)))
    res[[length(res) + 1L]] <- data.table(scenario = s, replicate = r, n_jagwas = n,
                                          n_direct = length(d), verdict = "FAIL: row count")
    next
  }

  rel <- abs(x - d) / pmax(1, d)
  i   <- which.max(rel)
  cj  <- sum(x >= cut58); cd <- sum(d >= cut58)
  ok  <- max(rel) < tol && abs(cj - cd) <= 1   # a value within rounding of 7.30103 may flip
  res[[length(res) + 1L]] <- data.table(
    scenario = s, replicate = r, n_jagwas = n, n_direct = length(d),
    max_rel_diff = max(rel), row_of_max = i, jagwas_at_max = x[i], direct_at_max = d[i],
    cor = cor(x, d), jagwas_hits_5e8 = cj, direct_hits_5e8 = cd,
    verdict = if (ok) "PASS" else "FAIL")
  msg(sprintf("%-20s %s  %s  max rel diff %.2e  cor %.8f  hits at 5e-8: JAGWAS %d / direct %d",
              s, tag, if (ok) "PASS" else "FAIL", max(rel), cor(x, d), cj, cd))
}

if (length(res)) {
  out <- rbindlist(res, fill = TRUE)
  dir.create(file.path(opt$outbase, "qc"), recursive = TRUE, showWarnings = FALSE)
  f <- file.path(opt$outbase, "qc", "jagwas_vs_direct.txt")
  fwrite(out, f, sep = "\t")
  msg("")
  msg(sprintf("%d of %d replicates PASS. Wrote %s", sum(out$verdict == "PASS"), nrow(out), f))
} else {
  msg("No replicate had both a JAGWAS output and a direct_neglog10p.f32.")
}
