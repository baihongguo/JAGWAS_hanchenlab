#!/usr/bin/env Rscript
# 14_summarise_maf.R
#
# Pools the per-replicate MAF-bin counts written by 09_maf_and_univariate.R
# into Supplementary Table 4: the empirical type I error rate at --threshold
# (default 5e-8) in each MAF bin, for the JAGWAS statistic of the 128 traits and
# for the single-trait statistics, per scenario.
#
#   <inbase>/results/maf_table_<threshold>.txt   one row per scenario, statistic
#                                                and MAF bin, with a replicate-
#                                                bootstrap 95% CI (as in 07) and
#                                                how many replicates the hits
#                                                come from (n_rep_hit, max_rep)
#
# 09 is run in chunks of replicates, each with its own --out_tag (for example
# _reps1-65), so the chunk files do not overwrite each other. Every file whose
# name starts with maf_multivariate_<scenario><--tag_prefix> or
# maf_univariate_<scenario><--tag_prefix> is read; a replicate that appears in
# two chunks is counted once. The log lists the replicates found per scenario
# and stops if any of 1..--n_rep is missing.
#
# CHECK. Summed over the bins, the JAGWAS hits must equal those of 05 for the
# same replicates, which 07 pools into type1_table_all.txt; the log compares the
# two when that table exists.
#
# Usage:
#   Rscript 14_summarise_maf.R
#   Rscript 14_summarise_maf.R --tag_prefix=_reps --n_rep=257 --threshold=5e-8

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  inbase     = "/data5/bguo1/JAGWAS_null_input",
  scenarios  = "all",
  k          = "128",
  n_rep      = "257",
  threshold  = "5e-8",
  tag_prefix = "_reps",
  B          = "10000",      # bootstrap resamples of whole replicates, as in 07
  seed       = "20260919"
)
opt   <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k     <- as_int(opt$k, "k")
n_rep <- as_int(opt$n_rep, "n_rep")
thr   <- as_num(opt$threshold, "threshold")
B     <- as_int(opt$B, "B")
set.seed(as_int(opt$seed, "seed"))
idx   <- matrix(sample.int(n_rep, n_rep * B, replace = TRUE), nrow = B)   # the same resamples for every bin
scen  <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])

resdir <- file.path(opt$inbase, "results")
log_open(file.path(opt$inbase, "logs", sprintf("14_summarise_maf_%s.log", opt$threshold)))
on.exit(log_close(), add = TRUE)

read_chunks <- function(kind, s) {
  pat <- sprintf("^maf_%s_%s%s.*\\.txt$", kind, s, gsub("([.+-])", "\\\\\\1", opt$tag_prefix))
  fs  <- list.files(resdir, pattern = pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  d <- rbindlist(lapply(fs, fread))
  d[, `:=`(observed = as.numeric(observed), n_variants = as.numeric(n_variants))]
  dup <- d[, .N, by = .(replicate, maf_bin, threshold)][N > 1]
  if (nrow(dup)) msg(sprintf("   %s %s: %d replicates appear in more than one chunk; counted once",
                             s, kind, uniqueN(dup$replicate)))
  unique(d, by = c("replicate", "maf_bin", "threshold"))
}

bin_lower <- function(b) as.numeric(sub("^[\\[(]([^,]+),.*$", "\\1", b))

fmt_rate <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  if (x >= 0.01) return(sprintf("%.3g", x))
  e <- floor(log10(x)); m <- x / 10^e
  if (round(m, 2) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%.3g x 10^%d", m, e)
}

out <- list()
for (s in scen) {
  for (kind in c("multivariate", "univariate")) {
    d <- read_chunks(kind, s)
    if (is.null(d)) { msg(s, " ", kind, ": no chunk files, skipping."); next }
    reps <- sort(unique(d$replicate))
    miss <- setdiff(seq_len(n_rep), reps)
    msg(sprintf("%s %s: %d replicates (%d-%d)", s, kind, length(reps), min(reps), max(reps)))
    if (length(miss))
      stop(sprintf("%s %s: %d replicates missing, e.g. %s", s, kind, length(miss),
                   paste(head(miss, 10), collapse = ", ")), call. = FALSE)
    d <- d[abs(threshold - thr) < thr * 1e-6]
    if (!nrow(d)) stop("Threshold ", opt$threshold, " not in the ", kind, " files of ", s, call. = FALSE)
    per_test <- if (kind == "multivariate") 1 else k   # 09 counts single-trait hits over all k traits
    setorder(d, maf_bin, replicate)
    tab <- d[, {
      o <- observed; n <- n_variants * per_test       # one value per replicate, replicates 1..n_rep
      boot <- rowSums(matrix(o[idx], nrow = B)) / rowSums(matrix(n[idx], nrow = B))
      list(n_rep = uniqueN(replicate), tests = sum(n), observed = sum(o),
           ci_boot_lo = unname(quantile(boot, 0.025)), ci_boot_hi = unname(quantile(boot, 0.975)),
           n_rep_hit = sum(o > 0), max_rep = max(o))    # how concentrated the hits are
    }, by = maf_bin]
    tab[, `:=`(scenario = s, statistic = if (kind == "multivariate") "JAGWAS" else "single-trait",
               rate = observed / tests)]
    out[[length(out) + 1L]] <- tab
  }
}
if (!length(out)) stop("Nothing to summarise.")

res <- rbindlist(out)
res[, lower := bin_lower(maf_bin)]
setorder(res, statistic, scenario, lower)
res[, lower := NULL]
setcolorder(res, c("scenario", "statistic", "maf_bin", "n_rep", "tests", "observed", "rate",
                   "ci_boot_lo", "ci_boot_hi", "n_rep_hit", "max_rep"))
f1 <- file.path(resdir, sprintf("maf_table_%s.txt", opt$threshold))
fwrite(res, f1, sep = "\t")

# JAGWAS hits summed over the bins against 05's counts for the same replicates.
f07 <- file.path(resdir, "type1_table_all.txt")
if (file.exists(f07)) {
  # 07's counts pass the 32-bit integer range, so read them as doubles (integer64 prints as 0 in sprintf)
  t07 <- fread(f07, integer64 = "double")[abs(threshold - thr) < thr * 1e-6]
  for (s in intersect(scen, t07$scenario)) {
    o09 <- res[scenario == s & statistic == "JAGWAS", sum(observed)]
    o07 <- as.numeric(t07[scenario == s, observed])
    if (length(o07) == 1)
      msg(sprintf("CHECK %s: %.0f JAGWAS hits summed over the MAF bins, %.0f in type1_table_all.txt%s",
                  s, o09, o07, if (isTRUE(abs(o09 - o07) <= 1)) "" else "  ** DIFFERENT **"))
  }
}

msg("")
msg(sprintf("Empirical type I error rate at P < %s, pooled over the replicates:", opt$threshold))
for (st in unique(res$statistic)) {
  msg(st)
  w <- dcast(res[statistic == st], maf_bin ~ scenario, value.var = "rate")
  tests <- res[statistic == st, .(tests = tests[1]), by = maf_bin]
  w <- merge(w, tests, by = "maf_bin")
  w <- w[order(bin_lower(maf_bin))]
  msg(sprintf("   %-15s %16s %s", "MAF bin", "tests",
              paste(sprintf("%-20s", setdiff(names(w), c("maf_bin", "tests"))), collapse = "")))
  for (i in seq_len(nrow(w))) {
    rates <- vapply(setdiff(names(w), c("maf_bin", "tests")), function(cn) fmt_rate(w[[cn]][i]), "")
    msg(sprintf("   %-15s %16s %s", w$maf_bin[i], format(w$tests[i], big.mark = ",", scientific = FALSE),
                paste(sprintf("%-20s", rates), collapse = "")))
  }
}
msg("")
msg("With 95% bootstrap intervals (whole replicates resampled), and how many replicates the hits come from:")
for (i in seq_len(nrow(res))) with(res[i], msg(sprintf("   %-12s %-18s %-15s %s (%s to %s); %.0f hits in %d replicates, at most %.0f in one",
                                                       statistic, scenario, maf_bin, fmt_rate(rate),
                                                       fmt_rate(ci_boot_lo), fmt_rate(ci_boot_hi),
                                                       observed, n_rep_hit, max_rep)))
msg("")
msg("Wrote ", f1)
