#!/usr/bin/env Rscript
# 07_summarise_type1.R
#
# Pools every replicate into the tables for the response letter:
#
#   <inbase>/results/type1_table_<variants>.txt        observed vs expected JAGWAS
#                                                      hits at each threshold, with
#                                                      95% CIs, per scenario
#   <inbase>/results/calibration_table_<variants>.txt  JAGWAS lambda_GC, plus
#                                                      univariate lambda_GC and E[T]/k
#                                                      for --variants=all
#
#   --variants=all     every variant: 05's direct_counts.txt, direct_hist.txt.gz
#   --variants=common  MAF > 1%:      10's direct_counts_common.txt, direct_hist_common.txt.gz
#   --variants=rare    MAF <= 1%:     10's *_rare files
#
# The QQ plots are drawn by 11_plot_qq.R.
#
# CONFIDENCE INTERVALS. Two are given for each observed/expected ratio:
#   * Poisson exact -- treats every p-value as independent. Too narrow here:
#     SNPs are in LD, so hits cluster within a replicate. For reference only.
#   * Replicate bootstrap -- resamples whole replicates, which are independent
#     by construction, so the LD clustering inside each replicate is carried
#     along automatically. This is the interval to report.
#
# TRAIT SUBSETS. Folders from 12_trait_subsets.R end in _k<number>; that number
# is used as the degrees of freedom for those folders, and --k for the others.
# minP folders from 13_minp_pass.R (..._minpB_k<k>, ..._minpE_k<k>) work the same
# way; they have no chi-square statistic, so their lambda_GC columns are NA.
#
# Usage:
#   Rscript 07_summarise_type1.R --variants=all
#   Rscript 07_summarise_type1.R --variants=common
#   Rscript 07_summarise_type1.R --variants=common --tag=k16 \
#       --scenarios=S1_normal_k16,S2_lognormal_k16,S3_cauchy_k16,S4_adjfill_k16

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  inbase    = "/data5/bguo1/JAGWAS_null_input",
  scenarios = "all",
  variants  = "all",        # all | common | rare
  k         = "128",
  B         = "10000",      # bootstrap resamples of replicates
  seed      = "20260919",
  tag       = ""            # appended to the output names; set automatically to
                            # "subset" when --scenarios is not "all", so a
                            # partial run never overwrites the main tables
)
opt  <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k    <- as_int(opt$k, "k")
B    <- as_int(opt$B, "B")
seed <- as_int(opt$seed, "seed")
scen <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
if (!opt$variants %in% c("all", "common", "rare"))
  stop("--variants must be all, common or rare.")
all_var  <- opt$variants == "all"
sfx      <- if (all_var) "" else paste0("_", opt$variants)
f_counts <- paste0("direct_counts", sfx, ".txt")
f_hist   <- paste0("direct_hist", sfx, ".txt.gz")
tag <- opt$tag
if (!nzchar(tag) && !identical(opt$scenarios, "all")) tag <- "subset"
sfx_out <- paste0(opt$variants, if (nzchar(tag)) paste0("_", tag) else "")

resdir <- file.path(opt$inbase, "results")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$inbase, "logs", sprintf("07_summarise_type1_%s.log", sfx_out)))
on.exit(log_close(), add = TRUE)
msg("Variants: ", opt$variants)

# Histogram layout -- must match 05 and 10.
TNB <- 40001L; TBW <- 0.01   # T on [0, 400)
LNB <- 10001L; LBW <- 0.01   # -log10 p on [0, 100)

read_kv <- function(f) { d <- fread(f, sep = "\t", colClasses = "character"); setNames(d$value, d$key) }
read_gz <- function(f) fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE)

# lambda_GC of the df-df statistic from a histogram of T.
lambda_T <- function(bins, counts, df) {
  th <- numeric(TNB); th[bins] <- counts
  if (sum(th) == 0) return(NA_real_)                 # minP folders (13) have no T histogram
  b <- which(cumsum(th) >= sum(th) / 2)[1]
  if (is.na(b) || b == TNB) NA_real_ else (b - 0.5) * TBW / qchisq(0.5, df = df)
}

type1 <- list(); calib <- list()

for (s in scen) {
  sdir <- file.path(opt$inbase, s)
  reps <- sort(list.files(sdir, pattern = "^rep[0-9]{3}$"))
  done <- file.path(sdir, reps, if (all_var) "qc_summary.txt" else f_counts)
  reps <- reps[file.exists(done) & file.exists(file.path(sdir, reps, f_hist))]
  if (!length(reps)) { msg(s, ": no completed replicates, skipping."); next }
  R <- length(reps)
  ks <- scenario_k(s, k)
  msg(sprintf("%s: %d completed replicates, k = %d", s, R, ks))

  # --- per-replicate inputs ---------------------------------------------------
  cnt <- rbindlist(lapply(reps, function(rp) {
    d <- fread(file.path(sdir, rp, f_counts)); d[, replicate := rp]; d }))
  setorder(cnt, threshold, replicate)
  # Counts summed over all replicates can pass the integer limit (S4 at p < 0.05
  # is ~2.2e9), so work in doubles throughout.
  cnt[, `:=`(observed = as.numeric(observed), expected = as.numeric(expected))]

  hl <- lapply(reps, function(rp) read_gz(file.path(sdir, rp, f_hist)))
  lgc_rep <- vapply(hl, function(h) { ht <- h[type == "T"]; lambda_T(ht$bin, ht$count, ks) }, 0)
  hs <- rbindlist(hl)[, .(count = sum(count)), by = .(type, bin)]
  rm(hl)
  N <- hs[type == "lp", sum(count)]
  lgc_pool <- lambda_T(hs[type == "T", bin], hs[type == "T", count], ks)

  # --- type I error: pooled counts, Poisson and replicate-bootstrap CIs -------
  set.seed(seed)
  idx <- matrix(sample.int(R, R * B, replace = TRUE), nrow = B)
  tab <- cnt[, {
    O <- sum(observed); E <- sum(expected)
    boot <- if (.N == R) {
      rowSums(matrix(observed[idx], nrow = B)) / rowSums(matrix(expected[idx], nrow = B))
    } else NA_real_                                  # threshold missing in some replicate
    list(n_rep = .N, observed = O, expected = E, ratio = O / E,
         ci_poisson_lo = if (O == 0) 0 else qchisq(0.025, 2 * O) / 2 / E,
         ci_poisson_hi = qchisq(0.975, 2 * (O + 1)) / 2 / E,
         ci_boot_lo = unname(quantile(boot, 0.025, na.rm = TRUE)),
         ci_boot_hi = unname(quantile(boot, 0.975, na.rm = TRUE)))
  }, by = threshold]
  tab[, `:=`(scenario = s, k = ks, variants = opt$variants)]
  type1[[s]] <- tab

  # --- calibration ------------------------------------------------------------
  # Univariate lambda_GC and E[T]/k are recorded by 05 over all variants only.
  if (all_var) {
    qs <- rbindlist(lapply(reps, function(rp) {
      kv <- read_kv(file.path(sdir, rp, "qc_summary.txt"))
      data.table(E_T_over_k = as.numeric(kv["E_T_over_k"]),
                 pred_ratio = as.numeric(kv["predicted_type1_ratio_5e-8"]))
    }))
    lgc_uni <- unlist(lapply(reps, function(rp) {
      f <- file.path(sdir, rp, "qc_traits.txt")       # absent in minP folders (13)
      if (file.exists(f)) fread(f)$lambda_gc else NA_real_
    }))
  }
  calib[[s]] <- data.table(
    scenario = s, k = ks, variants = opt$variants, n_rep = R, n_pvalues = N,
    jagwas_lambda_gc_pooled   = lgc_pool,
    jagwas_lambda_gc_rep_mean = mean(lgc_rep, na.rm = TRUE),
    jagwas_lambda_gc_rep_sd   = sd(lgc_rep, na.rm = TRUE),
    univariate_lambda_gc_median = if (all_var) median(lgc_uni, na.rm = TRUE) else NA_real_,
    univariate_lambda_gc_q025   = if (all_var) unname(quantile(lgc_uni, 0.025, na.rm = TRUE)) else NA_real_,
    univariate_lambda_gc_q975   = if (all_var) unname(quantile(lgc_uni, 0.975, na.rm = TRUE)) else NA_real_,
    E_T_over_k_mean          = if (all_var) mean(qs$E_T_over_k) else NA_real_,
    predicted_ratio_5e8_mean = if (all_var) mean(qs$pred_ratio) else NA_real_)

  r58 <- tab[threshold == 5e-8]
  msg(sprintf("   %.3g p-values | JAGWAS lambda_GC %.4f | P < 5e-8: %d observed, %.1f expected, ratio %.3f (bootstrap 95%% CI %.3f-%.3f)",
              N, lgc_pool, as.integer(r58$observed), r58$expected, r58$ratio,
              r58$ci_boot_lo, r58$ci_boot_hi))
}

if (!length(type1)) stop("Nothing to summarise.")

t1 <- rbindlist(type1)
setcolorder(t1, c("scenario", "k", "variants", "threshold"))
f1 <- file.path(resdir, sprintf("type1_table_%s.txt", sfx_out))
f2 <- file.path(resdir, sprintf("calibration_table_%s.txt", sfx_out))
fwrite(t1, f1, sep = "\t")
fwrite(rbindlist(calib), f2, sep = "\t")

msg("")
msg("Wrote ", f1)
msg("Wrote ", f2)
msg("Report the replicate-bootstrap interval; the Poisson one ignores LD clustering.")
