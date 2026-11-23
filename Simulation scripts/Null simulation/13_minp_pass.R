#!/usr/bin/env Rscript
# 13_minp_pass.R
#
# Null calibration of the two minP tests, from the same stored z-scores and the
# same trait subsets as 12_trait_subsets.R (k = 16, 64, 128):
#
#   p_min                 = 2 * Phi(-max_i |z_i|) over the k traits
#   Bonferroni minP       p = min(1, k * p_min)        (folder suffix _minpB)
#   effective-tests minP  p = min(1, M_eff * p_min)    (folder suffix _minpE)
#
# M_eff is the Li & Ji (2005) estimator applied to the replicate's R block for
# those traits (cov2cor first, so the overlap-adjusted fill matrix works too),
# as one would estimate it from one's own data. For k = 128 it is about 77.
# Testing p < alpha is the same as testing p_min < alpha / k (or / M_eff).
#
# One pass over tstat.f32 per replicate writes, for each method and k, the files
# 07 and 11 read, in <inbase>/<out_scenario>_<method>_k<k>/rep###/:
#   direct_counts.txt, direct_hist.txt.gz                         (all variants)
#   direct_counts_common.txt, direct_hist_common.txt.gz,
#   direct_counts_rare.txt,   direct_hist_rare.txt.gz              (MAF classes)
#   qc_summary.txt                                                  (written last)
# There is no chi-square statistic, so no T histogram is written and 07 and 11
# report no lambda_GC for these folders.
#
# Usage:
#   Rscript 13_minp_pass.R --scenario=S1_normal --rep_start=1 --rep_end=257
#   Rscript 13_minp_pass.R --scenario=S4_normal_miss1pct --rhat_prefix=Rhat_adjfill \
#       --out_scenario=S4_adjfill --rep_start=1 --rep_end=257

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  tstat_root   = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current",
  inbase       = "/data5/bguo1/JAGWAS_null_input",
  af           = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current/shared/variant_af.tsv",
  scenario     = "S1_normal",
  rep_start    = "1",
  rep_end      = "257",
  k            = "128",         # traits in tstat.f32 and Rhat
  subset_k     = "16,64,128",   # numbers of traits to test
  subset_seed  = "20261008",    # must match 12_trait_subsets.R
  n_variant    = "8931083",
  layout       = "variant",
  maf_cut      = "0.01",
  chunk        = "500000",
  threads      = "1",
  overwrite    = "false",
  rhat_prefix  = "Rhat",        # only used for M_eff
  out_scenario = ""             # folder base name; default = --scenario
)
opt       <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
out_base  <- if (nzchar(opt$out_scenario)) opt$out_scenario else opt$scenario
k_all     <- as_int(opt$k, "k")
ks        <- sort(unique(as.integer(trimws(strsplit(opt$subset_k, ",", fixed = TRUE)[[1]]))))
seed_base <- as_int(opt$subset_seed, "subset_seed")
nv        <- as_num(opt$n_variant, "n_variant")
chunk     <- as_int(opt$chunk, "chunk")
maf_cut   <- as_num(opt$maf_cut, "maf_cut")
rep_start <- as_int(opt$rep_start, "rep_start")
rep_end   <- as_int(opt$rep_end, "rep_end")
overwrite <- as_bool(opt$overwrite, "overwrite")
thr       <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
layout_variant <- identical(opt$layout, "variant")
if (anyNA(ks) || any(ks < 1 | ks > k_all)) stop("--subset_k must be integers in 1..", k_all)
if (!opt$scenario %in% ALL_SCENARIOS) stop("Unknown scenario '", opt$scenario, "'")

# Same thresholds and -log10 p histogram as 05, 10 and 12.
THR <- c(0.5, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8, 1e-9, 1e-10)
LNB <- 10001L; LBW <- 0.01
CLASSES <- c("common", "rare")
METHODS <- c("minpB", "minpE")
lthr <- -log10(THR)
i58  <- which(THR == 5e-8)

log_open(file.path(opt$inbase, "logs",
                   sprintf("13_minp_%s_rep%03d-%03d.log", out_base, rep_start, rep_end)))
on.exit(log_close(), add = TRUE)
msg("minP null calibration, k = ", paste(ks, collapse = ", "), ": ", opt$scenario,
    " replicates ", rep_start, "..", rep_end, "; results in ",
    file.path(opt$inbase, paste0(out_base, "_minp[BE]_k<k>")))

af <- fread(opt$af, select = c("variant_row", "minor_allele_frequency"), showProgress = FALSE)
if (nrow(af) != nv) stop(sprintf("%s has %d rows; the binary has %.0f variants.", opt$af, nrow(af), nv))
setorder(af, variant_row)
if (af$variant_row[1] != 1 || af$variant_row[nrow(af)] != nv)
  stop("variant_row is not 1..n_variant; cannot align the frequency file to the binary.")
cls <- ifelse(af$minor_allele_frequency > maf_cut, 1L, 2L)
rm(af)

read_chunk_t <- function(con, a, len) {
  if (layout_variant) {
    seek(con, where = (a - 1) * k_all * 4, origin = "start")
    matrix(readBin(con, "numeric", n = len * k_all, size = 4, endian = "little"),
           nrow = k_all, ncol = len)
  } else {
    Zt <- matrix(0, k_all, len)
    for (t in seq_len(k_all)) {
      seek(con, where = ((t - 1) * nv + (a - 1)) * 4, origin = "start")
      Zt[t, ] <- readBin(con, "numeric", n = len, size = 4, endian = "little")
    }
    Zt
  }
}
write_hist <- function(lh, f) {
  nz <- which(lh > 0)
  fwrite(data.table(type = "lp", bin = nz, count = lh[nz]), f, sep = "\t", compress = "gzip")
}
write_counts <- function(cnt, n, f) {
  e <- n * THR
  fwrite(data.table(threshold = THR, observed = cnt, expected = e, ratio = cnt / e), f, sep = "\t")
}

for (r in rep_start:rep_end) {
  tag  <- sprintf("rep%03d", r)
  fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
  frh  <- file.path(opt$inbase, opt$scenario, sprintf("%s_rep%03d.txt", opt$rhat_prefix, r))
  rdir <- function(m, kk) file.path(opt$inbase, sprintf("%s_%s_k%d", out_base, m, kk), tag)
  done <- unlist(lapply(METHODS, function(m) file.path(rdir(m, ks), "qc_summary.txt")))
  if (!overwrite && all(file.exists(done))) { msg(tag, ": done, skipping."); next }
  if (!file.exists(fbin)) { msg("MISSING ", fbin); next }
  if (!file.exists(frh))  { msg("MISSING ", frh);  next }
  t0 <- Sys.time()

  Rhat <- as.matrix(fread(frh, header = FALSE, showProgress = FALSE))
  dimnames(Rhat) <- NULL
  sels <- lapply(ks, function(kk) trait_subset(k_all, kk, seed_base, r))
  for (i in seq_along(ks)[-1]) stopifnot(all(sels[[i - 1]] %in% sels[[i]]))   # nested
  meff <- vapply(sels, function(s) li_ji_meff(cov2cor(Rhat[s, s, drop = FALSE])), 0)
  mult <- list(minpB = ks, minpE = meff)
  # traits added at each step, so max |z| is built up once over the nested subsets
  add <- lapply(seq_along(ks), function(i) if (i == 1) sels[[1]] else setdiff(sels[[i]], sels[[i - 1]]))

  acc <- lapply(METHODS, function(m) lapply(ks, function(kk) list(
    cnt = numeric(length(THR)), lh = numeric(LNB),
    ccnt = matrix(0, 2, length(THR)), clh = matrix(0, 2, LNB), nc = numeric(2))))
  names(acc) <- METHODS
  n_ok <- 0

  con <- file(fbin, "rb")
  for (a in seq(1, nv, by = chunk)) {
    b <- min(nv, a + chunk - 1); len <- as.integer(b - a + 1)
    Zt <- read_chunk_t(con, a, len)
    ok <- is.finite(colSums(Zt))                       # all 128 traits, as in 05 and 12
    A  <- abs(Zt[, ok, drop = FALSE]); rm(Zt)
    ci <- cls[(a:b)[ok]]
    n_ok <- n_ok + ncol(A)
    m <- numeric(ncol(A))
    for (i in seq_along(ks)) {
      for (t in add[[i]]) m <- pmax(m, A[t, ])
      lpmin <- log(2) + pnorm(-m, log.p = TRUE)          # natural log of p_min
      for (meth in METHODS) {
        lp <- pmax(0, -(log(mult[[meth]][i]) + lpmin) / log(10))   # -log10 of min(1, M p_min)
        lb <- as.integer(pmin(lp / LBW, LNB - 1)) + 1L
        s  <- acc[[meth]][[i]]
        s$cnt <- s$cnt + vapply(lthr, function(x) sum(lp >= x), 0)
        s$lh  <- s$lh + tabulate(lb, LNB)
        for (cl in 1:2) {
          w <- ci == cl; lpw <- lp[w]
          s$nc[cl] <- s$nc[cl] + sum(w)
          s$ccnt[cl, ] <- s$ccnt[cl, ] + vapply(lthr, function(x) sum(lpw >= x), 0)
          s$clh[cl, ]  <- s$clh[cl, ] + tabulate(lb[w], LNB)
        }
        acc[[meth]][[i]] <- s
      }
    }
  }
  close(con)

  for (meth in METHODS) for (i in seq_along(ks)) {
    s <- acc[[meth]][[i]]; kk <- ks[i]; d <- rdir(meth, kk)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(file.path(d, "qc_summary.txt"))) file.remove(file.path(d, "qc_summary.txt"))
    write_hist(s$lh, file.path(d, "direct_hist.txt.gz"))
    write_counts(s$cnt, n_ok, file.path(d, "direct_counts.txt"))
    for (cl in 1:2) {
      write_hist(s$clh[cl, ], file.path(d, sprintf("direct_hist_%s.txt.gz", CLASSES[cl])))
      write_counts(s$ccnt[cl, ], s$nc[cl], file.path(d, sprintf("direct_counts_%s.txt", CLASSES[cl])))
    }
    summ <- data.table(
      field = c("scenario", "rhat_prefix", "replicate", "method", "k_tested", "multiplier",
                "li_ji_meff", "traits", "subset_seed", "n_variant_written",
                "direct_obs_5e-8", "direct_exp_5e-8", "direct_ratio_5e-8", "minutes"),
      value = c(opt$scenario, opt$rhat_prefix, r, meth, kk, mult[[meth]][i], meff[i],
                paste(sels[[i]], collapse = ","), seed_base + r, format(n_ok, scientific = FALSE),
                s$cnt[i58], n_ok * THR[i58], s$cnt[i58] / (n_ok * THR[i58]),
                round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)))
    setnames(summ, "field", "key")
    fwrite(summ, file.path(d, "qc_summary.txt"), sep = "\t")      # last: marks completion
  }
  for (i in seq_along(ks))
    msg(sprintf("%s k=%3d M_eff %.1f | P < 5e-8: Bonferroni %.0f, effective-tests %.0f, expected %.2f",
                tag, ks[i], meff[i], acc$minpB[[i]]$cnt[i58], acc$minpE[[i]]$cnt[i58], n_ok * 5e-8))
  msg(sprintf("%s: %.1f min", tag, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
msg("Done.")
