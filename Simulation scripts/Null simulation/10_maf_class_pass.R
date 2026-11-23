#!/usr/bin/env Rscript
# 10_maf_class_pass.R
#
# Per replicate: the JAGWAS statistic for every variant, split into common and
# rare variants, saved in the same formats as 05 so that 07 (tables) and 11
# (QQ plots) can run on either class:
#
#   <inbase>/<scenario>/rep###/direct_hist_common.txt.gz, direct_counts_common.txt
#   <inbase>/<scenario>/rep###/direct_hist_rare.txt.gz,   direct_counts_rare.txt
#
#   common = MAF > --maf_cut (default 0.01);  rare = MAF <= --maf_cut
#
# THE STATISTIC is 05's, T = z' Rhat^-1 z, computed as ||L^-1 z||^2 with
# Rhat = L L'. A triangular solve needs half the arithmetic of multiplying by
# Rhat^-1, and it works on each chunk exactly as stored in the variant-major
# binary (one column per variant), so nothing is transposed.
#
# CHECK AGAINST 05. Common plus rare must reproduce 05's direct_counts.txt for
# the same replicate at every threshold. The log reports it per replicate; 05
# was itself checked against the JAGWAS binary (06), so a match carries that
# validation over to these counts.
#
# Finished replicates are skipped, so the command can be re-run safely.
#
# OVERLAP-ADJUSTED MATRIX. As in 05, --rhat_prefix reads 04's Rhat_adjfill /
# Rhat_adjdrop matrices and --out_scenario writes to the folder 05 used for them,
# so the check against 05 still applies.
#
# Usage:
#   Rscript 10_maf_class_pass.R --scenario=S1_normal --rep_start=1 --rep_end=257
#   Rscript 10_maf_class_pass.R --scenario=S4_normal_miss1pct --rhat_prefix=Rhat_adjfill \
#       --out_scenario=S4_adjfill --rep_start=1 --rep_end=20

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
  k            = "128",
  n_variant    = "8931083",
  layout       = "variant",
  maf_cut      = "0.01",
  chunk        = "500000",
  threads      = "1",
  overwrite    = "false",
  rhat_prefix  = "Rhat",       # reads <inbase>/<scenario>/<prefix>_rep###.txt, as in 05
  out_scenario = ""            # results folder under <inbase>; default = --scenario
)
opt       <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
out_scen  <- if (nzchar(opt$out_scenario)) opt$out_scenario else opt$scenario
k         <- as_int(opt$k, "k")
nv        <- as_num(opt$n_variant, "n_variant")
chunk     <- as_int(opt$chunk, "chunk")
maf_cut   <- as_num(opt$maf_cut, "maf_cut")
rep_start <- as_int(opt$rep_start, "rep_start")
rep_end   <- as_int(opt$rep_end, "rep_end")
overwrite <- as_bool(opt$overwrite, "overwrite")
thr       <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
layout_variant <- identical(opt$layout, "variant")
if (!opt$scenario %in% ALL_SCENARIOS) stop("Unknown scenario '", opt$scenario, "'")

# Same thresholds and histogram layout as 05_prepare_jagwas_input.R.
THR <- c(0.5, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8, 1e-9, 1e-10)
TNB <- 40001L; TBW <- 0.01   # T on [0, 400)
LNB <- 10001L; LBW <- 0.01   # -log10 p on [0, 100)
CLASSES <- c("common", "rare")

log_open(file.path(opt$inbase, "logs",
                   sprintf("10_maf_class_%s_rep%03d-%03d.log", out_scen, rep_start, rep_end)))
on.exit(log_close(), add = TRUE)

msg("MAF-class pass: ", opt$scenario, " replicates ", rep_start, "..", rep_end)
if (out_scen != opt$scenario || opt$rhat_prefix != "Rhat")
  msg("  matrix: ", opt$rhat_prefix, "_rep###.txt; results in ", file.path(opt$inbase, out_scen))
af <- fread(opt$af, select = c("variant_row", "minor_allele_frequency"), showProgress = FALSE)
if (nrow(af) != nv) stop(sprintf("%s has %d rows; the binary has %.0f variants.", opt$af, nrow(af), nv))
setorder(af, variant_row)
if (af$variant_row[1] != 1 || af$variant_row[nrow(af)] != nv)
  stop("variant_row is not 1..n_variant; cannot align the frequency file to the binary.")
if (anyNA(af$minor_allele_frequency))
  stop(sum(is.na(af$minor_allele_frequency)), " variants have no MAF in ", opt$af)
cls <- ifelse(af$minor_allele_frequency > maf_cut, 1L, 2L)
msg(sprintf("  common (MAF > %g): %d variants   rare (MAF <= %g): %d variants",
            maf_cut, sum(cls == 1L), maf_cut, sum(cls == 2L)))
rm(af)

# One chunk as a k x len matrix, one column per variant: exactly how the
# variant-major binary stores it.
read_chunk_t <- function(con, a, len) {
  if (layout_variant) {
    seek(con, where = (a - 1) * k * 4, origin = "start")
    matrix(readBin(con, "numeric", n = len * k, size = 4, endian = "little"),
           nrow = k, ncol = len)
  } else {
    Zt <- matrix(0, k, len)
    for (t in seq_len(k)) {
      seek(con, where = ((t - 1) * nv + (a - 1)) * 4, origin = "start")
      Zt[t, ] <- readBin(con, "numeric", n = len, size = 4, endian = "little")
    }
    Zt
  }
}

for (r in rep_start:rep_end) {
  tag  <- sprintf("rep%03d", r)
  fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
  frh  <- file.path(opt$inbase, opt$scenario, sprintf("%s_rep%03d.txt", opt$rhat_prefix, r))
  rdir <- file.path(opt$inbase, out_scen, tag)
  outs <- file.path(rdir, c(sprintf("direct_hist_%s.txt.gz", CLASSES),
                            sprintf("direct_counts_%s.txt", CLASSES)))
  if (!overwrite && all(file.exists(outs))) { msg(tag, ": done, skipping."); next }
  if (!file.exists(fbin)) { msg("MISSING ", fbin); next }
  if (!file.exists(frh))  { msg("MISSING ", frh);  next }
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()

  Rhat <- as.matrix(fread(frh, header = FALSE, showProgress = FALSE))
  dimnames(Rhat) <- NULL
  L <- t(chol(Rhat))                      # lower triangular, Rhat = L L'

  cnt <- matrix(0, 2, length(THR)); nc <- numeric(2)
  th  <- matrix(0, 2, TNB); lh <- matrix(0, 2, LNB)
  lthr <- -log10(THR)

  con <- file(fbin, "rb")
  for (a in seq(1, nv, by = chunk)) {
    b <- min(nv, a + chunk - 1); len <- as.integer(b - a + 1)
    Zt  <- read_chunk_t(con, a, len)
    Tst <- colSums(forwardsolve(L, Zt)^2)
    ok  <- is.finite(Tst)                 # a non-finite z in any trait gives NaN/Inf here
    ci  <- cls[a:b][ok]
    Tst <- Tst[ok]
    lp  <- -pchisq(Tst, df = k, lower.tail = FALSE, log.p = TRUE) / log(10)
    tb  <- as.integer(pmin(pmax(Tst, 0) / TBW, TNB - 1)) + 1L
    lb  <- as.integer(pmin(lp / LBW, LNB - 1)) + 1L
    for (cl in 1:2) {
      sel <- ci == cl
      nc[cl] <- nc[cl] + sum(sel)
      lpc <- lp[sel]
      cnt[cl, ] <- cnt[cl, ] + vapply(lthr, function(x) sum(lpc >= x), 0)
      th[cl, ]  <- th[cl, ] + tabulate(tb[sel], TNB)
      lh[cl, ]  <- lh[cl, ] + tabulate(lb[sel], LNB)
    }
  }
  close(con)

  # Common + rare must reproduce 05's counts for this replicate.
  check <- "no 05 counts to compare"
  f05 <- file.path(rdir, "direct_counts.txt")
  if (file.exists(f05)) {
    d05 <- fread(f05)
    m <- match(THR, d05$threshold)
    dmax <- max(abs(colSums(cnt) - d05$observed[m]), na.rm = TRUE)
    check <- if (dmax <= 1) "matches 05" else sprintf("** DIFFERS from 05 by up to %.0f **", dmax)
  }

  for (cl in 1:2) {
    nz_t <- which(th[cl, ] > 0); nz_l <- which(lh[cl, ] > 0)
    fwrite(rbind(data.table(type = "T",  bin = nz_t, count = th[cl, nz_t]),
                 data.table(type = "lp", bin = nz_l, count = lh[cl, nz_l])),
           file.path(rdir, sprintf("direct_hist_%s.txt.gz", CLASSES[cl])),
           sep = "\t", compress = "gzip")
  }
  for (cl in 1:2) {   # written last: their presence marks the replicate as done
    e <- nc[cl] * THR
    fwrite(data.table(threshold = THR, observed = cnt[cl, ], expected = e, ratio = cnt[cl, ] / e),
           file.path(rdir, sprintf("direct_counts_%s.txt", CLASSES[cl])), sep = "\t")
  }

  i58 <- which(THR == 5e-8)
  msg(sprintf("%s: P < 5e-8  common %.0f vs %.2f expected, rare %.0f vs %.2f expected | %s | %.1f min",
              tag, cnt[1, i58], nc[1] * 5e-8, cnt[2, i58], nc[2] * 5e-8, check,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
msg("Done.")
