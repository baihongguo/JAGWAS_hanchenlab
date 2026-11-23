#!/usr/bin/env Rscript
# 09_maf_and_univariate.R
#
# Asks where the S2 and S3 false positives come from, on a few replicates:
#
#   * are the JAGWAS hits concentrated in low-frequency variants?
#   * are the PER-TRAIT z-scores inflated in the same places?
#
# If both are yes, the cause is the normal approximation to the per-SNP score
# statistic breaking down for skewed or heavy-tailed phenotypes at variants with
# few minor-allele carriers -- the problem saddlepoint methods (SAIGE,
# fastGWA-SPA) exist to fix. That is upstream of JAGWAS and applies identically
# to MOSTest, minP and every other method built on the same summary statistics.
# If the JAGWAS hits were inflated while the univariate z-scores were clean, the
# fault would be in the combination step instead.
#
# lambda_GC cannot answer this: it is a median statistic, so it stays at 1.00
# while the tail misbehaves (S1 rep001: lambda_GC 0.9996 with every univariate
# distribution clean; S2: body calibrated, tail inflated 361-fold).
#
# Costs one pass over a replicate's binary, about the same as 05's first pass.
# Three replicates per scenario answer the question above; Supplementary Table 4
# pools all 257, run in chunks with --out_tag so that the chunks write separate
# files, and 14_summarise_maf.R pools them.
#
# Usage:
#   Rscript 09_maf_and_univariate.R --scenario=S2_lognormal --reps=1-3
#   Rscript 09_maf_and_univariate.R --scenario=S2_lognormal --reps=1-65 --out_tag=_reps1-65

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  tstat_root = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current",
  inbase     = "/data5/bguo1/JAGWAS_null_input",
  af         = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current/shared/variant_af.tsv",
  scenario   = "S2_lognormal",
  reps       = "1-3",
  k          = "128",
  n_variant  = "8931083",
  layout     = "variant",
  maf_breaks = "0,0.001,0.005,0.01,0.05,0.1,0.5",
  chunk      = "500000",
  threads    = "1",
  out_tag    = ""            # appended to the log and result file names
)
opt    <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k      <- as_int(opt$k, "k")
nv     <- as_num(opt$n_variant, "n_variant")
chunk  <- as_int(opt$chunk, "chunk")
thr    <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
layout_variant <- identical(opt$layout, "variant")
brk    <- as.numeric(trimws(strsplit(opt$maf_breaks, ",", fixed = TRUE)[[1]]))
reps   <- unlist(lapply(strsplit(opt$reps, ",", fixed = TRUE)[[1]], function(x) {
  ab <- as.integer(strsplit(trimws(x), "-", fixed = TRUE)[[1]])
  if (length(ab) == 2) ab[1]:ab[2] else ab
}))

# Thresholds for the 128-df JAGWAS statistic and for single traits. The
# univariate set includes 5e-8/128, the Bonferroni line the single-trait
# workflow would apply across 128 traits.
MTHR <- c(1e-4, 1e-5, 1e-6, 5e-8, 1e-10)
UTHR <- c(1e-4, 1e-6, 5e-8, 5e-8 / 128)

resdir <- file.path(opt$inbase, "results")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$inbase, "logs",
                   sprintf("09_maf_univariate_%s%s.log", opt$scenario, opt$out_tag)))
on.exit(log_close(), add = TRUE)

# ---------------------------------------------------------------- MAF bins --

msg("Reading allele frequencies: ", opt$af)
af <- fread(opt$af, showProgress = FALSE)
need <- c("variant_row", "minor_allele_frequency")
if (!all(need %in% names(af))) stop("Expected columns ", paste(need, collapse = ", "), " in ", opt$af)
if (nrow(af) != nv) stop(sprintf("%s has %d rows; the binary has %.0f variants.", opt$af, nrow(af), nv))
setorder(af, variant_row)
if (af$variant_row[1] != 1 || af$variant_row[nrow(af)] != nv)
  stop("variant_row is not 1..n_variant; cannot align the frequency file to the binary.")

maf   <- af$minor_allele_frequency
bin   <- cut(maf, breaks = brk, include.lowest = TRUE, right = FALSE)
blab  <- levels(bin)
bin   <- as.integer(bin)
nb    <- length(blab)
n_bin <- tabulate(bin, nb)
rm(af)

msg("MAF bins:")
for (i in seq_len(nb)) msg(sprintf("   %-16s %10d variants", blab[i], n_bin[i]))

read_chunk <- function(con, a, len) {
  if (layout_variant) {
    seek(con, where = (a - 1) * k * 4, origin = "start")
    matrix(readBin(con, "numeric", n = len * k, size = 4, endian = "little"),
           nrow = len, ncol = k, byrow = TRUE)
  } else {
    Z <- matrix(0, len, k)
    for (t in seq_len(k)) {
      seek(con, where = ((t - 1) * nv + (a - 1)) * 4, origin = "start")
      Z[, t] <- readBin(con, "numeric", n = len, size = 4, endian = "little")
    }
    Z
  }
}

mout <- list(); uout <- list(); tout <- list()

for (r in reps) {
  tag  <- sprintf("rep%03d", r)
  fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
  frh  <- file.path(opt$inbase, opt$scenario, sprintf("Rhat_rep%03d.txt", r))
  if (!file.exists(fbin)) { msg("MISSING ", fbin); next }
  if (!file.exists(frh))  { msg("MISSING ", frh);  next }
  t0 <- Sys.time()

  Rhat <- as.matrix(fread(frh, header = FALSE, showProgress = FALSE))
  dimnames(Rhat) <- NULL
  eR   <- eigen(Rhat, symmetric = TRUE)
  Rinv <- eR$vectors %*% (t(eR$vectors) / eR$values)

  mcrit <- qchisq(MTHR, df = k, lower.tail = FALSE)   # on the T scale
  ucrit <- qchisq(UTHR, df = 1, lower.tail = FALSE)   # on the z^2 scale

  m_cnt <- matrix(0, nb, length(MTHR))   # JAGWAS hits by MAF bin
  u_cnt <- matrix(0, nb, length(UTHR))   # univariate hits by MAF bin, summed over traits
  t_cnt <- matrix(0, k,  length(UTHR))   # univariate hits by trait
  n_seen <- numeric(nb)

  con <- file(fbin, "rb")
  for (a in seq(1, nv, by = chunk)) {
    b <- min(nv, a + chunk - 1); len <- as.integer(b - a + 1)
    Z <- read_chunk(con, a, len)
    ok <- is.finite(rowSums(Z))
    gi <- (a:b)[ok]                       # variant rows kept, aligned to the MAF file
    if (!all(ok)) Z <- Z[ok, , drop = FALSE]
    bi <- bin[gi]
    n_seen <- n_seen + tabulate(bi, nb)

    Tst <- rowSums((Z %*% Rinv) * Z)
    for (j in seq_along(mcrit))
      m_cnt[, j] <- m_cnt[, j] + tabulate(bi[Tst >= mcrit[j]], nb)

    for (t in seq_len(k)) {
      z2 <- Z[, t]^2
      for (j in seq_along(ucrit)) {
        hit <- z2 >= ucrit[j]
        nh <- sum(hit)
        if (nh) {
          u_cnt[, j] <- u_cnt[, j] + tabulate(bi[hit], nb)
          t_cnt[t, j] <- t_cnt[t, j] + nh
        }
      }
    }
  }
  close(con)

  mout[[length(mout) + 1L]] <- data.table(
    scenario = opt$scenario, replicate = r,
    maf_bin = rep(blab, length(MTHR)), n_variants = rep(n_seen, length(MTHR)),
    threshold = rep(MTHR, each = nb), observed = as.vector(m_cnt),
    expected = as.vector(outer(n_seen, MTHR)))
  uout[[length(uout) + 1L]] <- data.table(
    scenario = opt$scenario, replicate = r,
    maf_bin = rep(blab, length(UTHR)), n_variants = rep(n_seen, length(UTHR)),
    threshold = rep(UTHR, each = nb), observed = as.vector(u_cnt),
    expected = as.vector(outer(n_seen * k, UTHR)))
  tout[[length(tout) + 1L]] <- data.table(
    scenario = opt$scenario, replicate = r,
    trait = rep(seq_len(k), length(UTHR)),
    threshold = rep(UTHR, each = k), observed = as.vector(t_cnt),
    expected = rep(sum(n_seen) * UTHR, each = k))

  msg(sprintf("%s %s (%.1f min)", opt$scenario, tag,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  msg("   JAGWAS P < 5e-8 by MAF bin (observed / expected):")
  j <- which(MTHR == 5e-8)
  for (i in seq_len(nb))
    msg(sprintf("      %-16s %10.0f / %8.2f = %.4g", blab[i], m_cnt[i, j],
                n_seen[i] * MTHR[j], m_cnt[i, j] / (n_seen[i] * MTHR[j])))
  msg("   single-trait P < 5e-8 by MAF bin, summed over 128 traits:")
  j <- which(UTHR == 5e-8)
  for (i in seq_len(nb))
    msg(sprintf("      %-16s %10.0f / %8.2f = %.4g", blab[i], u_cnt[i, j],
                n_seen[i] * k * UTHR[j], u_cnt[i, j] / (n_seen[i] * k * UTHR[j])))
}

if (length(mout)) {
  m <- rbindlist(mout); m[, ratio := observed / expected]
  u <- rbindlist(uout); u[, ratio := observed / expected]
  tt <- rbindlist(tout); tt[, ratio := observed / expected]
  f1 <- file.path(resdir, sprintf("maf_multivariate_%s%s.txt", opt$scenario, opt$out_tag))
  f2 <- file.path(resdir, sprintf("maf_univariate_%s%s.txt", opt$scenario, opt$out_tag))
  f3 <- file.path(resdir, sprintf("univariate_by_trait_%s%s.txt", opt$scenario, opt$out_tag))
  fwrite(m, f1, sep = "\t"); fwrite(u, f2, sep = "\t"); fwrite(tt, f3, sep = "\t")
  msg("")
  msg("Wrote ", f1)
  msg("Wrote ", f2)
  msg("Wrote ", f3)
  msg("")
  msg("Reading it: if both tables are inflated in the low-MAF bins and near 1 for")
  msg("MAF >= 0.01, the cause is the per-SNP normal approximation at variants with")
  msg("few carriers -- upstream of JAGWAS, and shared by every method using these")
  msg("summary statistics.")
}
msg("Done.")
