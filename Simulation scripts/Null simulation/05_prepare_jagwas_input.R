#!/usr/bin/env Rscript
# 05_prepare_jagwas_input.R
#
# Converts the torchGWAS binary t-statistics into JAGWAS input, one replicate
# at a time:
#
#   <outdir>/<scenario>/rep###/QT0.fastGWA.gz ... QT127.fastGWA.gz
#   <outdir>/<scenario>/rep###/qc_summary.txt, qc_traits.txt
#   <outdir>/<scenario>/rep###/direct_counts.txt, direct_hits.txt
#
# FILE LAYOUT. JAGWAS reads columns by POSITION; the header line is skipped
# unread. With --score_test 0 --beta_se 1 it parses AF = column 7, BETA =
# column 8, SE = column 9, and echoes columns 1-7 of the FIRST file only into
# its output. So every file is
#
#   CHR  SNP  POS  A1  A2  N  AF1  BETA  SE
#
#   BETA = the torchGWAS t-statistic and SE = 1, so JAGWAS's BETA/SE is t.
#   AF1  = one constant, identical in every file. JAGWAS drops any variant whose
#          AF differs between files, so it must match exactly, and it must lie
#          strictly inside (MAF, 1 - MAF).
#   CHR SNP POS A1 A2 N are filled only in the first file -- the only one
#          JAGWAS echoes -- and are NA elsewhere. Without --variants, SNP is the
#          variant's 1-based index in tstat.f32, so every JAGWAS output row can
#          be traced back to the binary.
#
# WHY NOT --beta_se 0 WITH A Z COLUMN. In its z-score branch JAGWAS computes
# z' R z rather than z' R^-1 z (it multiplies by `cor` where the other two
# branches use `cor_i = inv(cor)`). --beta_se 1 with SE = 1 is just as simple
# and goes through the correct branch -- the one the real UDIP analysis used.
# These files also work unchanged with --beta_se 0 once that bug is fixed.
#
# GZIP. JAGWAS opens its inputs with std::ifstream and cannot decompress;
# run_jagwas_null.sh streams each file through zcat into a named pipe instead.
#
# THE CORRELATION MATRIX. <outdir>/<scenario>/Rhat_rep###.txt must already
# exist; 04_estimate_scenario_R.R makes it. This script only reads it, for the
# QC and the direct statistic.
#
# ORDER. Phenotype column QT001 = torchGWAS trait 1 = QT0.fastGWA.gz = row and
# column 1 of Rhat_rep###.txt. Pass the files to JAGWAS as QT0 ... QT127.
#
# BINARY LAYOUT. tstat.f32 is VARIANT-major (confirmed with Ziqian, 2026-09-20):
# the 128 traits of variant 1 come first, then the 128 of variant 2, so
# value (variant v, trait t) sits at index (v-1)*n_trait + t. Reading it as
# trait-major instead mixes unrelated variants and traits into each "variant";
# their z-scores are then uncorrelated, and because tr(Rhat^-1) = 8210 against
# k = 128, every variant clears every significance threshold. Use
# --layout=trait only if a future file really is trait-major; 08_check_layout.R
# tells the two apart.
#
# TWO PASSES over tstat.f32:
#   1. variant chunks, read contiguously: drop variants non-finite in any trait
#      (from ALL files, keeping them aligned), QC, and the direct statistic. If
#      the per-trait files are wanted, this pass also lays down a transposed
#      scratch copy (4.6 GB, deleted afterwards), so the original is read once.
#   2. trait by trait, from that scratch copy: each file in one gzip stream.
#
# QC AND DIRECT STATISTIC (pass 1):
#   * E[T]/k compares the z-scores' empirical cross-trait covariance with Rhat.
#     S1-S3 should give ~1.00; S4 ~1.63, which IS the scenario. ~1.00 for S1
#     also confirms that the binary is trait-major and the trait order is right.
#   * T = z' Rhat^-1 z is computed directly and type I error counts recorded at
#     a grid of thresholds. Check against JAGWAS on a few replicates; if they
#     agree, --write_sumstats=false skips pass 2 entirely.
#
# Usage:
#   Rscript 05_prepare_jagwas_input.R --scenario=S1_normal --rep_start=1 --rep_end=1 \
#     --variants=<torchGWAS variant list>

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  tstat_root     = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current",
  outdir         = "/data5/bguo1/JAGWAS_null_input",
  scenario       = "S1_normal",
  rep_start      = "1",
  rep_end        = "1",
  k              = "128",
  n_variant      = "8931083",
  layout         = "variant",     # "variant" (v-1)*k + t, or "trait" (t-1)*nv + v
  n              = "22250",       # analysed sample size, written to the N column
  af             = "0.3",         # placeholder allele frequency
  variants       = "",            # .bim / .pvar / header'd table, same order as tstat.f32
  file_prefix    = "QT",
  first_index    = "0",           # QT0 ... QT127, matching the real JAGWAS runs
  gzip           = "true",
  write_sumstats = "true",
  save_direct_p  = "false",       # write -log10 p as float32 for comparison with JAGWAS
  hit_p          = "5e-8",
  max_hits       = "200000",
  chunk          = "500000",      # variants per pass-1 chunk (x 128 traits x 8 bytes)
  threads        = "0",
  overwrite      = "false",
  rhat_prefix    = "Rhat",        # reads <outdir>/<scenario>/<prefix>_rep###.txt;
                                  # Rhat_adjfill / Rhat_adjdrop for 04's overlap-adjusted matrices
  out_scenario   = ""             # results folder under <outdir>; default = --scenario.
                                  # Use a new name (e.g. S4_adjfill) to keep the naive results
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
out_scen <- if (nzchar(opt$out_scenario)) opt$out_scenario else opt$scenario

k          <- as_int(opt$k, "k")
nv         <- as_num(opt$n_variant, "n_variant")     # double: byte offsets exceed 2^31
n_nominal  <- as_num(opt$n, "n")
af         <- as_num(opt$af, "af")
rep_start  <- as_int(opt$rep_start, "rep_start")
rep_end    <- as_int(opt$rep_end, "rep_end")
first_idx  <- as_int(opt$first_index, "first_index")
chunk      <- as_int(opt$chunk, "chunk")
hit_p      <- as_num(opt$hit_p, "hit_p")
max_hits   <- as_int(opt$max_hits, "max_hits")
gz         <- as_bool(opt$gzip, "gzip")
write_ss   <- as_bool(opt$write_sumstats, "write_sumstats")
save_p     <- as_bool(opt$save_direct_p, "save_direct_p")
overwrite  <- as_bool(opt$overwrite, "overwrite")
threads    <- as_int(opt$threads, "threads")
if (threads > 0) setDTthreads(threads)
if (!(af > 0.001 && af < 0.999)) stop("--af must lie in (0.001, 0.999).")
if (!opt$layout %in% c("variant", "trait")) stop("--layout must be 'variant' or 'trait'.")
layout_variant <- identical(opt$layout, "variant")
if (!opt$scenario %in% ALL_SCENARIOS)
  stop("Unknown scenario '", opt$scenario, "'. Valid: ", paste(ALL_SCENARIOS, collapse = ", "))

THR <- c(0.5, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8, 1e-9, 1e-10)
NB  <- 20001L   # z^2 histogram for lambda_GC: bins of 1e-4 on [0, 2), plus overflow
BW  <- 1e-4
# Histograms of the JAGWAS statistic, saved per replicate so 07 can pool them:
TNB <- 40001L; TBW <- 0.01   # T on [0, 400) in steps of 0.01 -> lambda_GC
LNB <- 10001L; LBW <- 0.01   # -log10 p on [0, 100) in steps of 0.01 -> QQ plot

file_names <- sprintf("%s%d.fastGWA%s", opt$file_prefix,
                      first_idx + seq_len(k) - 1L, if (gz) ".gz" else "")

log_open(file.path(opt$outdir, "logs",
                   sprintf("05_prepare_%s_rep%03d-%03d.log", out_scen, rep_start, rep_end)))
on.exit(log_close(), add = TRUE)

msg("JAGWAS input preparation: ", opt$scenario, " replicates ", rep_start, "..", rep_end)
if (out_scen != opt$scenario || opt$rhat_prefix != "Rhat")
  msg("  matrix: ", opt$rhat_prefix, "_rep###.txt; results in ", file.path(opt$outdir, out_scen))

# ------------------------------------------------------------------ helpers --

read_variants <- function(path) {
  if (!nzchar(path)) {
    msg("  No --variants file: SNP = index in tstat.f32; CHR, POS, A1, A2 = NA.")
    msg("  Real positions are only needed to map hits back to the genome.")
    return(NULL)
  }
  msg("  Reading variants: ", path)
  if (grepl("\\.bim$", path)) {
    # Alleles forced to character: a column of T alleles must not become logical.
    v <- fread(path, header = FALSE, select = c(1, 2, 4, 5, 6),
               colClasses = list(character = c(2, 5, 6)),
               col.names = c("CHR", "SNP", "POS", "A1", "A2"), showProgress = FALSE)
  } else {
    v <- fread(path, header = TRUE, showProgress = FALSE)
    nm <- toupper(sub("^#", "", names(v)))
    pick <- function(c) { i <- match(c, nm); i[!is.na(i)][1] }
    ic <- c(CHR = pick(c("CHR", "CHROM")), SNP = pick(c("SNP", "ID", "RSID")),
            POS = pick(c("POS", "BP")), A1 = pick(c("A1", "ALT")), A2 = pick(c("A2", "REF")))
    if (anyNA(ic)) stop("Cannot find CHR/SNP/POS/A1/A2 columns in ", path, call. = FALSE)
    v <- v[, ic, with = FALSE]
    setnames(v, names(ic))
  }
  if (nrow(v) != nv)
    stop(sprintf("%s has %d variants; the binary has %.0f.", path, nrow(v), nv), call. = FALSE)
  as.list(v)
}

# One chunk of `len` variants x k traits, from either layout.
read_chunk <- function(con, a, len) {
  if (layout_variant) {
    seek(con, where = (a - 1) * k * 4, origin = "start")
    matrix(readBin(con, what = "numeric", n = len * k, size = 4, endian = "little"),
           nrow = len, ncol = k, byrow = TRUE)
  } else {
    Z <- matrix(0, len, k)
    for (t in seq_len(k)) {
      seek(con, where = ((t - 1) * nv + (a - 1)) * 4, origin = "start")
      Z[, t] <- readBin(con, what = "numeric", n = len, size = 4, endian = "little")
    }
    Z
  }
}

sat_tail <- function(lam, q) {
  c1 <- sum(lam); c2 <- 2 * sum(lam^2)
  pchisq(q / (c2 / (2 * c1)), df = 2 * c1^2 / c2, lower.tail = FALSE)
}

vars <- read_variants(opt$variants)

# ------------------------------------------------------------ one replicate --

process_rep <- function(r) {
  tag  <- sprintf("rep%03d", r)
  fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
  rdir <- file.path(opt$outdir, out_scen, tag)
  done <- file.path(rdir, "qc_summary.txt")
  if (!overwrite && file.exists(done) &&
      (!write_ss || all(file.exists(file.path(rdir, file_names))))) {
    msg(tag, ": already prepared, skipping.")
    return(invisible(NULL))
  }
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(done)) file.remove(done)
  t0 <- Sys.time()

  if (!file.exists(fbin)) stop("Missing binary: ", fbin)
  if (file.size(fbin) != nv * k * 4)
    stop(sprintf("%s is %.0f bytes; expected %.0f (%.0f variants x %d traits x 4).",
                 fbin, file.size(fbin), nv * k * 4, nv, k))

  # --- Rhat, made beforehand by 04_estimate_scenario_R.R ---------------------
  frhat <- file.path(opt$outdir, opt$scenario, sprintf("%s_rep%03d.txt", opt$rhat_prefix, r))
  if (!file.exists(frhat))
    stop("Missing ", frhat, ": run 04_estimate_scenario_R.R first.", call. = FALSE)
  Rhat <- as.matrix(fread(frhat, header = FALSE, showProgress = FALSE))
  dimnames(Rhat) <- NULL
  if (!identical(dim(Rhat), c(k, k)))
    stop(sprintf("%s is %d x %d, expected %d x %d.", frhat, nrow(Rhat), ncol(Rhat), k, k))
  eR <- eigen(Rhat, symmetric = TRUE)
  if (min(eR$values) <= 0) stop("Rhat is not positive definite for ", tag)
  Rinv   <- eR$vectors %*% (t(eR$vectors) / eR$values)
  Rmhalf <- eR$vectors %*% (t(eR$vectors) / sqrt(eR$values))
  N_t <- rep(n_nominal, k)

  # --- pass 1: variant chunks x all traits -----------------------------------
  con <- file(fbin, "rb")
  on.exit(close(con), add = TRUE)
  pcon <- NULL
  if (save_p) {
    pcon <- file(file.path(rdir, "direct_neglog10p.f32"), "wb")
    on.exit(close(pcon), add = TRUE)
  }

  # Scratch trait-major copy, written as pass 1 goes, so pass 2 can read each
  # trait contiguously without re-reading the variant-major original.
  tmpf <- if (write_ss && layout_variant) file.path(rdir, "tstat_traitmajor.tmp") else NULL
  tmpcon <- NULL
  if (!is.null(tmpf)) {
    on.exit(unlink(tmpf), add = TRUE)
    tmpcon <- file(tmpf, "w+b")   # w+b: seeking is supported
  }

  keep_var <- logical(nv)
  S <- matrix(0, k, k); s1 <- numeric(k); n_ok <- 0
  hist <- matrix(0, k, NB)
  cnt <- numeric(length(THR)); lthr <- -log10(THR)
  hits <- list(); n_hits <- 0; n_hits_all <- 0
  thist <- numeric(TNB); lhist <- numeric(LNB)

  for (a in seq(1, nv, by = chunk)) {
    b <- min(nv, a + chunk - 1); len <- as.integer(b - a + 1)

    Z <- read_chunk(con, a, len)

    # A variant non-finite in ANY trait is dropped from ALL files.
    ok <- is.finite(rowSums(Z))
    keep_var[a:b] <- ok

    # Transposed before the dropped rows are removed, so trait t keeps a fixed
    # offset in the scratch file and dropping stays a pass-2 concern.
    if (!is.null(tmpcon)) for (t in seq_len(k)) {
      seek(tmpcon, where = ((t - 1) * nv + (a - 1)) * 4, origin = "start")
      writeBin(Z[, t], tmpcon, size = 4)
    }

    if (!all(ok)) Z <- Z[ok, , drop = FALSE]
    vidx <- (a:b)[ok]
    n_ok <- n_ok + nrow(Z)

    S  <- S + crossprod(Z)
    s1 <- s1 + colSums(Z)
    for (t in seq_len(k))
      hist[t, ] <- hist[t, ] + tabulate(as.integer(pmin(Z[, t]^2 / BW, NB - 1)) + 1L, NB)

    Tst <- rowSums((Z %*% Rinv) * Z)
    lp  <- -pchisq(Tst, df = k, lower.tail = FALSE, log.p = TRUE) / log(10)
    cnt <- cnt + vapply(lthr, function(x) sum(lp >= x), 0)
    thist <- thist + tabulate(as.integer(pmin(pmax(Tst, 0) / TBW, TNB - 1)) + 1L, TNB)
    lhist <- lhist + tabulate(as.integer(pmin(lp / LBW, LNB - 1)) + 1L, LNB)
    h <- which(lp >= -log10(hit_p))
    n_hits_all <- n_hits_all + length(h)
    if (length(h) && n_hits < max_hits) {
      h <- head(h, max_hits - n_hits)
      g <- vidx[h]
      d <- data.table(index = g)
      if (!is.null(vars)) d[, `:=`(SNP = vars$SNP[g], CHR = vars$CHR[g], POS = vars$POS[g])]
      d[, `:=`(T = Tst[h], neglog10p = lp[h])]
      hits[[length(hits) + 1L]] <- d
      n_hits <- n_hits + length(h)
    }
    if (!is.null(pcon)) writeBin(lp, pcon, size = 4)
  }
  n_drop <- nv - n_ok
  if (!is.null(tmpcon)) { close(tmpcon); tmpcon <- NULL }

  # --- pass 2: one file per trait, each a single gzip stream -----------------
  if (write_ss) {
    src <- file(if (!is.null(tmpf)) tmpf else fbin, "rb")   # trait-major either way
    idx <- which(keep_var)
    n <- length(idx)
    na_col <- rep(NA, n)
    af_col <- rep.int(af, n)
    se_col <- rep.int(1L, n)
    id_first <- if (is.null(vars)) {
      list(CHR = na_col, SNP = idx, POS = na_col, A1 = na_col, A2 = na_col)
    } else {
      list(CHR = vars$CHR[idx], SNP = vars$SNP[idx], POS = vars$POS[idx],
           A1 = vars$A1[idx], A2 = vars$A2[idx])
    }
    id_rest <- list(CHR = na_col, SNP = na_col, POS = na_col, A1 = na_col, A2 = na_col,
                    N = na_col)
    for (t in seq_len(k)) {
      seek(src, where = (t - 1) * nv * 4, origin = "start")
      z <- readBin(src, what = "numeric", n = as.integer(nv), size = 4, endian = "little")[idx]
      cols <- c(if (t == 1) c(id_first, list(N = rep.int(as.integer(N_t[1]), n))) else id_rest,
                list(AF1 = af_col, BETA = signif(z, 8), SE = se_col))
      fwrite(cols, file.path(rdir, file_names[t]), sep = "\t", quote = FALSE, na = "NA",
             compress = if (gz) "gzip" else "none")
    }
    close(src)
    if (!is.null(tmpf)) unlink(tmpf)
  }

  # --- QC ---------------------------------------------------------------------
  C   <- S / n_ok
  M   <- Rmhalf %*% C %*% Rmhalf
  lam <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  pred <- sat_tail(lam, qchisq(5e-8, df = k, lower.tail = FALSE))
  dcor <- max(abs(cov2cor(C) - Rhat)[lower.tri(Rhat)])

  med <- apply(hist, 1, function(h) {
    bin <- which(cumsum(h) >= n_ok / 2)[1]
    if (bin == NB) NA_real_ else (bin - 0.5) * BW
  })
  lgc <- med / qchisq(0.5, 1)

  tbin  <- which(cumsum(thist) >= n_ok / 2)[1]
  T_med <- if (tbin == TNB) NA_real_ else (tbin - 0.5) * TBW
  lgc_T <- T_med / qchisq(0.5, df = k)

  # Printed before anything is written, so a later failure cannot lose them.
  i58 <- which(THR == 5e-8)
  msg(sprintf("%s: %.0f variants written (%.0f dropped as non-finite)", tag, n_ok, n_drop))
  msg(sprintf("   z-cov vs Rhat: lambda [%.4f, %.4f], E[T]/k = %.4f, max|corZ - Rhat| = %.4f",
              min(lam), max(lam), sum(lam) / k, dcor))
  msg(sprintf("   univariate lambda_GC median %.4f [%.4f, %.4f]",
              median(lgc, na.rm = TRUE), min(lgc, na.rm = TRUE), max(lgc, na.rm = TRUE)))
  msg(sprintf("   direct z'Rhat^-1 z: lambda_GC %.4f, P < 5e-8 %.0f observed vs %.2f expected (ratio %.4g)",
              lgc_T, cnt[i58], n_ok * THR[i58], cnt[i58] / (n_ok * THR[i58])))
  naive_s4 <- opt$scenario == "S4_normal_miss1pct" && opt$rhat_prefix == "Rhat"
  if (!naive_s4 && abs(sum(lam) / k - 1) > 0.02)
    msg("   ** E[T]/k is far from 1: the matrix does not match the z-scores. Check the keep",
        " list and the matrix used (for S4: adjfill versus adjdrop), the binary layout and the",
        " trait order.")
  if (naive_s4)
    msg("   (S4: E[T]/k ~ 1.63 means the overlap mismatch is present, as intended. ~1.00 would",
        "\n    mean torchGWAS used complete cases only, and S4 would not test what we designed.)")

  nz_t <- which(thist > 0); nz_l <- which(lhist > 0)
  fwrite(rbind(data.table(type = "T",  bin = nz_t, count = thist[nz_t]),
               data.table(type = "lp", bin = nz_l, count = lhist[nz_l])),
         file.path(rdir, "direct_hist.txt.gz"), sep = "\t", compress = "gzip")

  fwrite(data.table(trait = seq_len(k), file = file_names, N = as.integer(N_t),
                    mean_z = s1 / n_ok, mean_z2 = diag(C), lambda_gc = lgc),
         file.path(rdir, "qc_traits.txt"), sep = "\t")

  expd <- n_ok * THR
  fwrite(data.table(threshold = THR, observed = cnt, expected = expd, ratio = cnt / expd),
         file.path(rdir, "direct_counts.txt"), sep = "\t")
  if (length(hits)) fwrite(rbindlist(hits), file.path(rdir, "direct_hits.txt"), sep = "\t")

  summ <- data.table(
    field = c("scenario", "rhat_prefix", "replicate", "n_variant_binary", "n_variant_written",
            "n_variant_dropped_nonfinite", "n_analysed_min", "n_analysed_max",
            "Rhat_lambda_min", "Rhat_condition", "tr_Rhat_inv",
            "zcov_vs_Rhat_lambda_min", "zcov_vs_Rhat_lambda_max", "E_T_over_k",
            "predicted_type1_ratio_5e-8", "max_abs_corZ_minus_Rhat",
            "lambda_gc_median", "lambda_gc_min", "lambda_gc_max",
            "jagwas_median_T", "jagwas_lambda_gc",
            "direct_obs_5e-8", "direct_exp_5e-8", "direct_ratio_5e-8",
            "direct_hits_recorded", "direct_hits_total", "minutes"),
    value = c(opt$scenario, opt$rhat_prefix, r, format(nv, scientific = FALSE),
              format(n_ok, scientific = FALSE), format(n_drop, scientific = FALSE),
              min(N_t), max(N_t),
              min(eR$values), max(eR$values) / min(eR$values), sum(1 / eR$values),
              min(lam), max(lam), sum(lam) / k, pred / 5e-8, dcor,
              median(lgc, na.rm = TRUE), min(lgc, na.rm = TRUE), max(lgc, na.rm = TRUE),
              T_med, lgc_T,
              cnt[i58], expd[i58], cnt[i58] / expd[i58], n_hits, n_hits_all,
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)))
  # data.table() reads a `key` argument as the table's key rather than a column,
  # so the column is built as `field` and renamed here.
  setnames(summ, "field", "key")
  fwrite(summ, done, sep = "\t")   # written last: its presence marks completion

  msg(sprintf("   %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

for (r in rep_start:rep_end) process_rep(r)
msg("Done.")
