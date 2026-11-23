#!/usr/bin/env Rscript
# 12_trait_subsets.R
#
# The null simulation for fewer phenotypes, from the existing 128-trait summary
# statistics: per replicate, a random subset of k of the 128 traits is tested
# with T = z_S' Rhat_SS^-1 z_S ~ chi2_k, where z_S are those traits' z-scores
# and Rhat_SS the matching block of the replicate's Rhat from 04. No new GWAS
# is needed: the subset block of Rhat is exactly the matrix 04 would estimate
# from those traits alone (for Rhat_adjfill too, entry by entry).
#
# One pass over tstat.f32 per replicate serves every k in --subset_k and both
# MAF classes, and writes, for each k, what 05 and 10 write for all 128 traits:
#
#   <inbase>/<out_scenario>_k<k>/rep###/
#     qc_summary.txt, qc_traits.txt, direct_counts.txt, direct_hist.txt.gz,
#     direct_hits.txt                                            (as 05)
#     direct_counts_common.txt, direct_hist_common.txt.gz,
#     direct_counts_rare.txt,   direct_hist_rare.txt.gz         (as 10)
#
# so 07 and 11 run on these folders unchanged; they read k from the _k suffix.
#
# WHICH TRAITS. trait_subset() in common.R: one random permutation of the 128
# traits per replicate (seed --subset_seed + replicate); the first k are kept.
# The 16 are therefore inside the 64, and all scenarios use the same subset for
# a given replicate. qc_traits.txt lists the traits used.
#
# VARIANTS. A variant non-finite in any of the 128 traits is dropped, as in 05,
# so every k is tested on exactly the variants of the 128-trait analysis.
#
# CHECK. --subset_k=128 runs the full set through this script (folder _k128)
# and compares the result with 05's direct_counts.txt and 10's class counts for
# the same replicate; the log says "matches 05/10". Do this on one replicate
# before the full run.
#
# Usage:
#   Rscript 12_trait_subsets.R --scenario=S1_normal --subset_k=16,64,128 --rep_start=1 --rep_end=1
#   Rscript 12_trait_subsets.R --scenario=S1_normal --rep_start=2 --rep_end=257
#   Rscript 12_trait_subsets.R --scenario=S4_normal_miss1pct --rhat_prefix=Rhat_adjfill \
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
  subset_k     = "16,64",       # numbers of traits to test
  subset_seed  = "20261008",    # replicate r uses seed subset_seed + r
  n_variant    = "8931083",
  layout       = "variant",
  maf_cut      = "0.01",
  hit_p        = "5e-8",
  max_hits     = "200000",
  chunk        = "500000",
  threads      = "1",
  overwrite    = "false",
  rhat_prefix  = "Rhat",        # reads <inbase>/<scenario>/<prefix>_rep###.txt, as in 05
  out_scenario = ""             # folder base name; default = --scenario. _k<k> is appended
)
opt       <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
out_base  <- if (nzchar(opt$out_scenario)) opt$out_scenario else opt$scenario
k_all     <- as_int(opt$k, "k")
ks        <- sort(unique(as.integer(trimws(strsplit(opt$subset_k, ",", fixed = TRUE)[[1]]))))
seed_base <- as_int(opt$subset_seed, "subset_seed")
nv        <- as_num(opt$n_variant, "n_variant")
chunk     <- as_int(opt$chunk, "chunk")
maf_cut   <- as_num(opt$maf_cut, "maf_cut")
hit_p     <- as_num(opt$hit_p, "hit_p")
max_hits  <- as_int(opt$max_hits, "max_hits")
rep_start <- as_int(opt$rep_start, "rep_start")
rep_end   <- as_int(opt$rep_end, "rep_end")
overwrite <- as_bool(opt$overwrite, "overwrite")
thr       <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
layout_variant <- identical(opt$layout, "variant")
if (anyNA(ks) || any(ks < 1 | ks > k_all)) stop("--subset_k must be integers in 1..", k_all)
if (!opt$scenario %in% ALL_SCENARIOS) stop("Unknown scenario '", opt$scenario, "'")
if (grepl("_k[0-9]+$", out_base)) stop("--out_scenario must not end in _k<number>; that suffix is added.")

# Same thresholds and histogram layout as 05 and 10.
THR <- c(0.5, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 5e-8, 1e-8, 1e-9, 1e-10)
NB  <- 20001L; BW  <- 1e-4    # z^2 on [0, 2), univariate lambda_GC
TNB <- 40001L; TBW <- 0.01    # T on [0, 400)
LNB <- 10001L; LBW <- 0.01    # -log10 p on [0, 100)
CLASSES <- c("common", "rare")
lthr <- -log10(THR)
i58  <- which(THR == 5e-8)

log_open(file.path(opt$inbase, "logs",
                   sprintf("12_subsets_%s_k%s_rep%03d-%03d.log", out_base,
                           paste(ks, collapse = "-"), rep_start, rep_end)))
on.exit(log_close(), add = TRUE)

msg("Trait subsets k = ", paste(ks, collapse = ", "), " of ", k_all, ": ", opt$scenario,
    " replicates ", rep_start, "..", rep_end)
msg("  matrix: ", opt$rhat_prefix, "_rep###.txt; results in ",
    file.path(opt$inbase, paste0(out_base, "_k<k>")))

af <- fread(opt$af, select = c("variant_row", "minor_allele_frequency"), showProgress = FALSE)
if (nrow(af) != nv) stop(sprintf("%s has %d rows; the binary has %.0f variants.", opt$af, nrow(af), nv))
setorder(af, variant_row)
if (af$variant_row[1] != 1 || af$variant_row[nrow(af)] != nv)
  stop("variant_row is not 1..n_variant; cannot align the frequency file to the binary.")
if (anyNA(af$minor_allele_frequency))
  stop(sum(is.na(af$minor_allele_frequency)), " variants have no MAF in ", opt$af)
cls <- ifelse(af$minor_allele_frequency > maf_cut, 1L, 2L)
rm(af)

# One chunk as a k_all x len matrix, one column per variant (the binary's own order).
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

sat_tail <- function(lam, q) {
  c1 <- sum(lam); c2 <- 2 * sum(lam^2)
  pchisq(q / (c2 / (2 * c1)), df = 2 * c1^2 / c2, lower.tail = FALSE)
}
hist_median <- function(h, n, w) {
  b <- which(cumsum(h) >= n / 2)[1]
  if (is.na(b) || b == length(h)) NA_real_ else (b - 0.5) * w
}
write_hist <- function(th, lh, f) {
  nz_t <- which(th > 0); nz_l <- which(lh > 0)
  fwrite(rbind(data.table(type = "T",  bin = nz_t, count = th[nz_t]),
               data.table(type = "lp", bin = nz_l, count = lh[nz_l])),
         f, sep = "\t", compress = "gzip")
}
write_counts <- function(cnt, n, f) {
  e <- n * THR
  fwrite(data.table(threshold = THR, observed = cnt, expected = e, ratio = cnt / e), f, sep = "\t")
}
compare_counts <- function(cnt, f) {
  if (!file.exists(f)) return(NA_real_)
  d <- fread(f)
  max(abs(cnt - d$observed[match(THR, d$threshold)]), na.rm = TRUE)
}

for (r in rep_start:rep_end) {
  tag  <- sprintf("rep%03d", r)
  fbin <- file.path(opt$tstat_root, opt$scenario, tag, "sumstats", "tstat.f32")
  frh  <- file.path(opt$inbase, opt$scenario, sprintf("%s_rep%03d.txt", opt$rhat_prefix, r))
  rdirs <- setNames(file.path(opt$inbase, sprintf("%s_k%d", out_base, ks), tag), ks)
  done  <- file.path(rdirs, "qc_summary.txt")
  if (!overwrite && all(file.exists(done))) { msg(tag, ": done, skipping."); next }
  if (!file.exists(fbin)) { msg("MISSING ", fbin); next }
  if (!file.exists(frh))  { msg("MISSING ", frh);  next }
  if (file.size(fbin) != nv * k_all * 4)
    stop(sprintf("%s is %.0f bytes; expected %.0f.", fbin, file.size(fbin), nv * k_all * 4))
  t0 <- Sys.time()

  Rhat <- as.matrix(fread(frh, header = FALSE, showProgress = FALSE))
  dimnames(Rhat) <- NULL
  if (!identical(dim(Rhat), c(k_all, k_all))) stop(frh, " is not ", k_all, " x ", k_all)

  # Per k: its traits, the Cholesky factor of its block of Rhat, and accumulators.
  S <- lapply(ks, function(kk) {
    sel <- trait_subset(k_all, kk, seed_base, r)
    Rs  <- Rhat[sel, sel, drop = FALSE]
    ev  <- eigen(Rs, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) stop(sprintf("%s: Rhat block for k = %d is not positive definite.", tag, kk))
    list(k = kk, sel = sel, Rs = Rs, ev = ev, L = t(chol(Rs)),
         cnt = numeric(length(THR)), th = numeric(TNB), lh = numeric(LNB),
         ccnt = matrix(0, 2, length(THR)), cth = matrix(0, 2, TNB), clh = matrix(0, 2, LNB),
         nc = numeric(2), hits = list(), n_hits = 0, n_hits_all = 0)
  })
  names(S) <- ks
  u <- sort(unique(unlist(lapply(S, `[[`, "sel"))))     # traits used by any k
  pos <- lapply(S, function(s) match(s$sel, u))          # their rows within u
  Szz <- matrix(0, length(u), length(u)); s1 <- numeric(length(u))
  uh  <- matrix(0, length(u), NB)
  n_ok <- 0

  con <- file(fbin, "rb")
  for (a in seq(1, nv, by = chunk)) {
    b <- min(nv, a + chunk - 1); len <- as.integer(b - a + 1)
    Zt <- read_chunk_t(con, a, len)
    ok <- is.finite(colSums(Zt))                        # all 128 traits, as in 05
    vidx <- (a:b)[ok]
    Zu <- Zt[u, ok, drop = FALSE]
    rm(Zt)
    ci <- cls[vidx]
    n_ok <- n_ok + length(vidx)

    Szz <- Szz + tcrossprod(Zu)
    s1  <- s1 + rowSums(Zu)
    for (j in seq_along(u))
      uh[j, ] <- uh[j, ] + tabulate(as.integer(pmin(Zu[j, ]^2 / BW, NB - 1)) + 1L, NB)

    for (i in seq_along(S)) {
      s   <- S[[i]]
      Tst <- colSums(forwardsolve(s$L, Zu[pos[[i]], , drop = FALSE])^2)
      lp  <- -pchisq(Tst, df = s$k, lower.tail = FALSE, log.p = TRUE) / log(10)
      tb  <- as.integer(pmin(pmax(Tst, 0) / TBW, TNB - 1)) + 1L
      lb  <- as.integer(pmin(lp / LBW, LNB - 1)) + 1L
      s$cnt <- s$cnt + vapply(lthr, function(x) sum(lp >= x), 0)
      s$th  <- s$th + tabulate(tb, TNB)
      s$lh  <- s$lh + tabulate(lb, LNB)
      for (cl in 1:2) {
        w <- ci == cl
        lpw <- lp[w]
        s$nc[cl] <- s$nc[cl] + sum(w)
        s$ccnt[cl, ] <- s$ccnt[cl, ] + vapply(lthr, function(x) sum(lpw >= x), 0)
        s$cth[cl, ]  <- s$cth[cl, ] + tabulate(tb[w], TNB)
        s$clh[cl, ]  <- s$clh[cl, ] + tabulate(lb[w], LNB)
      }
      h <- which(lp >= -log10(hit_p))
      s$n_hits_all <- s$n_hits_all + length(h)
      if (length(h) && s$n_hits < max_hits) {
        h <- head(h, max_hits - s$n_hits)
        s$hits[[length(s$hits) + 1L]] <- data.table(index = vidx[h], T = Tst[h], neglog10p = lp[h])
        s$n_hits <- s$n_hits + length(h)
      }
      S[[i]] <- s
    }
  }
  close(con)
  n_drop <- nv - n_ok
  lgc_u  <- apply(uh, 1, hist_median, n = n_ok, w = BW) / qchisq(0.5, 1)
  Cu     <- Szz / n_ok

  for (i in seq_along(S)) {
    s <- S[[i]]; kk <- s$k; rdir <- rdirs[[as.character(kk)]]
    dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(file.path(rdir, "qc_summary.txt"))) file.remove(file.path(rdir, "qc_summary.txt"))

    # z-covariance against the Rhat block: E[T]/k ~ 1 when they match (05's QC).
    C    <- Cu[pos[[i]], pos[[i]], drop = FALSE]
    eR   <- eigen(s$Rs, symmetric = TRUE)
    Rmh  <- eR$vectors %*% (t(eR$vectors) / sqrt(eR$values))
    M    <- Rmh %*% C %*% Rmh
    lam  <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    pred <- sat_tail(lam, qchisq(5e-8, df = kk, lower.tail = FALSE))
    dcor <- if (kk > 1) max(abs(cov2cor(C) - s$Rs)[lower.tri(s$Rs)]) else 0
    lgc  <- lgc_u[pos[[i]]]
    lgc_T <- hist_median(s$th, n_ok, TBW) / qchisq(0.5, df = kk)

    write_hist(s$th, s$lh, file.path(rdir, "direct_hist.txt.gz"))
    for (cl in 1:2) {
      write_hist(s$cth[cl, ], s$clh[cl, ], file.path(rdir, sprintf("direct_hist_%s.txt.gz", CLASSES[cl])))
      write_counts(s$ccnt[cl, ], s$nc[cl], file.path(rdir, sprintf("direct_counts_%s.txt", CLASSES[cl])))
    }
    write_counts(s$cnt, n_ok, file.path(rdir, "direct_counts.txt"))
    if (length(s$hits)) fwrite(rbindlist(s$hits), file.path(rdir, "direct_hits.txt"), sep = "\t")
    fwrite(data.table(trait = s$sel, mean_z = s1[pos[[i]]] / n_ok, mean_z2 = diag(C), lambda_gc = lgc),
           file.path(rdir, "qc_traits.txt"), sep = "\t")

    summ <- data.table(
      field = c("scenario", "rhat_prefix", "replicate", "k_tested", "traits", "subset_seed",
                "n_variant_binary", "n_variant_written", "n_variant_dropped_nonfinite",
                "Rhat_lambda_min", "Rhat_condition", "tr_Rhat_inv",
                "zcov_vs_Rhat_lambda_min", "zcov_vs_Rhat_lambda_max", "E_T_over_k",
                "predicted_type1_ratio_5e-8", "max_abs_corZ_minus_Rhat",
                "lambda_gc_median", "lambda_gc_min", "lambda_gc_max",
                "jagwas_median_T", "jagwas_lambda_gc",
                "direct_obs_5e-8", "direct_exp_5e-8", "direct_ratio_5e-8",
                "direct_hits_recorded", "direct_hits_total", "minutes"),
      value = c(opt$scenario, opt$rhat_prefix, r, kk, paste(s$sel, collapse = ","), seed_base + r,
                format(nv, scientific = FALSE), format(n_ok, scientific = FALSE),
                format(n_drop, scientific = FALSE),
                min(s$ev), max(s$ev) / min(s$ev), sum(1 / s$ev),
                min(lam), max(lam), sum(lam) / kk, pred / 5e-8, dcor,
                median(lgc, na.rm = TRUE), min(lgc, na.rm = TRUE), max(lgc, na.rm = TRUE),
                lgc_T * qchisq(0.5, df = kk), lgc_T,
                s$cnt[i58], n_ok * THR[i58], s$cnt[i58] / (n_ok * THR[i58]), s$n_hits, s$n_hits_all,
                round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)))
    setnames(summ, "field", "key")
    fwrite(summ, file.path(rdir, "qc_summary.txt"), sep = "\t")   # last: marks completion

    check <- ""
    if (kk == k_all) {
      d05 <- compare_counts(s$cnt, file.path(opt$inbase, out_base, tag, "direct_counts.txt"))
      d10 <- max(vapply(1:2, function(cl) compare_counts(s$ccnt[cl, ],
                   file.path(opt$inbase, out_base, tag, sprintf("direct_counts_%s.txt", CLASSES[cl]))), 0))
      check <- sprintf(" | vs 05: %s, vs 10: %s",
                       if (is.na(d05)) "no file" else if (d05 <= 1) "matches 05" else sprintf("** differs by %.0f **", d05),
                       if (is.na(d10)) "no file" else if (d10 <= 1) "matches 10" else sprintf("** differs by %.0f **", d10))
    }
    msg(sprintf("%s k=%3d: E[T]/k %.4f, lambda_GC %.4f, P < 5e-8 %.0f vs %.2f expected (common %.0f, rare %.0f)%s",
                tag, kk, sum(lam) / kk, lgc_T, s$cnt[i58], n_ok * THR[i58],
                s$ccnt[1, i58], s$ccnt[2, i58], check))
    if (opt$rhat_prefix != "Rhat" || opt$scenario != "S4_normal_miss1pct")
      if (abs(sum(lam) / kk - 1) > 0.02) msg("   ** E[T]/k is far from 1: the Rhat block does not match the z-scores.")
  }
  msg(sprintf("%s: %.1f min", tag, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
msg("Done.")
