#!/usr/bin/env Rscript
# 04_estimate_scenario_R.R
#
# Writes the correlation matrix JAGWAS uses (--cor_matrix) for every replicate:
#
#   <outdir>/<scenario>/Rhat_rep001.txt ... Rhat_rep257.txt
#
# HOW RHAT IS ESTIMATED. The Pearson correlation of that replicate's own 128
# simulated phenotypes, among the individuals torchGWAS analysed (--keep).
# S4 has missing values, so each pair of traits uses the individuals observed
# for both (pairwise-complete). This mirrors the real pipeline, where R is the
# residual correlation of the phenotypes that went into the GWAS, and it is the
# right matrix: with the phenotypes fixed and the genotypes varying across SNPs,
# the z-scores' cross-trait covariance is exactly this sample correlation. The
# population R, or a matrix from a different set of individuals, is not -- and
# R's condition number (4791) amplifies any difference.
#
# S3_cauchy has no population correlation (no finite variance), but for a given
# replicate the sample correlation is still exactly the z-scores' covariance.
# It just varies a lot between replicates; the summary shows by how much.
#
# OVERLAP ADJUSTMENT (--overlap_adjust). With phenotype-specific missingness
# (S4) each trait's GWAS uses a different subset of people, so the z-scores'
# correlation is r_ij * n_ij / sqrt(n_i n_j), not the phenotypic r_ij -- the
# cause of S4's inflation. The adjusted matrix is computed in one step: each
# trait standardized on its own observed values, missing entries set to 0,
# Ytil = that matrix. Only people observed for both traits contribute to an
# off-diagonal entry of Ytil'Ytil, which builds in the overlap factor exactly,
# and the result is a valid (positive semi-definite) matrix by construction.
#   none  Pearson correlation, pairwise-complete when values are missing
#         (the naive matrix)                               -> Rhat_rep###.txt
#   drop  cov2cor(Ytil'Ytil): GWAS that drop each trait's
#         missing individuals (fastGWA, PLINK, GMMAT)      -> Rhat_adjdrop_rep###.txt
#   fill  Ytil'Ytil / N: GWAS that fill missing values and
#         scale by the full N; diagonal n_i / N            -> Rhat_adjfill_rep###.txt
# The two adjusted versions differ only on the diagonal (1 versus n_i/N), and
# which one is right depends on how the GWAS software handled missing values.
# Write both; 05's QC (E[T]/k ~ 1.00) shows which matches the z-scores. With
# complete data all three are identical.
#
# FORMAT. Bare numbers, space-separated, no header or row names, %.17g so every
# double round-trips exactly. JAGWAS loads the matrix with Armadillo's
# auto_detect, which cannot parse labels, and it reports a failed load without
# stopping. Row/column i = phenotype column QT00i = QT{i-1}.fastGWA.gz.
#
# Usage:
#   Rscript 04_estimate_scenario_R.R --pheno_dir=<dir holding S1_normal/, ...> \
#     --keep=<IDs torchGWAS analysed>
#   Rscript 04_estimate_scenario_R.R ... --scenarios=S4_normal_miss1pct --overlap_adjust=fill

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  pheno_dir = "/HGCNT95FS/ADDLIE/work/Test/JAGWAS_null/pheno",   # <pheno_dir>/<scenario>/rep###.pheno
  outdir    = "/HGCNT95FS/ADDLIE/work/Test/JAGWAS_null/JAGWAS_null_input",
  keep      = "",            # IDs torchGWAS analysed (FID IID, or IID only)
  n         = "22250",       # expected number of analysed individuals, for a check
  k         = "128",
  scenarios = "all",
  rep_start = "1",
  rep_end   = "257",
  header    = "false",       # phenotype files have no header
  overlap_adjust = "none",   # none | drop | fill
  overwrite = "false"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k         <- as_int(opt$k, "k")
n_expect  <- as_num(opt$n, "n")
rep_start <- as_int(opt$rep_start, "rep_start")
rep_end   <- as_int(opt$rep_end, "rep_end")
has_hdr   <- as_bool(opt$header, "header")
overwrite <- as_bool(opt$overwrite, "overwrite")
adj       <- opt$overlap_adjust
if (!adj %in% c("none", "drop", "fill")) stop("--overlap_adjust must be none, drop or fill.")
prefix    <- if (adj == "none") "Rhat" else paste0("Rhat_adj", adj)

scen <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
bad <- setdiff(scen, ALL_SCENARIOS)
if (length(bad)) stop("Unknown scenario(s): ", paste(bad, collapse = ", "))

log_open(file.path(opt$outdir, "logs",
                   sprintf("04_%s_rep%03d-%03d.log", prefix, rep_start, rep_end)))
on.exit(log_close(), add = TRUE)

msg("Per-replicate Rhat from simulated phenotypes (", prefix, ", overlap adjustment: ", adj, ")")
msg("  phenotypes : ", opt$pheno_dir)
msg("  scenarios  : ", paste(scen, collapse = ", "))
msg("  replicates : ", rep_start, "..", rep_end)

keep_ids <- NULL
if (nzchar(opt$keep)) {
  kp <- fread(opt$keep, header = FALSE, colClasses = "character", showProgress = FALSE)
  if (toupper(sub("^#", "", kp[[1]][1])) %in% c("FID", "IID")) kp <- kp[-1]
  keep_ids <- unique(if (ncol(kp) >= 2) kp[[2]] else kp[[1]])
  msg(sprintf("  keep list  : %d individuals (%s)", length(keep_ids), opt$keep))
} else {
  msg("  ** No --keep list: using every phenotyped individual, not the ", n_expect,
      " torchGWAS analysed.")
  msg("     The z-scores' correlation is the sample correlation of the ANALYSED")
  msg("     individuals, so this is a small but avoidable mismatch.")
}

# Armadillo raw_ascii.
write_matrix_plain <- function(M, out) {
  writeLines(apply(M, 1, function(x) paste(sprintf("%.17g", x), collapse = " ")), out)
}

summ <- vector("list", 0L)
warned_n <- FALSE
t0 <- Sys.time()

for (s in scen) {
  for (r in rep_start:rep_end) {
    fin  <- file.path(opt$pheno_dir, s, sprintf("rep%03d.pheno", r))
    sdir <- file.path(opt$outdir, s)
    fout <- file.path(sdir, sprintf("%s_rep%03d.txt", prefix, r))
    if (!file.exists(fin)) { msg("  MISSING: ", fin); next }
    if (file.exists(fout) && !overwrite) next

    dt <- fread(fin, header = has_hdr, na.strings = c("NA", "-9"),
                colClasses = list(character = 1:2), showProgress = FALSE)
    if (ncol(dt) != k + 2L)
      stop(sprintf("%s has %d columns, expected %d.", fin, ncol(dt), k + 2L))

    sel <- if (is.null(keep_ids)) rep(TRUE, nrow(dt)) else dt[[2]] %in% keep_ids
    if (!is.null(keep_ids) && sum(sel) != n_expect && !warned_n) {
      msg(sprintf("  ** %d phenotyped individuals are in the keep list; expected %.0f.",
                  sum(sel), n_expect))
      warned_n <- TRUE
    }
    X <- as.matrix(dt[sel, -(1:2), with = FALSE])
    storage.mode(X) <- "double"

    if (adj == "none") {
      Rhat <- if (anyNA(X)) cor(X, use = "pairwise.complete.obs") else cor(X)
      Rhat <- (Rhat + t(Rhat)) / 2
      diag(Rhat) <- 1
    } else {
      # Each trait standardized on its own observed values, missing set to 0,
      # so that column t has sum of squares n_t and only people observed for
      # both traits contribute to an off-diagonal entry of the cross-product.
      Yt <- X
      for (tt in seq_len(k)) {
        o <- !is.na(X[, tt])
        v <- X[o, tt] - mean(X[o, tt])
        Yt[o, tt] <- v / sqrt(mean(v^2))
        Yt[!o, tt] <- 0
      }
      G <- crossprod(Yt)
      Rhat <- if (adj == "drop") cov2cor(G) else G / nrow(Yt)
      Rhat <- (Rhat + t(Rhat)) / 2
    }
    ev <- eigen(Rhat, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) msg(sprintf("  ** %s rep%03d: Rhat is not positive definite.", s, r))

    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    write_matrix_plain(Rhat, fout)

    n_obs <- colSums(!is.na(X))
    summ[[length(summ) + 1L]] <- data.table(
      scenario = s, replicate = r, n_used = nrow(X),
      n_obs_min = min(n_obs), n_obs_max = max(n_obs),
      lambda_min = min(ev), lambda_max = max(ev),
      condition = max(ev) / min(ev),
      tr_Rinv = if (min(ev) > 0) sum(1 / ev) else NA_real_)
  }
  msg(sprintf("  %s done (%.1f min elapsed)", s,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

if (length(summ)) {
  out <- rbindlist(summ)
  dir.create(file.path(opt$outdir, "qc"), recursive = TRUE, showWarnings = FALSE)
  f <- file.path(opt$outdir, "qc", sprintf("%s_summary_rep%03d-%03d.txt", prefix, rep_start, rep_end))
  fwrite(out, f, sep = "\t")
  msg("")
  for (s in unique(out$scenario)) {
    o <- out[scenario == s]
    msg(sprintf("%-20s n_used %d  lambda_min %.5g..%.5g  tr(Rinv) %.0f..%.0f  (%d reps)",
                s, o$n_used[1], min(o$lambda_min), max(o$lambda_min),
                min(o$tr_Rinv, na.rm = TRUE), max(o$tr_Rinv, na.rm = TRUE), nrow(o)))
  }
  msg("Wrote ", f)
}
msg(sprintf("Total %.1f min. Done.", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
