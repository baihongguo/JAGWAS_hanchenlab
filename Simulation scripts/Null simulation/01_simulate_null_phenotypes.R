#!/usr/bin/env Rscript
# 01_simulate_null_phenotypes.R
#
# Simulates null phenotypes for the JAGWAS type-I-error calibration study.
# 128 traits x n individuals x 257 replicates, under four scenarios:
#
#   S1_normal            multivariate normal with the observed residual
#                        correlation matrix R
#   S2_lognormal         log-normal marginals, moderate skew, with the PEARSON
#                        correlation matched exactly to R
#   S3_cauchy            standard Cauchy marginals through the Gaussian copula
#                        (rank correlation identical to S1; Pearson undefined)
#   S4_normal_miss1pct   S1 with 1% of values missing, independently for each
#                        trait and individual
#
# The phenotypes are independent of genotype by construction, so every SNP is a
# null SNP and the entire genome contributes to the empirical null.
#
# COMMON RANDOM NUMBERS. Within a replicate all four scenarios are built from
# one draw of iid standard normals. S3 is an exact marginal transform of S1 and
# S4 is S1 plus a missingness mask, so differences between scenarios reflect the
# marginal distribution and the missingness, not Monte Carlo noise.
#
# Usage (all four scenarios, all replicates):
#   Rscript 01_simulate_null_phenotypes.R
#
# Usage (SLURM array over replicates):
#   Rscript 01_simulate_null_phenotypes.R --rep_start=$i --rep_end=$i
#
# Usage (one scenario only):
#   Rscript 01_simulate_null_phenotypes.R --scenarios=S3_cauchy

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

# Pinned explicitly so that a collaborator regenerating these phenotypes on a
# different machine gets byte-identical output from the same seeds, regardless
# of their R version's defaults. Mersenne-Twister + Inversion are stable across
# R versions and platforms.
RNGkind("Mersenne-Twister", "Inversion")

defaults <- list(
  corr       = "/HGCNT95FS/ADDLIE/work/Test/new_GWAS/residuals/T1_pooled_correlation_matrix.txt",
  samples    = "/HGCNT95FS/ADDLIE/work/Test/new_GWAS/T1_new_discovery.txt",
  outdir     = "/HGCNT95FS/ADDLIE/work/Test/JAGWAS_null",
  k          = "128",         # number of traits
  nrep       = "257",         # replicates per scenario
  rep_start  = "1",
  rep_end    = "0",           # 0 => nrep
  scenarios  = "all",         # comma-separated subset of SCENARIOS, or "all"
  sigma_log  = "0.5",         # log-normal SD on the log scale; skewness ~ 1.75
  miss_rate  = "0.01",        # S4 per-trait, per-individual missingness
  apply_rint = "false",       # rank-inverse-normal transform before writing
  seed_base  = "20260902",
  digits     = "6",           # significant digits written out
  header     = "false",       # GCTA/fastGWA phenotype files carry NO header
  na_string  = "NA",          # GCTA accepts "NA" or "-9"
  threads    = "0",           # data.table threads; 0 => leave default
  overwrite  = "false"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

k          <- as_int(opt$k, "k")
nrep       <- as_int(opt$nrep, "nrep")
rep_start  <- as_int(opt$rep_start, "rep_start")
rep_end    <- as_int(opt$rep_end, "rep_end"); if (rep_end == 0L) rep_end <- nrep
sigma_log  <- as_num(opt$sigma_log, "sigma_log")
miss_rate  <- as_num(opt$miss_rate, "miss_rate")
apply_rint <- as_bool(opt$apply_rint, "apply_rint")
seed_base  <- as_int(opt$seed_base, "seed_base")
digits     <- as_int(opt$digits, "digits")
write_hdr  <- as_bool(opt$header, "header")
overwrite  <- as_bool(opt$overwrite, "overwrite")
threads    <- as_int(opt$threads, "threads")
if (threads > 0) setDTthreads(threads)

scen <- if (identical(opt$scenarios, "all")) SCENARIOS else
  trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
bad <- setdiff(scen, ALL_SCENARIOS)
if (length(bad)) stop("Unknown scenario(s): ", paste(bad, collapse = ", "),
                      ". Valid: ", paste(ALL_SCENARIOS, collapse = ", "))
if (rep_start < 1L || rep_end > nrep || rep_start > rep_end)
  stop("Invalid replicate range: ", rep_start, "..", rep_end, " (nrep = ", nrep, ")")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$outdir, "logs",
                   sprintf("01_simulate_rep%04d-%04d.log", rep_start, rep_end)))
on.exit(log_close(), add = TRUE)

msg("JAGWAS null phenotype simulation")
msg("  scenarios   : ", paste(scen, collapse = ", "))
msg("  replicates  : ", rep_start, "..", rep_end, " of ", nrep)
msg("  traits (k)  : ", k)
msg("  sigma_log   : ", sigma_log)
msg("  miss_rate   : ", miss_rate)
msg("  apply_rint  : ", apply_rint)
msg("  seed_base   : ", seed_base)
msg("  outdir      : ", opt$outdir)

# ------------------------------------------------------------------- inputs --

msg("Reading sample IDs: ", opt$samples)
ids <- fread(opt$samples, header = TRUE, colClasses = "character",
             select = c("FID", "IID"), data.table = TRUE)
if (!all(c("FID", "IID") %in% names(ids)))
  stop("Sample file must have FID and IID columns.")
n <- nrow(ids)
if (anyDuplicated(ids$IID)) stop("Duplicate IIDs in ", opt$samples)
msg(sprintf("  n = %d individuals", n))

msg("Reading correlation matrix: ", opt$corr)
R <- read_corr_matrix(opt$corr, k = k)
lmin <- min_eigen(R)
msg(sprintf("  lambda_min = %.6g", lmin))
if (lmin <= 1e-10) {
  msg("  R is (near-)singular; projecting to the nearest positive-definite")
  msg("  correlation matrix for SAMPLING ONLY. Keep testing with the original R.")
  R_samp <- project_to_corr(R)
  msg(sprintf("  max |R_sampling - R| = %.3g", max(abs(R_samp - R))))
} else {
  R_samp <- R
}
U_norm <- chol(R_samp)                      # upper triangular, t(U) %*% U = R

# --------------------------------------- latent correlation for the log-normal --

rho_ln <- lognormal_latent_rho(R, sigma_log)
if (is.null(rho_ln) || any(!is.finite(rho_ln)) || min(rho_ln) < -1) {
  stop(sprintf(paste0("sigma_log = %.3f is too large: the target Pearson ",
                      "correlations are unattainable with log-normal marginals.\n",
                      "  Run 00_check_correlation_matrix.R for the maximum feasible value."),
               sigma_log))
}
diag(rho_ln) <- 1
rho_ln <- (rho_ln + t(rho_ln)) / 2
if (min_eigen(rho_ln) <= 1e-10) {
  msg("  Log-normal latent matrix is not PD; projecting to nearest PD correlation matrix.")
  rho_ln_pd <- project_to_corr(rho_ln)
  msg(sprintf("  max |rho_pd - rho| = %.3g", max(abs(rho_ln_pd - rho_ln))))
  rho_ln <- rho_ln_pd
}
U_ln <- chol(rho_ln)
msg(sprintf("Log-normal latent correlation: range [%.4f, %.4f], max shift from R = %.4f",
            min(rho_ln[lower.tri(rho_ln)]), max(rho_ln[lower.tri(rho_ln)]),
            max(abs(rho_ln - R))))
msg(sprintf("Log-normal marginal skewness = %.3f",
            (exp(sigma_log^2) + 2) * sqrt(exp(sigma_log^2) - 1)))

# --------------------------------------------------------- output scaffolding --

for (s in scen) dir.create(file.path(opt$outdir, "pheno", s),
                           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "qc"), recursive = TRUE, showWarnings = FALSE)

trait_names <- sprintf("QT%03d", seq_len(k))
colmap <- data.table(mpheno = seq_len(k), file_column = seq_len(k) + 2L,
                     trait = trait_names)
fwrite(colmap, file.path(opt$outdir, "qc", "phenotype_columns.txt"), sep = "\t")

# Latent matrices actually used, for the record.
fwrite(as.data.table(R_samp), file.path(opt$outdir, "qc", "R_sampling.txt"),
       sep = "\t", col.names = FALSE)
fwrite(as.data.table(rho_ln), file.path(opt$outdir, "qc", "R_latent_lognormal.txt"),
       sep = "\t", col.names = FALSE)

pheno_path <- function(s, r) file.path(opt$outdir, "pheno", s, sprintf("rep%03d.pheno", r))

write_pheno <- function(X, s, r) {
  f <- pheno_path(s, r)
  if (file.exists(f) && !overwrite) { msg("  skip (exists): ", f); return(invisible(FALSE)) }
  dt <- data.table(FID = ids$FID, IID = ids$IID)
  Xs <- signif(X, digits)
  dt <- cbind(dt, as.data.table(Xs))
  setnames(dt, c("FID", "IID", trait_names))
  fwrite(dt, f, sep = "\t", na = opt$na_string, quote = FALSE,
         col.names = write_hdr)
  invisible(TRUE)
}

# ------------------------------------------------------------------ the loop --

manifest <- vector("list", 0L)
t0 <- Sys.time()

for (r in rep_start:rep_end) {

  # One draw of iid normals per replicate, shared by every scenario.
  set.seed(seed_base + r)
  E <- matrix(rnorm(n * k), nrow = n, ncol = k)

  Z_norm <- E %*% U_norm                    # cov = R_samp
  X1 <- NULL

  for (s in scen) {
    X <- switch(s,
      "S1_normal" = {
        if (is.null(X1)) X1 <- Z_norm
        X1
      },
      "S2_lognormal" = {
        # Separate latent matrix: exp() distorts correlation, and the latent rho
        # above is chosen so the OUTPUT Pearson correlation equals R exactly.
        exp(sigma_log * (E %*% U_ln))
      },
      "S3_cauchy" = {
        # Pure marginal transform of S1 through the Gaussian copula: the rank
        # correlation is identical to S1, so any difference in calibration is
        # attributable to the tails alone.
        cauchy_from_z(Z_norm)
      },
      "S4_normal_miss1pct" = {
        if (is.null(X1)) X1 <- Z_norm
        Xm <- X1
        set.seed(seed_base + 500000L + r)   # mask stream kept separate from E
        Xm[matrix(runif(n * k) < miss_rate, nrow = n, ncol = k)] <- NA_real_
        Xm
      },
      # Opt-in S5_t<df>: same Gaussian copula as S3, so rank correlation is
      # again identical to S1 and only the tail weight changes.
      t_from_z(Z_norm, df = t_df_from_name(s))
    )

    # Guard: the Cauchy transform would return +/-Inf if |z| exceeded ~38, which
    # cannot happen with 7.6e8 standard normal draws but would silently poison
    # the GWAS if it did. Fail loudly rather than write Inf.
    bad <- !is.finite(X) & !is.na(X)
    if (any(bad)) {
      stop(sprintf("Non-finite values in %s replicate %d (%d cells). Investigate before proceeding.",
                   s, r, sum(bad)))
    }

    if (apply_rint) X <- rint(X)

    if (write_pheno(X, s, r)) {
      manifest[[length(manifest) + 1L]] <- data.table(
        scenario = s, replicate = r, file = pheno_path(s, r),
        n = n, k = k, seed_E = seed_base + r,
        seed_mask = if (s == "S4_normal_miss1pct") seed_base + 500000L + r else NA_integer_,
        n_missing = sum(is.na(X))
      )
    }
  }

  if (r %% 10L == 0L || r == rep_end) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    done <- r - rep_start + 1L
    msg(sprintf("  replicate %d/%d  (%.1f min elapsed, ~%.1f min remaining)",
                r, rep_end, el, el / done * (rep_end - r)))
  }
}

if (length(manifest)) {
  mf <- rbindlist(manifest)
  mfile <- file.path(opt$outdir, "qc",
                     sprintf("manifest_rep%04d-%04d.txt", rep_start, rep_end))
  fwrite(mf, mfile, sep = "\t")
  msg("Wrote manifest: ", mfile)
  msg(sprintf("Generated %d phenotype files.", nrow(mf)))
}

msg(sprintf("Total elapsed: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
msg("Done.")
