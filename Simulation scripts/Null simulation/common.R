# common.R -- shared helpers for the JAGWAS null-calibration simulation.
# Sourced by 00_check_correlation_matrix.R, 01_simulate_null_phenotypes.R,
# 02_validate_simulation.R.

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. install.packages('data.table')")
  }
  library(data.table)
})

# ---------------------------------------------------------------- arguments --

parse_args <- function(args, defaults) {
  out <- defaults
  for (a in args) {
    if (!grepl("^--[^=]+=", a)) {
      stop("Malformed argument '", a, "'. Use --key=value.", call. = FALSE)
    }
    k <- sub("^--([^=]+)=.*$", "\\1", a)
    v <- sub("^--[^=]+=", "", a)
    if (!k %in% names(defaults)) {
      stop("Unknown option --", k, ". Valid options: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    }
    out[[k]] <- v
  }
  out
}

as_num <- function(x, what) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) stop("Option --", what, " must be numeric, got '", x, "'", call. = FALSE)
  v
}

as_int <- function(x, what) as.integer(as_num(x, what))

as_bool <- function(x, what) {
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (x %in% c("false", "f", "no", "n", "0")) return(FALSE)
  stop("Option --", what, " must be true/false, got '", x, "'", call. = FALSE)
}

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
}

# ------------------------------------------------------------------ logging --

.LOGCON <- NULL

log_open <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  .LOGCON <<- file(path, open = "wt")
}

log_close <- function() {
  if (!is.null(.LOGCON)) { close(.LOGCON); .LOGCON <<- NULL }
}

msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(line, "\n", sep = "")
  if (!is.null(.LOGCON)) { writeLines(line, .LOGCON); flush(.LOGCON) }
  invisible(NULL)
}

# ----------------------------------------------------------- matrix reading --

# Reads a square correlation matrix, tolerating: header row or not, row-name
# column or not. Validates dimension, symmetry and unit diagonal.
read_corr_matrix <- function(path, k = 128L, tol = 1e-6) {
  if (!file.exists(path)) stop("Correlation matrix not found: ", path, call. = FALSE)
  df <- data.table::fread(path, header = "auto", data.table = FALSE)

  # Drop a leading row-name column if present.
  if (!is.numeric(df[[1]])) df <- df[, -1, drop = FALSE]

  M <- as.matrix(df)
  storage.mode(M) <- "double"
  dimnames(M) <- NULL

  if (nrow(M) != k || ncol(M) != k) {
    stop(sprintf(paste0("Expected a %d x %d matrix in %s but parsed %d x %d.\n",
                        "  Check for stray header/row-name columns, or pass --k=<dim>."),
                 k, k, path, nrow(M), ncol(M)), call. = FALSE)
  }
  if (anyNA(M)) stop("Correlation matrix contains NA.", call. = FALSE)

  asym <- max(abs(M - t(M)))
  if (asym > tol) {
    warning(sprintf("Matrix is asymmetric (max |R - R'| = %.3g); symmetrising.", asym))
  }
  M <- (M + t(M)) / 2

  ddev <- max(abs(diag(M) - 1))
  if (ddev > tol) {
    warning(sprintf("Diagonal departs from 1 by up to %.3g; rescaling to a correlation matrix.", ddev))
    d <- sqrt(diag(M)); M <- M / outer(d, d)
  }
  diag(M) <- 1
  attr(M, "max_asymmetry") <- asym
  attr(M, "max_diag_dev")  <- ddev
  M
}

# --------------------------------------------------- positive-definite fixes --

# Nearest correlation matrix. Uses Matrix::nearPD when available, otherwise an
# eigenvalue-clipping projection (adequate for k = 128).
project_to_corr <- function(M, eps = 1e-8, iter = 100L) {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    np <- Matrix::nearPD(M, corr = TRUE, keepDiag = TRUE, eig.tol = eps,
                         conv.tol = 1e-9, maxit = 1000L)
    return(as.matrix(np$mat))
  }
  for (i in seq_len(iter)) {
    e <- eigen(M, symmetric = TRUE)
    v <- pmax(e$values, eps)
    M <- e$vectors %*% (v * t(e$vectors))
    d <- sqrt(diag(M)); M <- M / outer(d, d)
    M <- (M + t(M)) / 2
    diag(M) <- 1
    if (min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) > eps / 2) break
  }
  M
}

min_eigen <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# --------------------------------------------- log-normal correlation mapping --

# For X_i = exp(sigma * Z_i) with corr(Z_i, Z_j) = rho, the Pearson correlation
# of the log-normals is  r = (exp(rho * s2) - 1) / (exp(s2) - 1),  s2 = sigma^2.
# Inverting gives the latent rho that reproduces a target Pearson r exactly.
lognormal_latent_rho <- function(r, sigma) {
  s2 <- sigma^2
  e  <- exp(s2) - 1
  arg <- 1 + r * e
  if (any(arg <= 0)) return(NULL)          # target unattainable
  log(arg) / s2
}

# Largest sigma for which every target correlation in R remains attainable
# (i.e. the implied latent rho stays >= -1). Returns NA if even tiny sigma fails.
lognormal_max_sigma <- function(rmin, upper = 3, tol = 1e-6) {
  feasible <- function(s) {
    rho <- lognormal_latent_rho(rmin, s)
    !is.null(rho) && is.finite(rho) && rho >= -1
  }
  if (!feasible(tol)) return(NA_real_)
  if (feasible(upper)) return(upper)
  lo <- tol; hi <- upper
  while (hi - lo > 1e-6) {
    mid <- (lo + hi) / 2
    if (feasible(mid)) lo <- mid else hi <- mid
  }
  lo
}

# --------------------------------------------------- marginal transformations --

# Gaussian -> standard Cauchy through the Gaussian copula, computed on the log
# scale so that |z| > 8 does not saturate pnorm() at 0 or 1 and produce +/-Inf.
cauchy_from_z <- function(z) {
  d  <- dim(z)
  zz <- as.vector(z)
  s  <- sign(zz); s[s == 0] <- 1
  lp <- pnorm(-abs(zz), log.p = TRUE)           # log P(Z < -|z|), always <= log(0.5)
  out <- -s * qcauchy(lp, log.p = TRUE)
  if (!is.null(d)) dim(out) <- d                # R's p/q functions may drop dim
  out
}

# Gaussian -> Student t with df degrees of freedom, through the Gaussian copula.
# Same log-scale treatment as cauchy_from_z. t_1 is Cauchy; t_3 has finite
# variance but infinite kurtosis; t_5 has finite kurtosis.
t_from_z <- function(z, df) {
  d  <- dim(z)
  zz <- as.vector(z)
  s  <- sign(zz); s[s == 0] <- 1
  lp <- pnorm(-abs(zz), log.p = TRUE)
  out <- -s * qt(lp, df = df, log.p = TRUE)
  if (!is.null(d)) dim(out) <- d
  out
}

# Rank-based inverse normal transform (Blom), NA-preserving, column-wise.
rint <- function(x) {
  apply(x, 2, function(v) {
    ok <- !is.na(v)
    out <- rep(NA_real_, length(v))
    n <- sum(ok)
    if (n > 0) out[ok] <- qnorm((rank(v[ok], ties.method = "average") - 0.375) / (n + 0.25))
    out
  })
}

# ----------------------------------------------------------------- utilities --

skewness <- function(v) {
  v <- v[is.finite(v)]; n <- length(v)
  if (n < 3) return(NA_real_)
  m <- mean(v); s <- sqrt(mean((v - m)^2))
  if (s == 0) return(NA_real_)
  mean((v - m)^3) / s^3
}

ex_kurtosis <- function(v) {
  v <- v[is.finite(v)]; n <- length(v)
  if (n < 4) return(NA_real_)
  m <- mean(v); s <- sqrt(mean((v - m)^2))
  if (s == 0) return(NA_real_)
  mean((v - m)^4) / s^4 - 3
}

# Li & Ji (2005) effective number of independent tests.
li_ji_meff <- function(R) {
  lam <- abs(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  sum(as.numeric(lam >= 1) + (lam - floor(lam)))
}

# The four scenarios generated by --scenarios=all.
SCENARIOS <- c("S1_normal", "S2_lognormal", "S3_cauchy", "S4_normal_miss1pct")

# Opt-in only: name them explicitly, e.g. --scenarios=S5_t3,S5_t5. These bracket
# the gap between S1 (normal) and S3 (Cauchy = t_1): t_3 has finite variance so
# the CLT still applies, t_5 additionally has finite kurtosis. If t_3 is
# calibrated and Cauchy is not, the failure is located at the finite-variance
# boundary rather than merely observed.
OPTIONAL_SCENARIOS <- c("S5_t3", "S5_t5")
ALL_SCENARIOS <- c(SCENARIOS, OPTIONAL_SCENARIOS)

t_df_from_name <- function(s) as.integer(sub("^S5_t", "", s))

# ----------------------------------------------------------- trait subsets --

# Result folders written by 12_trait_subsets.R end in _k<number> (S1_normal_k16,
# S4_adjfill_k64); the number is the degrees of freedom of the statistic in that
# folder. Folders without the suffix hold the full set of traits (default).
scenario_k <- function(s, default) {
  m <- regmatches(s, regexpr("_k[0-9]+$", s))
  if (length(m)) as.integer(sub("^_k", "", m)) else as.integer(default)
}
scenario_base <- function(s) sub("_k[0-9]+$", "", s)

# The traits tested in replicate r: one random permutation of 1..k_all per
# replicate, seeded by seed_base + r, of which the first k_sub are kept. The
# subsets are therefore nested (the 16 lie inside the 64), and every scenario
# uses the same subset for a given replicate, as with the common random numbers
# the phenotypes are built from.
trait_subset <- function(k_all, k_sub, seed_base, r) {
  if (k_sub >= k_all) return(seq_len(k_all))
  set.seed(seed_base + r)
  sort(sample.int(k_all)[seq_len(k_sub)])
}
