#!/usr/bin/env Rscript
# pw01_define_scenarios.R
#
# Defines the power-simulation design and writes <root>/design/scenarios.rds.
# Every random choice lives here -- causal and null allele frequencies, which
# traits each causal variant affects, random effect signs -- so pw02 only
# follows instructions and every method is scored against the same truth.
#
# DATA SETS. Each scenario has --n_ds simulated data sets (default 1000). Data
# set j has n = 22,250 unrelated individuals, one causal variant, --K null
# variants and 128 phenotypes:
#
#     Y = g_j beta_j' + E_j,    rows of E_j ~ MVN(0, R),    g_j standardized
#
# One causal variant per data set keeps every variant's signal exactly as
# specified. Putting all 1000 causal variants into one phenotype panel would,
# for concordant architectures, add a genetic covariance of order
# 1000 * beta beta' to the phenotypes and change the very correlation
# structure being studied. The noise E_j, the genotypes and the allele
# frequencies depend only on j, so every scenario is run on the same data sets
# and differences between scenarios are not Monte Carlo noise.
#
# SIGNAL. As Reviewer 1 asks, the total multivariate noncentrality is held
# fixed across architectures:
#     --signal=mahalanobis   n beta' R^-1 beta = ncp   (default: JAGWAS's own NCP)
#     --signal=marginal      n sum(beta^2)     = ncp   (sum of single-trait NCPs)
# With --ncp=0 the NCP is chosen so that JAGWAS has --target_power at alpha.
# Under the Mahalanobis scaling JAGWAS then has the same power in every
# scenario by construction, which doubles as a check on the whole pipeline;
# the result is how the other methods move around that line.
#
# ARCHITECTURES
#   m = 1, 2, 4, ..., 128 affected traits, drawn at random for each causal
#       variant ("random"); and
#   one correlated block: the --block_size traits most tightly correlated in R
#       (greedy on |r|), the same block for every causal variant.
# DIRECTIONS for m >= 2, equal magnitudes across the affected traits:
#   concordant   signs of the leading eigenvector of R on those traits:
#                effects agree with the phenotypic correlation
#   discordant   signs of the trailing eigenvector: effects run against it
#   mixed        independent random signs
#
# ALLELE FREQUENCIES are drawn from the real MAF >= --maf_min spectrum in
# variant_af.tsv (uniform on [maf_min, 0.5] if that file is unavailable).
#
# Usage:
#   Rscript pw01_define_scenarios.R --corr=<128 x 128 correlation matrix>

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

defaults <- list(
  root         = POWER_ROOT,
  corr         = file.path(POWER_ROOT, "T1_pooled_correlation_matrix.txt"),
  af           = "/data484_4/zxie3/bolong_torchgwas_full_tstats_current/shared/variant_af.tsv",
  k            = "128",
  n            = "22250",
  n_ds         = "1000",       # data sets (= causal variants) per scenario
  K            = "200",        # null variants per data set
  m_levels     = "1,2,4,8,16,32,64,128",
  directions   = "concordant,discordant,mixed",
  block_size   = "8",
  signal       = "mahalanobis",
  ncp          = "0",          # 0: set from --target_power
  target_power = "0.5",
  alpha        = "5e-8",
  maf_min      = "0.01",
  seed         = "20261001"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k       <- as_int(opt$k, "k")
n       <- as_int(opt$n, "n")
n_ds    <- as_int(opt$n_ds, "n_ds")
K       <- as_int(opt$K, "K")
alpha   <- as_num(opt$alpha, "alpha")
maf_min <- as_num(opt$maf_min, "maf_min")
seed    <- as_int(opt$seed, "seed")
bsize   <- as_int(opt$block_size, "block_size")
m_lev   <- as.integer(strsplit(opt$m_levels, ",", fixed = TRUE)[[1]])
dirs    <- trimws(strsplit(opt$directions, ",", fixed = TRUE)[[1]])
if (!opt$signal %in% c("mahalanobis", "marginal")) stop("--signal must be mahalanobis or marginal.")
if (any(m_lev < 1 | m_lev > k)) stop("--m_levels must lie in 1..", k)
if (length(setdiff(dirs, c("concordant", "discordant", "mixed"))))
  stop("--directions must be from concordant, discordant, mixed.")

dir.create(file.path(opt$root, "design"), recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$root, "logs", "pw01_define_scenarios.log"))

# ------------------------------------------------------ correlation and NCP --

R <- read_corr_matrix(opt$corr, k = k)
ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
if (min(ev) <= 0) stop("R is not positive definite.")
Rinv <- solve(R)
meff <- li_ji_meff(R)
qcrit <- qchisq(alpha, df = k, lower.tail = FALSE)
jag_power <- function(ncp) pchisq(qcrit, df = k, ncp = ncp, lower.tail = FALSE)

ncp <- as_num(opt$ncp, "ncp")
if (ncp == 0) {
  if (opt$signal != "mahalanobis")
    stop("--ncp=0 (derived from --target_power) is defined for --signal=mahalanobis only;",
         " give --ncp explicitly.")
  ncp <- uniroot(function(l) jag_power(l) - as_num(opt$target_power, "target_power"),
                 c(1e-6, 5000), tol = 1e-10)$root
}

msg(sprintf("R: %d x %d, lambda_min %.4g, condition %.0f, tr(R^-1) %.0f, Li & Ji M_eff %.2f",
            k, k, min(ev), max(ev) / min(ev), sum(diag(Rinv)), meff))
msg(sprintf("Signal: %s NCP = %.3f per causal variant (JAGWAS analytic power at %g: %.3f)",
            opt$signal, ncp, alpha, if (opt$signal == "mahalanobis") jag_power(ncp) else NA))

# ---------------------------------------------------------------- the block --

find_block <- function(R, size) {
  A <- abs(R); diag(A) <- 0
  S <- as.integer(which(A == max(A), arr.ind = TRUE)[1, ])
  while (length(S) < size) {
    rest <- setdiff(seq_len(nrow(R)), S)
    S <- c(S, rest[which.max(rowMeans(A[rest, S, drop = FALSE]))])
  }
  sort(S)
}
block <- find_block(R, bsize)
off <- function(M) M[upper.tri(M)]
msg(sprintf("Correlated block (%d traits): %s", bsize, paste(block, collapse = ",")))
msg(sprintf("  mean |r| within block %.3f (range %.3f to %.3f); mean |r| over all pairs %.3f",
            mean(abs(off(R[block, block]))), min(off(R[block, block])),
            max(off(R[block, block])), mean(abs(off(R)))))

# ------------------------------------------------------------ effect vectors --

dir_signs <- function(S, dir) {
  m <- length(S)
  if (m == 1) return(1)
  if (dir == "mixed") return(sample(c(-1, 1), m, replace = TRUE))
  e <- eigen(R[S, S, drop = FALSE], symmetric = TRUE)
  v <- if (dir == "concordant") e$vectors[, 1] else e$vectors[, m]
  s <- sign(v); s[s == 0] <- 1
  s * s[1]                      # only the relative signs matter
}

scale_beta <- function(d) {
  q <- if (opt$signal == "mahalanobis") drop(crossprod(d, Rinv %*% d)) else sum(d^2)
  sqrt(ncp / (n * q)) * d
}

scen <- list()
for (m in m_lev) for (d in if (m == 1) "single" else dirs)
  scen[[length(scen) + 1L]] <- data.table(name = sprintf("m%03d_%s", m, d), arch = "random",
                                          m = m, direction = d)
for (d in dirs)
  scen[[length(scen) + 1L]] <- data.table(name = sprintf("block%02d_%s", bsize, d), arch = "block",
                                          m = bsize, direction = d)
scen <- rbindlist(scen)

beta <- list()
for (i in seq_len(nrow(scen))) {
  s <- scen[i]
  set.seed(seed + 1000L * i)
  B <- matrix(0, n_ds, k)
  for (j in seq_len(n_ds)) {
    S <- if (s$arch == "block") block else sort(sample.int(k, s$m))
    dvec <- numeric(k)
    dvec[S] <- dir_signs(S, if (s$direction == "single") "concordant" else s$direction)
    B[j, ] <- scale_beta(dvec)
  }
  beta[[s$name]] <- B
  maha <- n * rowSums((B %*% Rinv) * B)
  scen[i, `:=`(ncp_mahalanobis = mean(maha),
               ncp_marginal_total = mean(n * rowSums(B^2)),
               ncp_marginal_max = mean(n * apply(B^2, 1, max)),
               jagwas_power_analytic = mean(jag_power(maha)))]
}

# ------------------------------------------------------- allele frequencies --

set.seed(seed)
pool <- NULL
if (nzchar(opt$af) && file.exists(opt$af)) {
  pool <- fread(opt$af, select = "minor_allele_frequency", showProgress = FALSE)[[1]]
  pool <- pool[is.finite(pool) & pool >= maf_min & pool <= 0.5]
  msg(sprintf("MAF pool: %d real variants with MAF >= %g (%s)", length(pool), maf_min, opt$af))
} else {
  msg(sprintf("** %s not found: MAF drawn uniformly on [%g, 0.5].", opt$af, maf_min))
}
draw_maf <- function(x) if (length(pool)) sample(pool, x, replace = TRUE) else runif(x, maf_min, 0.5)
maf_causal <- draw_maf(n_ds)
maf_null   <- matrix(draw_maf(n_ds * K), n_ds, K)

# ------------------------------------------------------------------- output --

design <- list(created = Sys.time(), k = k, n = n, n_ds = n_ds, K = K, alpha = alpha,
               signal = opt$signal, ncp = ncp, qcrit = qcrit, meff = meff, R = R,
               block = block, scen = scen, beta = beta,
               maf_causal = maf_causal, maf_null = maf_null,
               seed_data = seed + 7L, opt = opt)
f <- file.path(opt$root, "design", "scenarios.rds")
saveRDS(design, f)
fwrite(scen, file.path(opt$root, "design", "scenarios.txt"), sep = "\t")

msg("")
msg(sprintf("%-22s %4s %-10s %9s %9s %9s %8s", "scenario", "m", "direction", "NCP_maha",
            "NCP_marg", "max_marg", "JAGWAS"))
for (i in seq_len(nrow(scen))) {
  s <- scen[i]
  msg(sprintf("%-22s %4d %-10s %9.1f %9.1f %9.2f %8.3f", s$name, s$m, s$direction,
              s$ncp_mahalanobis, s$ncp_marginal_total, s$ncp_marginal_max,
              s$jagwas_power_analytic))
}
msg("")
msg("NCP_marg is the summed single-trait NCP and max_marg the largest single one --")
msg("the quantity that drives minP. Wrote ", f)
