#!/usr/bin/env Rscript
# pz01_zscore_power.R
#
# Power of JAGWAS and Bonferroni minP for 128 phenotypes, computed from simulated
# z-scores, with the total multivariate signal held fixed as Reviewer 1 asks.
#
# WHY Z-SCORES. Both tests use only the k single-phenotype z-scores of a variant
# and the correlation matrix R. For a causal variant z ~ N(delta, R), where
# delta_j = sqrt(2 f (1 - f) N) * beta_j / sd_j is the expected z-score of
# phenotype j (beta_j its per-allele effect). In the null simulation the z-scores
# computed from real genotypes behaved as N(0, R) (E[T]/k = 1.000), so power can
# be computed by drawing z directly, without simulating genotypes or running GWAS.
#
# THE TESTS
#   JAGWAS        T = z' R^-1 z on k df, as in the released software. Its power is
#                 exact: P(chi2_k(lambda) > q_k), lambda = delta' R^-1 delta.
#   JAGWAS (generalized inverse)
#                 T+ = z' R+ z on df = rank(R), with R+ the Moore-Penrose inverse
#                 from the eigendecomposition, eigenvalues below --eig_tol set to 0.
#                 Identical to JAGWAS when no eigenvalue is below the tolerance.
#   minP          Bonferroni: significant when max_j |z_j| > qnorm(1 - alpha/(2k)),
#                 i.e. when the smallest single-phenotype p-value is below alpha/k.
#                 Power by Monte Carlo.
#
# FIXED TOTAL SIGNAL. delta is rescaled so that lambda = delta' R+ delta is the
# same in every setting, at the values that give JAGWAS (k df) --target_power at
# --alpha (80%: 127.7 and 50%: 105.1 for k = 128 at 5e-8). JAGWAS power is then
# the same in every setting with the same df, by construction; what changes is
# minP, through the largest expected single-phenotype z-score max_j |delta_j|.
# By Cauchy-Schwarz max_j delta_j^2 <= lambda, with equality only when delta is
# proportional to a column of R.
#
# COVARIANCE SCENARIOS (--scenarios)
#   empirical     the 128 x 128 T1 UDIP correlation matrix (--corr).
#   nearsingular  64 UDIPs (a random block of --corr, seed --seed) and a near-
#                 duplicate of each, y' = (y + e)/sqrt(1 + v) with var(e) =
#                 --dup_noise: 64 eigenvalues of about dup_noise/2, far below
#                 --eig_tol, so R has numerical rank 64. Effects are consistent:
#                 a near-duplicate has the same effect as its UDIP, so the m
#                 affected phenotypes are m/2 UDIPs and their near-duplicates.
#   signflip      the empirical matrix after reversing the sign of each phenotype
#                 with probability 1/2 (R -> D R D, drawn anew for each causal
#                 set), with effects of equal sign in the new coding, i.e. mixed
#                 directions on the original coding. If the effects were flipped
#                 together with the phenotypes (a pure recoding), both tests would
#                 give exactly the scenario-1 results; the script checks this.
#                 Labelled "Mixed effect directions" in the figure.
# EFFECT PATTERNS (--patterns), for m directly affected phenotypes (--m_levels)
# drawn at random for each causal set:
#   direct        delta_j equal on the m phenotypes and 0 on the others.
#   correlated    delta = R w, with w equal on the m phenotypes and 0 elsewhere:
#                 the other phenotypes respond through their correlation with
#                 the affected ones.
# In the eigenvector coordinates of R (eigenvalues mu_i, eigenvectors v_i), the
# share of lambda on eigenvector i is proportional to (v_i'w)^2 / mu_i for direct
# effects, which therefore load on the low-variance eigenvectors, and to
# mu_i (v_i'w)^2 for correlated effects, which load on the high-variance ones and
# favor minP. The near-singular scenario needs an even m (m/2 UDIPs and their
# near-duplicates), so odd m are skipped there.
#
# OUTPUT (--out)
#   zpower_matrix_summary.tsv  eigenvalues, condition number, rank, tr(R^-1)/k,
#                              Li-Ji M_eff, the number of retained eigenvalues
#                              below each --lambda_share_below cutoff (n_eig_lt_*),
#                              and the mean of T/df that a 1% mismatch (z
#                              covariance 0.99 R + 0.01 I) would give
#   zpower_results.tsv         power by scenario, pattern, m and target, and the
#                              share of lambda on the retained eigenvalues below
#                              each cutoff (lambda_share_eig_lt_*)
#   zpower_null_check.tsv      type I error of the three tests with delta = 0
#   zpower_power_<target>.png  power against m (.pdf with --plot_format=pdf, or
#                              when R has no PNG device, as on headless servers)
#
# Usage:
#   Rscript pz01_zscore_power.R --corr=<128 x 128 correlation matrix> --out=<folder>
# Needs common.R from null_simulation/ in the same folder. With the defaults, about
# 50 minutes with R's reference BLAS (10 with an optimized BLAS); --n_sets=50
# --n_draws=500 --n_null=100000 gives a test run of about 2 minutes.

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
if (!file.exists(file.path(.here, "common.R")))
  stop("common.R not found next to this script: copy it from null_simulation/.")
suppressPackageStartupMessages(source(file.path(.here, "common.R")))

defaults <- list(
  corr         = "/data5/bguo1/JAGWAS_power/T1_pooled_correlation_matrix.txt",
  out          = "/data5/bguo1/JAGWAS_power/zscore",
  k            = "128",
  scenarios    = "empirical,nearsingular,signflip",
  patterns     = "direct,correlated",
  m_levels     = "1,4,8,16,32,64",
  target_power = "0.8,0.5",
  alpha        = "5e-8",
  n_sets       = "1000",     # causal sets (random choices of the affected phenotypes)
  n_draws      = "1000",     # z-score draws per causal set
  n_null       = "1000000",  # null draws per scenario
  eig_tol      = "1e-5",
  dup_noise    = "1e-6",
  lambda_share_below = "0.01,0.1",   # eigenvalue cutoffs for the share of lambda
  seed        = "20261028",
  plot         = "true",
  plot_format  = "auto"      # png, pdf, or auto (png if this R can write it, else pdf)
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k        <- as_int(opt$k, "k")
alpha    <- as_num(opt$alpha, "alpha")
n_sets   <- as_int(opt$n_sets, "n_sets")
n_draws  <- as_int(opt$n_draws, "n_draws")
n_null   <- as_int(opt$n_null, "n_null")
eig_tol  <- as_num(opt$eig_tol, "eig_tol")
dup_v    <- as_num(opt$dup_noise, "dup_noise")
seed     <- as_int(opt$seed, "seed")
m_lev    <- as.integer(strsplit(opt$m_levels, ",", fixed = TRUE)[[1]])
targets  <- as.numeric(strsplit(opt$target_power, ",", fixed = TRUE)[[1]])
scens    <- trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
pats     <- trimws(strsplit(opt$patterns, ",", fixed = TRUE)[[1]])
share_lab <- trimws(strsplit(opt$lambda_share_below, ",", fixed = TRUE)[[1]])
share_tau <- suppressWarnings(as.numeric(share_lab))
if (anyNA(share_tau)) stop("--lambda_share_below must be numbers, e.g. 0.01,0.1.")
share_cols <- paste0("lambda_share_eig_lt_", share_lab)
if (length(setdiff(scens, c("empirical", "nearsingular", "signflip"))))
  stop("--scenarios must be from empirical, nearsingular, signflip.")
if (length(setdiff(pats, c("direct", "correlated"))))
  stop("--patterns must be from direct, correlated.")
if (any(m_lev < 1 | m_lev > k)) stop("--m_levels must lie in 1..", k)
if (k %% 2 != 0 && "nearsingular" %in% scens) stop("The near-singular scenario needs an even k.")
fmt <- if (opt$plot_format == "auto") { if (capabilities("png")) "png" else "pdf" } else opt$plot_format
if (!fmt %in% c("png", "pdf")) stop("--plot_format must be png, pdf or auto.")
if (fmt == "png" && as_bool(opt$plot, "plot") && !capabilities("png"))
  stop("This R cannot write PNG files: use --plot_format=pdf.")

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(opt$out) || file.access(opt$out, 2) != 0)
  stop("Cannot write to --out=", opt$out, ": its parent folder must exist and be writable on this server.")
log_open(file.path(opt$out, "pz01_zscore_power.log"))
on.exit(log_close(), add = TRUE)
msg("pz01_zscore_power: k = ", k, ", alpha = ", alpha, ", ", n_sets, " causal sets x ",
    n_draws, " draws per setting")

# --------------------------------------------------------------- the matrices --

R <- read_corr_matrix(opt$corr, k = k)

# Eigendecomposition with what every test needs: z = delta + E, E = U %*% M with
# U iid N(0, 1), so that the whitened coordinates of z are U + gamma.
decomp <- function(R) {
  e <- eigen(R, symmetric = TRUE)
  mu <- pmax(e$values, 0)
  list(R = R, V = e$vectors, mu = mu, keep = mu > eig_tol,
       M = t(e$vectors %*% diag(sqrt(mu))))
}

half <- k %/% 2
set.seed(seed)
base <- sort(sample.int(k, half))         # UDIPs carried into the near-singular scenario
R64 <- R[base, base]
a <- 1 / sqrt(1 + dup_v)
R_ns <- rbind(cbind(R64, a * R64), cbind(a * R64, (R64 + dup_v * diag(half)) / (1 + dup_v)))

DEC <- list()
if (any(c("empirical", "signflip") %in% scens)) DEC$empirical <- decomp(R)
if ("signflip" %in% scens) DEC$signflip <- DEC$empirical   # same eigenvalues; see below
if ("nearsingular" %in% scens) DEC$nearsingular <- decomp(R_ns)

summ <- rbindlist(lapply(names(DEC), function(s) {
  d <- DEC[[s]]; mu <- d$mu; r <- sum(d$keep)
  inv_all <- if (min(mu) > 0) mean(1 / mu) else Inf
  out <- data.table(scenario = s, k = k, rank = r, eig_tol = eig_tol,
             lambda_max = max(mu), lambda_min = min(mu), condition = max(mu) / min(mu),
             lambda_min_retained = min(mu[d$keep]),
             condition_retained = max(mu) / min(mu[d$keep]),
             mean_inv_eigenvalue = inv_all,              # tr(R^-1)/k, the mean VIF
             mean_inv_eigenvalue_retained = mean(1 / mu[d$keep]),
             li_ji_meff = li_ji_meff(d$R),
             ET_df_mismatch1pct_standard = 0.99 + 0.01 * inv_all,
             ET_df_mismatch1pct_ginv = 0.99 + 0.01 * mean(1 / mu[d$keep]))
  for (i in seq_along(share_tau))
    set(out, j = paste0("n_eig_lt_", share_lab[i]), value = sum(d$keep & mu < share_tau[i]))
  out
}))
fwrite(summ, file.path(opt$out, "zpower_matrix_summary.tsv"), sep = "\t")
for (i in seq_len(nrow(summ))) with(summ[i], msg(sprintf(
  "%-12s rank %3d | eigenvalues %.3g to %.3g, condition %.4g (retained part %.4g) | tr(R^-1)/k %.4g | M_eff %.1f",
  scenario, rank, lambda_min, lambda_max, condition, condition_retained,
  mean_inv_eigenvalue, li_ji_meff)))

# --------------------------------------------------------------- the tests --

cz <- qnorm(alpha / (2 * k), lower.tail = FALSE)          # minP threshold on max |z|
q_k <- qchisq(alpha, df = k, lower.tail = FALSE)
lam_for <- function(target, df) {
  q <- qchisq(alpha, df = df, lower.tail = FALSE)
  uniroot(function(l) pchisq(q, df = df, ncp = l, lower.tail = FALSE) - target,
          c(1e-6, 1e4), tol = 1e-10)$root
}
lambdas <- vapply(targets, lam_for, numeric(1), df = k)
msg(sprintf("minP: max |z| > %.4f. JAGWAS (%d df): T > %.2f. Fixed lambda = %s for JAGWAS power %s",
            cz, k, q_k, paste(sprintf("%.2f", lambdas), collapse = ", "),
            paste(targets, collapse = ", ")))

rowmax_abs <- function(Z) {
  m <- abs(Z[, 1])
  for (j in 2:ncol(Z)) m <- pmax(m, abs(Z[, j]))
  m
}

# --------------------------------------------------------- the null check --

null_rows <- list()
thr <- c(0.05, 1e-3, 1e-4)
for (s in names(DEC)) {
  if (s == "signflip") next                 # same matrix up to recoding: same result
  d <- DEC[[s]]; r <- sum(d$keep)
  set.seed(seed + 7L)
  cnt <- matrix(0, 3, length(thr), dimnames = list(c("jagwas", "jagwas_ginv", "minp"), NULL))
  done <- 0L
  while (done < n_null) {
    n <- min(n_draws * 10L, n_null - done)
    U <- matrix(rnorm(n * k), n, k)
    p1 <- pchisq(rowSums(U^2), df = k, lower.tail = FALSE)
    p2 <- pchisq(rowSums(U[, d$keep, drop = FALSE]^2), df = r, lower.tail = FALSE)
    p3 <- pmin(1, k * 2 * pnorm(-rowmax_abs(U %*% d$M)))
    for (j in seq_along(thr)) cnt[, j] <- cnt[, j] + c(sum(p1 < thr[j]), sum(p2 < thr[j]), sum(p3 < thr[j]))
    done <- done + n
  }
  for (meth in rownames(cnt)) for (j in seq_along(thr)) {
    x <- cnt[meth, j]
    ci <- if (x > 0) qbeta(c(0.025, 0.975), c(x, x + 1), c(n_null - x + 1, n_null - x)) else
      c(0, qbeta(0.975, 1, n_null))
    null_rows[[length(null_rows) + 1L]] <- data.table(
      scenario = s, method = meth, df = if (meth == "jagwas_ginv") r else if (meth == "jagwas") k else NA_integer_,
      alpha = thr[j], n = n_null, rejections = x, rate = x / n_null,
      ci_lo = ci[1], ci_hi = ci[2])
  }
}
null_dt <- rbindlist(null_rows)
fwrite(null_dt, file.path(opt$out, "zpower_null_check.tsv"), sep = "\t")
msg("Null check (empirical type I error, delta = 0):")
for (i in seq_len(nrow(null_dt))) with(null_dt[i], msg(sprintf(
  "   %-12s %-12s alpha %-6g rate %.3g (95%% CI %.3g-%.3g)", scenario, method, alpha, rate, ci_lo, ci_hi)))

# --------------------------------------------------------------- effects --

# Direction of the expected z-scores for one causal set, before scaling.
effect_dir <- function(Rs, w, pattern) if (pattern == "direct") w else drop(Rs %*% w)

# Checks once that flipping the phenotypes together with their effects changes
# nothing: (D R D, D delta, noise E D) gives the same T and max |z| as (R, delta, E).
check_recoding <- function(d) {
  set.seed(seed + 11L)
  D <- sample(c(-1, 1), k, replace = TRUE)
  delta <- drop(d$R %*% c(rep(1, 8), rep(0, k - 8))) * 0.3
  U <- matrix(rnorm(200 * k), 200, k)
  z <- U %*% d$M + rep(delta, each = 200)
  zf <- z * rep(D, each = 200)
  Rf <- d$R * outer(D, D)
  T1 <- rowSums((z %*% solve(d$R)) * z); T2 <- rowSums((zf %*% solve(Rf)) * zf)
  c(max_T_diff = max(abs(T1 - T2)), max_zmax_diff = max(abs(rowmax_abs(z) - rowmax_abs(zf))))
}
if ("empirical" %in% names(DEC)) {
  ck <- check_recoding(DEC$empirical)
  msg(sprintf("Recoding check (phenotypes and effects flipped together): max |T diff| %.2g, max |max|z| diff| %.2g",
              ck[1], ck[2]))
}

# ----------------------------------------------------------- power runs --

res <- list()
t0 <- Sys.time()
for (s in scens) {
  d <- DEC[[s]]; r <- sum(d$keep)
  q_r <- qchisq(alpha, df = r, lower.tail = FALSE)
  for (m in m_lev) {
    if (s == "nearsingular" && m %% 2 != 0) {
      msg(sprintf("%-12s m = %2d skipped: affected phenotypes come in UDIP/near-duplicate pairs", s, m))
      next
    }
    m_base <- if (s == "nearsingular") m %/% 2 else m
    per <- vector("list", n_sets)
    for (j in seq_len(n_sets)) {
      # The same causal set, signs and noise for a given (m, j) in every scenario on
      # 128 phenotypes, so the empirical and sign-flip results are paired. Seeding by
      # m keeps the draws for one m independent of the other levels run.
      set.seed(seed + 100003L * m + j)
      S <- sort(sample.int(if (s == "nearsingular") half else k, m_base))
      D <- sample(c(-1, 1), k, replace = TRUE)
      U <- matrix(rnorm(n_draws * k), n_draws, k)
      E <- U %*% d$M
      w <- numeric(k)
      if (s == "nearsingular") w[c(S, S + half)] <- 1 else w[S] <- 1
      out <- list()
      for (p in pats) {
        if (s == "signflip") {
          # Effects of equal sign on the phenotypes recoded by D (matrix D R D),
          # written on the original coding: delta = D * dir(D R D, w).
          Rf <- d$R * outer(D, D)
          dvec <- D * effect_dir(Rf, w, p)
        } else {
          dvec <- effect_dir(d$R, w, p)
        }
        g_unit <- drop(crossprod(d$V, dvec)) / sqrt(pmax(d$mu, .Machine$double.xmin))
        lam_unit <- sum(g_unit[d$keep]^2)
        sh <- vapply(share_tau, function(tau) sum(g_unit[d$keep & d$mu < tau]^2), numeric(1)) / lam_unit
        for (ti in seq_along(lambdas)) {
          sc <- sqrt(lambdas[ti] / lam_unit)
          delta <- sc * dvec; g <- sc * g_unit
          zmax <- rowmax_abs(E + rep(delta, each = n_draws))
          Tfull <- rowSums((U + rep(g, each = n_draws))^2)
          Tginv <- rowSums((U[, d$keep, drop = FALSE] + rep(g[d$keep], each = n_draws))^2)
          aff <- which(w != 0)
          row <- data.table(
            pattern = p, target = targets[ti], lambda = lambdas[ti],
            jagwas_exact = pchisq(q_k, df = k, ncp = sum(g^2), lower.tail = FALSE),
            jagwas_ginv_exact = pchisq(q_r, df = r, ncp = lambdas[ti], lower.tail = FALSE),
            jagwas_sim = mean(Tfull > q_k), jagwas_ginv_sim = mean(Tginv > q_r),
            minp = mean(zmax > cz),
            max_abs_delta = max(abs(delta)), mean_abs_delta_affected = mean(abs(delta[aff])))
          row[, (share_cols) := as.list(sh)]
          out[[length(out) + 1L]] <- row
        }
      }
      per[[j]] <- rbindlist(out)[, set := j]
    }
    per <- rbindlist(per)
    agg <- per[, c(list(df = k, rank = r,
                     jagwas = mean(jagwas_exact), jagwas_ginv = mean(jagwas_ginv_exact),
                     jagwas_sim = mean(jagwas_sim), jagwas_ginv_sim = mean(jagwas_ginv_sim),
                     minp = mean(minp), minp_se = sd(minp) / sqrt(.N),
                     minp_p10 = quantile(minp, 0.1), minp_p90 = quantile(minp, 0.9),
                     max_abs_delta = mean(max_abs_delta),
                     mean_abs_delta_affected = mean(mean_abs_delta_affected), n_sets = .N),
                   lapply(.SD, mean)),
               by = .(pattern, target, lambda), .SDcols = share_cols]
    agg[, `:=`(scenario = s, m = m)]
    res[[length(res) + 1L]] <- agg
    for (i in seq_len(nrow(agg))) with(agg[i], msg(sprintf(
      "%-12s %-10s m = %2d, power %.2f: JAGWAS %.3f (sim %.3f) | g-inverse %.3f | minP %.3f (SE %.3f) | max E|z| %.2f | lambda share %s",
      s, pattern, m, target, jagwas, jagwas_sim, jagwas_ginv, minp, minp_se, max_abs_delta,
      paste(sprintf("eig < %s: %.2f", share_lab, unlist(agg[i, share_cols, with = FALSE])), collapse = ", "))))
  }
}
res <- rbindlist(res)
setcolorder(res, c("scenario", "pattern", "m", "target", "lambda", "df", "rank"))
fwrite(res, file.path(opt$out, "zpower_results.tsv"), sep = "\t")
msg(sprintf("Power runs: %.1f min. Wrote %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
            file.path(opt$out, "zpower_results.tsv")))

# ----------------------------------------------------------------- figure --

if (as_bool(opt$plot, "plot")) {
  lab_s <- c(empirical = "Empirical UDIP correlation",
             nearsingular = "Nearly singular (64 UDIPs + near-duplicates)",
             signflip = "Mixed effect directions")
  lab_p <- c(direct = "Equal effects",          # labels as in the manuscript (pz02 redraws from results)
             correlated = "Weighted effects")
  col_j <- "black"; col_m <- "#D55E00"
  for (ti in seq_along(targets)) {
    f <- file.path(opt$out, sprintf("zpower_power_%02d.%s", round(100 * targets[ti]), fmt))
    h_in <- (1050 * length(pats) + 360) / 300
    if (fmt == "png") png(f, width = 10, height = h_in, units = "in", res = 300) else
      pdf(f, width = 10, height = h_in)
    par(mfrow = c(length(pats), length(scens)), mar = c(4, 4, 2.6, 0.6), oma = c(0, 0, 4.4, 0),
        mgp = c(2.4, 0.7, 0), las = 1, cex = 0.8)
    for (p in pats) for (s in scens) {
      x <- res[scenario == s & pattern == p & target == targets[ti]][order(m)]
      plot(NA, xlim = range(log2(m_lev)), ylim = c(0, 1), xaxt = "n",
           xlab = "Associated phenotypes (m)", ylab = "Power",
           main = sprintf("%s\n%s", lab_s[[s]], lab_p[[p]]), font.main = 1, cex.main = 0.9)
      axis(1, at = log2(x$m), labels = x$m)         # odd m are absent from the near-singular panels
      abline(h = seq(0, 1, 0.2), col = "grey90", lwd = 0.8)
      lines(log2(x$m), x$jagwas, col = col_j, lwd = 2); points(log2(x$m), x$jagwas, pch = 16, col = col_j)
      if (x$rank[1] < k) {
        lines(log2(x$m), x$jagwas_ginv, col = col_j, lwd = 2, lty = 2)
        points(log2(x$m), x$jagwas_ginv, pch = 17, col = col_j)
      }
      lines(log2(x$m), x$minp, col = col_m, lwd = 2); points(log2(x$m), x$minp, pch = 15, col = col_m)
      ok <- 1.96 * x$minp_se > 0.005         # skip error bars too short to draw
      if (any(ok)) arrows(log2(x$m[ok]), pmax(0, x$minp[ok] - 1.96 * x$minp_se[ok]), log2(x$m[ok]),
                          pmin(1, x$minp[ok] + 1.96 * x$minp_se[ok]),
                          angle = 90, code = 3, length = 0.03, col = col_m)
    }
    top_in <- grconvertY(1, "nic", "ndc")    # top of the panel area, in device units
    line_ndc <- (1 - top_in) / 4.4
    mtext(sprintf("Total signal fixed at lambda = %.1f (JAGWAS power %.0f%% with %d df); alpha = %g",
                  lambdas[ti], 100 * targets[ti], k, alpha), outer = TRUE, line = 2.9, cex = 0.85)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    legend(0.5, top_in + 2.1 * line_ndc, xjust = 0.5, yjust = 1, horiz = TRUE, bty = "n", cex = 0.8,
           legend = paste0(c(sprintf("JAGWAS (%d df)", k), "JAGWAS, generalized inverse (df = rank)",
                             "Bonferroni minP"), "    "),   # trailing spaces keep the keys apart
           col = c(col_j, col_j, col_m), lty = c(1, 2, 1), pch = c(16, 17, 15), lwd = 2)
    dev.off()
    msg("Wrote ", f)
  }
}
msg("Done.")
