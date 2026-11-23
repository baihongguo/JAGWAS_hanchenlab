# power_common.R -- helpers shared by the pw*.R power-simulation scripts.
# Source after common.R (which supplies argument parsing, logging, data.table).

POWER_ROOT <- "/data5/bguo1/JAGWAS_power"

load_design <- function(root) {
  f <- file.path(root, "design", "scenarios.rds")
  if (!file.exists(f)) stop("Missing ", f, ": run pw01_define_scenarios.R first.", call. = FALSE)
  readRDS(f)
}

# Correlation matrix of the phenotypes Y = g beta' + E of one data set, from the
# data set's cross-products (everything centred, g standardized):
#     Y'Y = E'E + beta c' + c beta' + (g'g) beta beta',   c = E'g
# This is exactly cor(Y) -- pw02 --check=true verifies it -- without forming Y.
rhat_update <- function(SE, cvec, gg, beta) {
  cov2cor(SE + tcrossprod(beta, cvec) + tcrossprod(cvec, beta) + gg * tcrossprod(beta))
}

# Scenario-independent pieces of every simulated data set, written by pw02.
load_datasets <- function(root, design) {
  fs <- sort(list.files(file.path(root, "data"), pattern = "^ds_[0-9]+-[0-9]+\\.rds$",
                        full.names = TRUE))
  if (!length(fs)) stop("No simulated data sets in ", file.path(root, "data"),
                        ": run pw02_simulate_sumstats.R first.", call. = FALSE)
  parts <- lapply(fs, readRDS)
  ds <- unlist(lapply(parts, `[[`, "ds"))
  if (anyDuplicated(ds)) stop("Overlapping data-set ranges in ", file.path(root, "data"))
  k  <- design$k
  SE <- array(0, c(k, k, length(ds)))
  at <- 0L
  for (p in parts) { m <- length(p$ds); SE[, , at + seq_len(m)] <- p$SE; at <- at + m }
  pos <- integer(max(ds)); pos[ds] <- seq_along(ds)
  out <- list(ds = ds, pos = pos, SE = SE,
              cvec     = do.call(rbind, lapply(parts, `[[`, "cvec")),
              gg       = unlist(lapply(parts, `[[`, "gg")),
              z_null   = do.call(rbind, lapply(parts, `[[`, "z_null")),
              null_ds  = unlist(lapply(parts, `[[`, "null_ds")),
              null_snp = unlist(lapply(parts, `[[`, "null_snp")))
  if (length(ds) < design$n_ds)
    msg(sprintf("  ** only %d of %d data sets simulated so far", length(ds), design$n_ds))
  out
}

# Causal-variant summary statistics of one scenario, all ranges, ordered by data set.
load_sumstats <- function(root, scen) {
  fs <- sort(list.files(file.path(root, "sumstats", scen), pattern = "^ss_[0-9]+-[0-9]+\\.rds$",
                        full.names = TRUE))
  if (!length(fs)) stop("No summary statistics for ", scen, call. = FALSE)
  parts <- lapply(fs, readRDS)
  ds <- unlist(lapply(parts, `[[`, "ds"))
  o  <- order(ds)
  bind <- function(x) do.call(rbind, lapply(parts, `[[`, x))[o, , drop = FALSE]
  list(ds = ds[o], t = bind("t"), b = bind("b"), se = bind("se"),
       af = unlist(lapply(parts, `[[`, "af"))[o])
}

# What a method tests in one scenario: the causal variant of every data set, or
# (scen = "null") every null variant. Each row comes with the phenotypic
# correlation matrix of its own data set, the matrix every method is given.
test_inputs <- function(root, design, D, scen) {
  if (scen == "null") {
    ok <- is.finite(rowSums(D$z_null))
    list(z = D$z_null[ok, , drop = FALSE], ds = D$null_ds[ok], snp = D$null_snp[ok],
         rhat = function(j) cov2cor(D$SE[, , D$pos[j]]))
  } else {
    S <- load_sumstats(root, scen)
    B <- design$beta[[scen]]
    list(z = S$t, ds = S$ds, snp = rep(0L, length(S$ds)), af = S$af,
         rhat = function(j) {
           i <- D$pos[j]
           rhat_update(D$SE[, , i], D$cvec[i, ], D$gg[i], B[j, ])
         })
  }
}

# Applies test_fun(z_rows, Rhat) -> list(stat, p), one data set at a time.
by_dataset <- function(inp, test_fun) {
  grp <- split(seq_along(inp$ds), inp$ds)
  rbindlist(lapply(grp, function(r) {
    j   <- inp$ds[r[1]]
    res <- test_fun(inp$z[r, , drop = FALSE], inp$rhat(j))
    data.table(ds = j, snp = inp$snp[r], stat = res$stat, p = res$p)
  }))
}

scenario_names <- function(which, design) {
  all <- c("null", design$scen$name)
  if (identical(which, "all")) return(all)
  s <- trimws(strsplit(which, ",", fixed = TRUE)[[1]])
  bad <- setdiff(s, all)
  if (length(bad)) stop("Unknown scenario(s): ", paste(bad, collapse = ", "),
                        "\n  Valid: ", paste(all, collapse = ", "), call. = FALSE)
  s
}

pval_file <- function(root, method, scen) file.path(root, "pvalues", method, paste0(scen, ".txt"))

save_pvals <- function(dt, root, method, scen, secs, mem_mb) {
  f <- pval_file(root, method, scen)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, f, sep = "\t")
  fwrite(data.table(method = method, scenario = scen, n_tests = nrow(dt),
                    seconds = round(secs, 2), peak_mem_mb = round(mem_mb, 1)),
         sub("\\.txt$", ".runtime.txt", f), sep = "\t")
}

# Evaluates expr, returning its value with elapsed seconds and the peak memory
# R used meanwhile (MB).
timed <- function(expr) {
  invisible(gc(reset = TRUE))
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  secs <- proc.time()[["elapsed"]] - t0
  g <- gc()
  list(value = value, secs = secs, mem_mb = sum(g[, which(colnames(g) == "max used") + 1]))
}

# One-line power / calibration report for a finished scenario.
report <- function(method, scen, dt, alpha, secs) {
  if (scen == "null") {
    msg(sprintf("%-10s %-22s %7d null tests | P<1e-3: %.3f x expected | P<1e-4: %.3f x | %.1f s",
                method, scen, nrow(dt), mean(dt$p < 1e-3, na.rm = TRUE) / 1e-3,
                mean(dt$p < 1e-4, na.rm = TRUE) / 1e-4, secs))
  } else {
    msg(sprintf("%-10s %-22s power at %g: %.3f  (%d variants, %d NA) | %.1f s",
                method, scen, alpha, mean(dt$p < alpha, na.rm = TRUE), nrow(dt),
                sum(is.na(dt$p)), secs))
  }
}

# Standard boilerplate for a method script: parses the common options and loads
# the design and the simulated data sets.
method_setup <- function(method, extra = list()) {
  defaults <- c(list(root = POWER_ROOT, scenarios = "all", overwrite = "false",
                     threads = "1"), extra)
  opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
  thr <- as_int(opt$threads, "threads"); if (thr > 0) setDTthreads(thr)
  log_open(file.path(opt$root, "logs", sprintf("pw03_%s.log", method)))
  design <- load_design(opt$root)
  msg(sprintf("%s: %d-trait design, %d data sets, alpha %g", method, design$k, design$n_ds,
              design$alpha))
  D <- load_datasets(opt$root, design)
  list(opt = opt, design = design, D = D, overwrite = as_bool(opt$overwrite, "overwrite"),
       scen = scenario_names(opt$scenarios, design))
}
