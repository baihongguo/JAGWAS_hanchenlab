#!/usr/bin/env Rscript
# pw03_metausat.R -- metaUSAT on the power-simulation summary statistics,
# using the authors' own R code (Ray & Boehnke; github.com/RayDebashree/metaUSAT).
# metaUSAT adaptively combines metaMANOVA (z' R^-1 z, i.e. JAGWAS's statistic)
# with metaSSU (z'z) and reports one p-value per variant.
#
# Each variant is passed with the phenotypic correlation of its own data set,
# the same matrix the other methods get.
#
# metaUSAT evaluates its p-value by numerical integration, one variant at a
# time, so it is by far the slowest method here. Null calibration therefore
# uses a random subset of --null_n null variants rather than all of them.
#
# RUN THE SMOKE TEST FIRST. The function is called as metausat(z, R) with z the
# k z-scores of one variant; --smoke=true runs one variant and prints the
# returned object, so you can confirm the call and which element holds the
# p-value before committing to the full run:
#
#   Rscript pw03_metausat.R --src=/path/to/metaUSAT_v1.17.R --smoke=true
#   Rscript pw03_metausat.R --src=/path/to/metaUSAT_v1.17.R
#
# Writes <root>/pvalues/metausat/<scenario>.txt.

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

M <- "metausat"
S <- method_setup(M, list(src = "", fun = "metausat", null_n = "20000", smoke = "false",
                          seed = "11"))
if (!nzchar(S$opt$src) || !file.exists(S$opt$src))
  stop("Give the metaUSAT R source with --src=/path/to/metaUSAT_vX.R", call. = FALSE)
env <- new.env()
sys.source(S$opt$src, envir = env)
if (!exists(S$opt$fun, envir = env, inherits = FALSE))
  stop("No function '", S$opt$fun, "' in ", S$opt$src, ". Functions found: ",
       paste(ls(env)[vapply(ls(env), function(x) is.function(env[[x]]), TRUE)], collapse = ", "),
       call. = FALSE)
fn <- get(S$opt$fun, envir = env)
msg("metaUSAT function: ", S$opt$fun, "(", paste(names(formals(fn)), collapse = ", "), ")")

# The p-value element of metaUSAT's result, by name.
p_of <- function(res) {
  if (is.numeric(res) && length(res) == 1) return(res)
  nm <- names(res)
  hit <- intersect(c("p.metausat", "pval.metausat", "p.value", "pvalue", "p"), nm)
  if (!length(hit)) hit <- grep("^p.*usat", nm, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) stop("Cannot find the metaUSAT p-value among: ", paste(nm, collapse = ", "),
                         call. = FALSE)
  as.numeric(res[[hit[1]]])[1]
}

n_fail <- 0L
metausat_test <- function(z, R) {
  p <- vapply(seq_len(nrow(z)), function(i) {
    res <- tryCatch(fn(z[i, ], R), error = function(e) e)
    if (inherits(res, "error")) { n_fail <<- n_fail + 1L; return(NA_real_) }
    p_of(res)
  }, 0)
  list(stat = NA_real_, p = p)
}

if (as_bool(S$opt$smoke, "smoke")) {
  scen <- S$design$scen$name[1]
  inp  <- test_inputs(S$opt$root, S$design, S$D, scen)
  res  <- fn(inp$z[1, ], inp$rhat(inp$ds[1]))
  cat("metaUSAT returned:\n"); str(res)
  cat(sprintf("\nExtracted p-value: %g\n", p_of(res)))
  quit(save = "no")
}

for (scen in S$scen) {
  if (!S$overwrite && file.exists(pval_file(S$opt$root, M, scen))) next
  inp <- test_inputs(S$opt$root, S$design, S$D, scen)
  if (scen == "null") {
    set.seed(as_int(S$opt$seed, "seed"))
    keep <- sort(sample.int(length(inp$ds), min(length(inp$ds), as_int(S$opt$null_n, "null_n"))))
    inp$z <- inp$z[keep, , drop = FALSE]; inp$ds <- inp$ds[keep]; inp$snp <- inp$snp[keep]
  }
  n_fail <- 0L
  r <- timed(by_dataset(inp, metausat_test))
  save_pvals(r$value, S$opt$root, M, scen, r$secs, r$mem_mb)
  report(M, scen, r$value, S$design$alpha, r$secs)
  if (n_fail) msg(sprintf("  ** %d metaUSAT calls failed (p set to NA)", n_fail))
}
msg("Done.")
