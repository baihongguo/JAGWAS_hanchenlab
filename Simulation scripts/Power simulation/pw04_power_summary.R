#!/usr/bin/env Rscript
# pw04_power_summary.R -- collects every method's p-values into the results:
#
#   <root>/results/power_table.txt        power at alpha per scenario x method,
#                                         with exact (Clopper-Pearson) 95% CIs
#   <root>/results/power_by_maf.txt       the same by causal MAF bin
#   <root>/results/calibration_table.txt  null variants: observed / expected at
#                                         1e-2 ... 1e-5
#   <root>/results/runtime_table.txt      seconds and peak R memory per method
#   <root>/results/power.png              power by architecture (PDF if no PNG)
#
# Under --signal=mahalanobis JAGWAS should sit at the design's target power in
# every scenario (dashed line): that is the check that the simulation, the
# summary statistics and the correlation matrices are all consistent.
#
# Usage:
#   Rscript pw04_power_summary.R

.here <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})()
suppressPackageStartupMessages({
  source(file.path(.here, "common.R"))
  source(file.path(.here, "power_common.R"))
})

opt    <- parse_args(commandArgs(trailingOnly = TRUE), list(root = POWER_ROOT))
design <- load_design(opt$root)
alpha  <- design$alpha
resdir <- file.path(opt$root, "results")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$root, "logs", "pw04_power_summary.log"))

# Fixed order: each method keeps its colour, marker and line in every figure.
METHODS <- c(jagwas    = "JAGWAS",           multiabel = "MultiABEL",
             mostest   = "MOSTest",          metausat  = "metaUSAT",
             minp_bonf = "minP (Bonferroni)", minp_eff = "minP (effective tests)")
COL <- c(jagwas = "#2a78d6", multiabel = "#eb6834", mostest = "#1baf7a",
         metausat = "#eda100", minp_bonf = "#e87ba4", minp_eff = "#008300")
PCH <- c(jagwas = 19, multiabel = 2, mostest = 15, metausat = 17, minp_bonf = 18, minp_eff = 4)
LTY <- c(jagwas = 1, multiabel = 2, mostest = 1, metausat = 1, minp_bonf = 1, minp_eff = 3)
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"

methods <- names(METHODS)[dir.exists(file.path(opt$root, "pvalues", names(METHODS)))]
if (!length(methods)) stop("No p-values under ", file.path(opt$root, "pvalues"))
msg("Methods found: ", paste(METHODS[methods], collapse = ", "))

maf_bin <- cut(design$maf_causal, c(0, 0.05, 0.5), right = FALSE, include.lowest = TRUE,
               labels = c("MAF 0.01-0.05", "MAF 0.05-0.5"))

ci <- function(x, n) if (n > 0) binom.test(x, n)$conf.int else c(NA_real_, NA_real_)

# ------------------------------------------------------------------ power --

pw <- list(); pmaf <- list()
for (mth in methods) for (i in seq_len(nrow(design$scen))) {
  s <- design$scen[i]
  f <- pval_file(opt$root, mth, s$name)
  if (!file.exists(f)) next
  d <- fread(f)[snp == 0]
  hit <- d$p < alpha
  x <- sum(hit, na.rm = TRUE); n <- sum(!is.na(hit)); b <- ci(x, n)
  pw[[length(pw) + 1L]] <- data.table(scenario = s$name, arch = s$arch, m = s$m,
    direction = s$direction, method = mth, power = x / n, ci_lo = b[1], ci_hi = b[2],
    n = n, n_na = sum(is.na(hit)))
  bins <- maf_bin[d$ds]
  for (lv in levels(maf_bin)) {
    h <- hit[bins == lv]; x2 <- sum(h, na.rm = TRUE); n2 <- sum(!is.na(h)); b2 <- ci(x2, n2)
    pmaf[[length(pmaf) + 1L]] <- data.table(scenario = s$name, method = mth, maf_bin = lv,
      power = if (n2) x2 / n2 else NA_real_, ci_lo = b2[1], ci_hi = b2[2], n = n2)
  }
}
pw <- rbindlist(pw); pmaf <- rbindlist(pmaf)
fwrite(pw, file.path(resdir, "power_table.txt"), sep = "\t")
fwrite(pmaf, file.path(resdir, "power_by_maf.txt"), sep = "\t")

wide <- dcast(pw, scenario ~ method, value.var = "power")
setcolorder(wide, c("scenario", intersect(names(METHODS), names(wide))))
msg("")
msg("Power at ", alpha, ":")
msg(paste(sprintf("%-18s", c("scenario", METHODS[setdiff(names(wide), "scenario")])), collapse = ""))
for (i in seq_len(nrow(wide)))
  msg(paste(c(sprintf("%-18s", wide$scenario[i]),
              sprintf("%-18.3f", unlist(wide[i, -1, with = FALSE]))), collapse = ""))

target <- if (design$signal == "mahalanobis") mean(design$scen$jagwas_power_analytic) else NA_real_
if (!is.na(target) && "jagwas" %in% methods) {
  j <- pw[method == "jagwas"]
  msg(sprintf("JAGWAS check: design target %.3f; simulated power ranges %.3f to %.3f across scenarios",
              target, min(j$power), max(j$power)))
}

# ------------------------------------------------------- calibration, time --

cal <- list()
for (mth in methods) {
  f <- pval_file(opt$root, mth, "null")
  if (!file.exists(f)) next
  p <- fread(f)$p; p <- p[!is.na(p)]
  for (thr in c(1e-2, 1e-3, 1e-4, 1e-5))
    cal[[length(cal) + 1L]] <- data.table(method = mth, threshold = thr, n_null = length(p),
      observed = sum(p < thr), expected = length(p) * thr, ratio = sum(p < thr) / (length(p) * thr))
}
cal <- rbindlist(cal)
if (nrow(cal)) fwrite(cal, file.path(resdir, "calibration_table.txt"), sep = "\t")

rt <- rbindlist(lapply(list.files(file.path(opt$root, "pvalues"), pattern = "\\.runtime\\.txt$",
                                  recursive = TRUE, full.names = TRUE), fread))
if (nrow(rt)) fwrite(rt, file.path(resdir, "runtime_table.txt"), sep = "\t")

# ------------------------------------------------------------------ figure --

draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.8, 1), oma = c(5, 0, 2, 0),
      mgp = c(2.6, 0.7, 0), las = 1, col.axis = INK2, col.lab = INK, fg = AXIS, bg = "white")
  frame_panel <- function(xlim, main, xlab) {
    plot(NA, xlim = xlim, ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = xlab, ylab = "Power")
    abline(h = seq(0, 1, 0.2), col = GRID, lwd = 0.8)
    axis(2, at = seq(0, 1, 0.2), col = AXIS, col.ticks = AXIS)
    if (!is.na(target)) abline(h = target, lty = 2, col = MUTED, lwd = 1)
    title(main = main, col.main = INK, cex.main = 1.1, line = 1)
  }
  single <- pw[direction == "single"]
  titles <- c(concordant = "Concordant effects", discordant = "Discordant effects",
              mixed = "Mixed-sign effects")
  for (d in intersect(names(titles), pw[arch == "random", unique(direction)])) {
    sub <- rbind(single, pw[arch == "random" & direction == d])
    xs  <- sort(unique(log2(sub$m)))
    frame_panel(range(xs) + c(-0.3, 0.3), titles[[d]], "Number of affected traits")
    axis(1, at = xs, labels = 2^xs, col = AXIS, col.ticks = AXIS)
    for (mth in methods) {
      q <- sub[method == mth][order(m)]
      lines(log2(q$m), q$power, col = COL[[mth]], lwd = 1.6, lty = LTY[[mth]])
      points(log2(q$m), q$power, col = COL[[mth]], pch = PCH[[mth]], cex = 1)
    }
    box(col = AXIS)
  }
  b <- pw[arch == "block"]
  if (nrow(b)) {
    dl  <- intersect(c("concordant", "discordant", "mixed"), unique(b$direction))
    off <- seq(-0.3, 0.3, length.out = length(methods))
    frame_panel(c(0.5, length(dl) + 0.5),
                sprintf("One correlated block (%d traits)", length(design$block)), "Effect direction")
    axis(1, at = seq_along(dl), labels = dl, col = AXIS, col.ticks = AXIS)
    for (i in seq_along(methods)) {
      q <- b[method == methods[i]]
      x <- match(q$direction, dl) + off[i]
      segments(x, q$ci_lo, x, q$ci_hi, col = COL[[methods[i]]], lwd = 1.2)
      points(x, q$power, col = COL[[methods[i]]], pch = PCH[[methods[i]]], cex = 1)
    }
    box(col = AXIS)
  }
  mtext(sprintf("Power at alpha = %g, %d causal variants per scenario; %s noncentrality fixed at %.1f",
                alpha, design$n_ds, design$signal, design$ncp),
        side = 3, outer = TRUE, line = 0.4, col = INK, cex = 0.9)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = c(METHODS[methods], if (!is.na(target)) "JAGWAS expected"),
         col = c(COL[methods], if (!is.na(target)) MUTED),
         pch = c(PCH[methods], if (!is.na(target)) NA),
         lty = c(LTY[methods], if (!is.na(target)) 2), lwd = 1.6,
         ncol = 4, bty = "n", text.col = INK, cex = 0.9, inset = 0.01)
}

f_png <- file.path(resdir, "power.png")
opened <- tryCatch({
  png(f_png, width = 2000, height = 1950, res = 250,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  TRUE
}, error = function(e) FALSE)
if (opened) { draw(); invisible(dev.off()) } else {
  f_png <- sub("\\.png$", ".pdf", f_png)
  pdf(f_png, width = 8, height = 7.8); draw(); invisible(dev.off())
}

msg("")
msg("Wrote ", file.path(resdir, "power_table.txt"))
msg("Wrote ", file.path(resdir, "power_by_maf.txt"))
if (nrow(cal)) msg("Wrote ", file.path(resdir, "calibration_table.txt"))
if (nrow(rt))  msg("Wrote ", file.path(resdir, "runtime_table.txt"))
msg("Wrote ", f_png)
