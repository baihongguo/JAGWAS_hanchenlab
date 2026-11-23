#!/usr/bin/env Rscript
# pz02_plot_power.R
#
# Redraws the power figure of pz01_zscore_power.R from its zpower_results.tsv,
# without rerunning the simulation, e.g. after a change of labels. The drawing
# code is pz01's; the two effect patterns are labelled as in the manuscript:
# equal effects on the associated phenotypes ("direct" in the results file) and
# weighted effects, delta = R w ("correlated").
#
# Usage:
#   Rscript pz02_plot_power.R --results=zscore_T1_v2/zpower_results.tsv --out=zscore_T1_v2
#   Rscript pz02_plot_power.R --results=... --out=... --format=pdf --targets=0.8

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(results = "zpower_results.tsv", out = ".", format = "png", res = "300",
            targets = "all", alpha = "5e-8")
for (a in args) {
  if (!grepl("^--[^=]+=", a)) stop("Malformed argument '", a, "'. Use --key=value.")
  k_ <- sub("^--([^=]+)=.*$", "\\1", a)
  if (!k_ %in% names(opt)) stop("Unknown option --", k_, ". Valid options: ", paste(names(opt), collapse = ", "))
  opt[[k_]] <- sub("^--[^=]+=", "", a)
}
if (!opt$format %in% c("png", "pdf")) stop("--format must be png or pdf.")

res <- fread(opt$results)
k <- max(res$df)
alpha <- as.numeric(opt$alpha)
scens <- intersect(c("empirical", "nearsingular", "signflip"), unique(res$scenario))
pats <- intersect(c("direct", "correlated"), unique(res$pattern))
m_lev <- sort(unique(res$m))
tl <- unique(res[, .(target, lambda)])[order(-target)]
if (opt$targets != "all") tl <- tl[target %in% as.numeric(strsplit(opt$targets, ",")[[1]])]

lab_s <- c(empirical = "Empirical UDIP correlation",
           nearsingular = "Nearly singular (64 UDIPs + near-duplicates)",
           signflip = "Mixed effect directions")
lab_p <- c(direct = "Equal effects", correlated = "Weighted effects")
col_j <- "black"; col_m <- "#D55E00"
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

for (ti in seq_len(nrow(tl))) {
  f <- file.path(opt$out, sprintf("zpower_power_%02d.%s", round(100 * tl$target[ti]), opt$format))
  h_in <- (1050 * length(pats) + 360) / 300
  if (opt$format == "png") png(f, width = 10, height = h_in, units = "in", res = as.numeric(opt$res)) else
    pdf(f, width = 10, height = h_in)
  par(mfrow = c(length(pats), length(scens)), mar = c(4, 4, 2.6, 0.6), oma = c(0, 0, 4.4, 0),
      mgp = c(2.4, 0.7, 0), las = 1, cex = 0.8)
  for (p in pats) for (s in scens) {
    x <- res[scenario == s & pattern == p & target == tl$target[ti]][order(m)]
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
                tl$lambda[ti], 100 * tl$target[ti], k, alpha), outer = TRUE, line = 2.9, cex = 0.85)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  legend(0.5, top_in + 2.1 * line_ndc, xjust = 0.5, yjust = 1, horiz = TRUE, bty = "n", cex = 0.8,
         legend = paste0(c(sprintf("JAGWAS (%d df)", k), "JAGWAS, generalized inverse (df = rank)",
                           "Bonferroni minP"), "    "),   # trailing spaces keep the keys apart
         col = c(col_j, col_j, col_m), lty = c(1, 2, 1), pch = c(16, 17, 15), lwd = 2)
  dev.off()
  cat("Wrote", f, "\n")
}
