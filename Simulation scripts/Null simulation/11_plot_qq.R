#!/usr/bin/env Rscript
# 11_plot_qq.R
#
# QQ plots of the pooled JAGWAS null p-values, one panel per scenario, with a
# grey 95% band.
#
#   --variants=all     every variant   (05's direct_hist.txt.gz)
#   --variants=common  MAF > 1%        (10's direct_hist_common.txt.gz)
#   --variants=rare    MAF <= 1%       (10's direct_hist_rare.txt.gz)
#
# Writes <inbase>/results/qq_<variants>.png (a PDF if PNG is unavailable) and
# qq_points_<variants>.txt, the points as plotted.
#
# THE GREY BAND. For N independent uniform p-values, the i-th smallest follows
# Beta(i, N - i + 1), which here is Gamma(i, 1)/N to any visible precision. The
# band is its 2.5% and 97.5% quantiles on the -log10 scale. It assumes the
# p-values are independent; SNPs in LD are not, so near the top of the plot the
# real spread is wider than the band. The replicate-bootstrap intervals in
# type1_table_*.txt are the LD-aware inference; the band is the conventional
# visual reference.
#
# THE POINTS. The p-values are held as histograms of -log10 p in steps of 0.01,
# pooled over replicates. Each non-empty step is drawn once, at its mid-point,
# against the expected quantile of its middle rank. The histograms stop at 100,
# so p-values below 1e-100 are drawn as triangles along y = 100.
#
# OTHER RESULT FOLDERS. --scenarios lists the folders to plot, in panel order,
# e.g. S4_adjfill (05 run with 04's overlap-adjusted matrix) in place of
# S4_normal_miss1pct. Such panels keep the scenario title and say "overlap-
# adjusted R" next to lambda_GC. A non-default list adds --tag to the output
# names ("custom" if no tag is given), so the original plots are kept.
#
# TRAIT SUBSETS. Folders from 12_trait_subsets.R end in _k<number>, used as the
# degrees of freedom for lambda_GC (--k for the others). Their panels take the
# scenario's title; when the figure is not 128 phenotypes throughout, every panel
# title also gives its k.
#
# Usage:
#   Rscript 11_plot_qq.R --variants=all
#   Rscript 11_plot_qq.R --variants=common
#   Rscript 11_plot_qq.R --variants=common \
#       --scenarios=S1_normal,S2_lognormal,S3_cauchy,S4_adjfill --tag=S4adj
#   Rscript 11_plot_qq.R --variants=common --tag=k16 \
#       --scenarios=S1_normal_k16,S2_lognormal_k16,S3_cauchy_k16,S4_adjfill_k16
#   --ncol=4 with 12 folders (16, 64, 128 phenotypes) gives one row per number of phenotypes.

suppressPackageStartupMessages(source(file.path(
  (function() {
    a <- commandArgs(trailingOnly = FALSE)
    f <- grep("^--file=", a, value = TRUE)
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  })(), "common.R")))

defaults <- list(
  inbase    = "/data5/bguo1/JAGWAS_null_input",
  variants  = "all",      # all | common | rare
  k         = "128",
  scenarios = "S1_normal,S2_lognormal,S3_cauchy,S4_normal_miss1pct",   # panel order
  tag       = "",         # appended to the output names; "custom" when --scenarios
                          # is changed and no tag is given
  ncol      = "2"         # panels per row, filled row by row
)
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
k   <- as_int(opt$k, "k")
ncol_p <- as_int(opt$ncol, "ncol")
if (!opt$variants %in% c("all", "common", "rare"))
  stop("--variants must be all, common or rare.")
f_hist <- if (opt$variants == "all") "direct_hist.txt.gz" else
  sprintf("direct_hist_%s.txt.gz", opt$variants)
scen <- trimws(strsplit(opt$scenarios, ",", fixed = TRUE)[[1]])
tag  <- opt$tag
if (!nzchar(tag) && !identical(opt$scenarios, defaults$scenarios)) tag <- "custom"
sfx_out <- paste0(opt$variants, if (nzchar(tag)) paste0("_", tag) else "")

TITLES <- c(S1_normal          = "Normal",
            S2_lognormal       = "Log-normal",
            S3_cauchy          = "Cauchy",
            S4_normal_miss1pct = "Normal with missing",
            S4_adjfill         = "Normal with missing",
            S4_adjdrop         = "Normal with missing")
# Shown before lambda_GC, to mark panels tested with 04's overlap-adjusted matrix.
NOTES  <- c(S4_adjfill = "overlap-adjusted R",
            S4_adjdrop = "overlap-adjusted R (drop)")
# minP folders from 13_minp_pass.R (<scenario>_minpB_k<k>, _minpE_k<k>): the method
# replaces the note, and there is no lambda_GC.
MINP   <- c(minpB = "Bonferroni minP", minpE = "Effective-tests minP")
OUTER <- switch(opt$variants,
                all    = NULL,
                common = "Common variants (MAF > 1%)",
                rare   = quote("Rare variants (MAF" <= "1%)"))

# Histogram layout -- must match 05 and 10.
TNB <- 40001L; TBW <- 0.01   # T on [0, 400)
LNB <- 10001L; LBW <- 0.01   # -log10 p on [0, 100); bin LNB holds everything beyond
TOP <- (LNB - 1) * LBW       # 100

resdir <- file.path(opt$inbase, "results")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
log_open(file.path(opt$inbase, "logs", sprintf("11_plot_qq_%s.log", sfx_out)))
on.exit(log_close(), add = TRUE)

read_gz <- function(f) fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE)

lambda_T <- function(th, df) {
  if (sum(th) == 0) return(NA_real_)                 # minP folders have no T histogram
  b <- which(cumsum(th) >= sum(th) / 2)[1]
  if (is.na(b) || b == TNB) NA_real_ else (b - 0.5) * TBW / qchisq(0.5, df = df)
}

# 95% band for the order statistics of N independent uniforms, on a log-spaced
# grid of ranks. Gamma quantiles for small ranks; the normal limit beyond 1e6,
# where the two agree to far below plotting resolution.
null_band <- function(N) {
  i <- unique(c(1, round(exp(seq(0, log(N), length.out = 800)))))
  lo <- hi <- numeric(length(i))
  s <- i < 1e6
  lo[s]  <- qgamma(0.025, shape = i[s]);  hi[s]  <- qgamma(0.975, shape = i[s])
  lo[!s] <- i[!s] - 1.959964 * sqrt(i[!s]); hi[!s] <- i[!s] + 1.959964 * sqrt(i[!s])
  list(x  = -log10((i - 0.5) / N),
       up = -log10(lo / N),
       dn = -log10(pmin(hi / N, 1)))
}

panels <- list()
for (s in scen) {
  sdir <- file.path(opt$inbase, s)
  reps <- sort(list.files(sdir, pattern = "^rep[0-9]{3}$"))
  reps <- reps[file.exists(file.path(sdir, reps, f_hist))]
  if (!length(reps)) { msg(s, ": no ", f_hist, " files, skipping."); next }

  hs <- rbindlist(lapply(reps, function(rp) read_gz(file.path(sdir, rp, f_hist))))
  hs <- hs[, .(count = sum(count)), by = .(type, bin)]
  lh <- numeric(LNB); lh[hs[type == "lp", bin]] <- hs[type == "lp", count]
  th <- numeric(TNB); th[hs[type == "T",  bin]] <- hs[type == "T",  count]
  N  <- sum(lh)

  # Non-empty bins, most significant first, with the ranks each one covers.
  bins  <- rev(which(lh > 0))
  last  <- cumsum(lh[bins])
  first <- last - lh[bins] + 1
  pts <- data.table(expected = -log10(((first + last) / 2 - 0.5) / N),
                    observed = (bins - 0.5) * LBW, beyond = FALSE)
  if (lh[LNB] > 0) {
    # Everything past the histogram's top: spread over its rank range at y = 100.
    rk  <- unique(round(exp(seq(0, log(lh[LNB]), length.out = 40))))
    pts <- rbind(data.table(expected = -log10((rk - 0.5) / N), observed = TOP, beyond = TRUE),
                 pts[bins != LNB])
  }

  ks <- scenario_k(s, k); sb <- scenario_base(s)
  mth <- regmatches(sb, regexpr("minp[BE]$", sb))
  if (length(mth)) sb <- sub("_minp[BE]$", "", sb)
  panels[[s]] <- list(pts = pts, N = N, lgc = lambda_T(th, ks), R = length(reps), k = ks,
                      title = if (sb %in% names(TITLES)) TITLES[[sb]] else s,
                      note  = if (length(mth)) MINP[[mth]] else
                              if (sb %in% names(NOTES)) NOTES[[sb]] else NA_character_)
  msg(sprintf("%-20s k = %3d  %3d replicates  %.3g p-values  lambda_GC %.4f%s",
              s, ks, length(reps), N, panels[[s]]$lgc,
              if (lh[LNB] > 0) sprintf("  (%.0f beyond 1e-100)", lh[LNB]) else ""))
}
if (!length(panels)) stop("Nothing to plot.")
if (any(vapply(panels, `[[`, 0L, "k") != k))
  for (s in names(panels)) panels[[s]]$title <- sprintf("%s, %d phenotypes", panels[[s]]$title, panels[[s]]$k)

np <- length(panels)
nc <- min(np, ncol_p); nr <- ceiling(np / nc)

draw <- function() {
  par(mfrow = c(nr, nc),
      mar = c(4.4, 4.8, 3.4, 1.2), oma = c(0, 0, if (is.null(OUTER)) 0 else 2.2, 0),
      mgp = c(2.9, 0.75, 0), las = 1, pty = "s")
  for (s in names(panels)) {
    p    <- panels[[s]]
    xmax <- -log10(0.5 / p$N)                       # expected value of the top rank
    ymax <- max(p$pts$observed, xmax)
    xlim <- c(0, xmax * 1.03); ylim <- c(0, ymax * 1.03)
    b    <- null_band(p$N)

    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
         xlab = expression(Expected ~ -log[10](italic(p))),
         ylab = expression(Observed ~ -log[10](italic(p))))
    polygon(c(b$x, rev(b$x)), c(b$dn, rev(b$up)), col = "grey82", border = NA)
    abline(0, 1, col = "grey30", lwd = 1)
    q <- p$pts[beyond == FALSE]
    points(q$expected, q$observed, pch = 20, cex = 0.45)
    q <- p$pts[beyond == TRUE]
    if (nrow(q)) points(q$expected, q$observed, pch = 17, cex = 0.6, xpd = NA)
    box()

    title(main = p$title, line = 1.8, cex.main = 1.25)
    lg <- sprintf("%.4f", p$lgc)
    sub_line <- if (is.na(p$lgc)) p$note else if (is.na(p$note)) bquote(lambda[GC] == .(lg)) else
                  bquote(.(p$note) * ",  " * lambda[GC] == .(lg))
    if (!is.null(sub_line) && !identical(sub_line, NA_character_))
      mtext(sub_line, side = 3, line = 0.4, cex = 0.9)
  }
  if (!is.null(OUTER)) mtext(OUTER, outer = TRUE, line = 0.6, cex = 1.15)
}

wd <- 1000 * nc                                     # 1000 px per panel, as before
h  <- 1000 * nr + if (is.null(OUTER)) 0 else 90
f_png <- file.path(resdir, sprintf("qq_%s.png", sfx_out))
opened <- tryCatch({
  png(f_png, width = wd, height = h, res = 250,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  TRUE
}, error = function(e) FALSE)
if (opened) {
  draw(); invisible(dev.off())
} else {
  f_png <- sub("\\.png$", ".pdf", f_png)
  pdf(f_png, width = 4 * nc, height = 4 * h / 1000); draw(); invisible(dev.off())
}

pts_out <- rbindlist(lapply(names(panels), function(s) panels[[s]]$pts[, scenario := s]))
setcolorder(pts_out, "scenario")
f_pts <- file.path(resdir, sprintf("qq_points_%s.txt", sfx_out))
fwrite(pts_out, f_pts, sep = "\t")

msg("Wrote ", f_png)
msg("Wrote ", f_pts)
