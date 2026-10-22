# =============================================================================
# steps/peak_region_size.R — peak width distribution (encode_lib_genomic.get_region_size_metrics)
#   Rscript peak_region_size.R <peak.narrowPeak.gz> <out_prefix>
#   writes <out_prefix>.peak_region_size.qc  (Min, quartiles, Max, Mean; pandas describe() = R quantile type 7)
#          <out_prefix>.peak_region_size.pdf (100-bin histogram line plot, numpy.histogram bins)
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
peak <- args[1]; prefix <- args[2]

con <- gzfile(peak)
lines <- readLines(con); close(con)
size <- if (length(lines) == 0) numeric(0) else
  vapply(strsplit(lines, "\t", fixed = TRUE), function(x) as.numeric(x[3]) - as.numeric(x[2]), numeric(1))

if (length(size) == 0) {
  stats <- rep(NA_real_, 6)
} else {
  q <- quantile(size, c(0, 0.25, 0.5, 0.75, 1), type = 7, names = FALSE)
  stats <- c(q, mean(size))
}
keys <- c("Min size", "25 percentile", "50 percentile (median)", "75 percentile", "Max size", "Mean")
writeLines(paste(keys, ifelse(is.na(stats), "NA", formatC(stats, digits = 15, format = "g")), sep = "\t"), paste0(prefix, ".peak_region_size.qc"))

pdf(paste0(prefix, ".peak_region_size.pdf"), width = 6.4, height = 4.8)
if (length(size) > 0) {
  lo <- min(size); hi <- max(size)
  if (lo == hi) { lo <- lo - 0.5; hi <- hi + 0.5 }                  # numpy widens a zero-width range
  h <- hist(size, breaks = seq(lo, hi, length.out = 101), right = FALSE, include.lowest = TRUE, plot = FALSE)
  plot(h$mids, h$counts, type = "l", col = "#1f77b4", xlab = "peak width (bp)", ylab = "count",
       main = paste("Peak width distribution for", basename(peak)), cex.main = 0.8)
} else {
  plot.new(); title(main = paste("No peaks in", basename(peak)), cex.main = 0.8)
}
invisible(dev.off())
