# Usage: Rscript count_reads.R samplesheet.csv peak_counts.tsv [threads]
library(DiffBind)

args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 3) as.integer(args[3]) else 1

dbo <- dba(sampleSheet = args[1])
dbo$config$cores <- cores
dbo <- dba.count(dbo, summits = FALSE, score = DBA_SCORE_READS)

counts <- dba.peakset(dbo, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.table(counts, args[2], sep = "\t", quote = FALSE, row.names = FALSE)

sessionInfo()
