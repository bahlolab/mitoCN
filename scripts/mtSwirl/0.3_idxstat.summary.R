args <- commandArgs(trailingOnly = TRUE)
idxstats <- args[1]
idx_summary <- args[2]

# idxstats <- "/vast/scratch/users/wang.lo/mtSwirl/output/PP-10874/PP-10874.stats.tsv"
# idx_summary <- "/vast/scratch/users/wang.lo/mtSwirl/output/PP-10874/PP-10874.idxstat.txt"
  
vec <- read.csv(idxstats, sep='\t', stringsAsFactors=F, col.names=c('chr','len','mapped_reads','unmapped_reads'), header=F)[1:22,]
total_mapped_reads <- sum(vec$mapped_reads)
total_unmapped_reads <- sum(vec$unmapped_reads)
genome_length <- sum(vec$len)

results_vec <- as.numeric(c(total_mapped_reads, total_unmapped_reads, genome_length))
names(results_vec) <- c("total_mapped_reads", "total_unmapped_reads", "genome_length")
df <- do.call(data.frame, as.list(results_vec))

write.table(df, sep ='\t', row.names = F, file = idx_summary, quote = F)
