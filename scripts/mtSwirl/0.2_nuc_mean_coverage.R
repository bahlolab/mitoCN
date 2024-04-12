args <- commandArgs(trailingOnly = TRUE)
final_stats_table <- args[1]
summary_stat <- args[2]
final_res <- args[3]

# final_stats_table <- "/vast/scratch/users/wang.lo/mtSwirl/output/res/BioFIND/BF-1001/BF-1001.stat.summary.txt"
# summary_stat <- "/vast/scratch/users/wang.lo/mtSwirl/output/res/BioFIND/BF-1001/BF-1001_mtanalysis_diagnostic_statistics.tsv"

df <- read.delim(final_stats_table, sep = "\t", header = T)[1,]
stat <- read.delim(summary_stat, sep = "\t", header = T)

nuc_mean_coverage <- ((df$total_mapped_reads - df$singletons - df$mate_diff_chr - df$duplicates)*as.numeric(df$READ_LENGTH))/df$genome_length
mt_mean_coverage <- stat$mean_coverage

mtcn <- round(2 * mt_mean_coverage / nuc_mean_coverage, 2)

res <- c(round(nuc_mean_coverage,2), mt_mean_coverage, mtcn)
names(res) <- c("nuc_mean_coverage", "mt_mean_coverage", "mtcn")
mtcn <- do.call(data.frame, as.list(res))

write.table(mtcn, sep ='\t', row.names = F, file = final_res, quote = F)
