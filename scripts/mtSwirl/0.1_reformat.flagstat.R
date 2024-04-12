args <- commandArgs(trailingOnly = TRUE)
flagstat.pre <- args[1]
flagstat.fmt <- args[2]

# flagstat.pre <- "/vast/scratch/users/wang.lo/mtSwirl/output/PP-10874/PP-10874.flagstat.txt"
# flagstat.fmt <- "/vast/scratch/users/wang.lo/mtSwirl/output/PP-10874/PP-10874.flagstat.tmp.txt"

vec <- readLines(flagstat.pre)
titles <- c('total', 'secondary', 'supplementary', 'duplicates', 'mapped', 'paired', 'read1', 'read2', 'properly_paired', 'with_itself_and_mate_mapped', 'singletons', 'mate_diff_chr', 'mate_diff_chr_mapq_5')
get_ele <- function(x) gregexpr('^[0-9]+',x)[[1]]
results_vec <- as.numeric(sapply(vec, function(x) substr(x, get_ele(x)[1], get_ele(x)[1] + attr(get_ele(x), 'match.length') - 1)))
names(results_vec) <- titles
df <- do.call(data.frame, as.list(results_vec))
write.table(df, sep ='\t', row.names = F, file = flagstat.fmt, quote = F)
