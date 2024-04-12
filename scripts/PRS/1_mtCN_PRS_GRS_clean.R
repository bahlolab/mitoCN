##### mtCN_adj GRS clean #####
rm(list = ls())
library(readxl)
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/PRS/")

GRS_all <- read_excel("GRS/mtCN_PRS/Nature_Supp.xlsx", sheet = "Table 2", skip = 1)
adj <- GRS_all[GRS_all$phenotype_id == "mtcn",]
raw <- GRS_all[GRS_all$phenotype_id == "raw_mtcn",]
head(adj)

adj$ALT_ALLELE <- sub(".*:", "", adj$variant_b37)
adj_mtCN_GRS <- adj[,c("rsid","ALT_ALLELE","beta")]
## correct alt allele swap issue in the GRS file
adj_mtCN_GRS$beta <- (-1)*adj_mtCN_GRS$beta
colnames(adj_mtCN_GRS) <- c("rsID",	"ALT_ALLELE",	"BETA")
write.table(adj_mtCN_GRS, "./GRS/mtCN_PRS/adj_mtCN_GRS.txt", row.names = F, quote = F, sep = "\t")

raw$ALT_ALLELE <- sub(".*:", "", raw$variant_b37)
raw_mtCN_GRS <- raw[,c("rsid","ALT_ALLELE","beta")]
raw_mtCN_GRS$beta <- (-1)*raw_mtCN_GRS$beta
colnames(raw_mtCN_GRS) <- c("rsID",	"ALT_ALLELE",	"BETA")
write.table(raw_mtCN_GRS, "./GRS/mtCN_PRS/raw_mtCN_GRS.txt", row.names = F, quote = F, sep = "\t")
