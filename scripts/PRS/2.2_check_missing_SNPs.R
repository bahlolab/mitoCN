### check missing SNPs ###
rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/PRS/")
library(readxl)

adj_bim <- read.table("./results/mtCN/mtCN_GRS_loci.bim")
raw_bim <- read.table("./results/mtCN/raw_mtCN_GRS_loci.bim")
GRS_all <- read_excel("GRS/mtCN_PRS/Nature_Supp.xlsx", sheet = "Table 2", skip = 1)
GRS_all$beta <- (-1)*GRS_all$beta
adj_mtcn <- GRS_all[GRS_all$phenotype_id == "mtcn",]
raw_mtcn <- GRS_all[GRS_all$phenotype_id == "raw_mtcn",]
ms_adj <- mtcn[!(mtcn$rsid %in% adj_bim$V2),]
ms_raw <- raw_mtcn[!(raw_mtcn$rsid %in% raw_bim$V2),]
write.table(ms_adj, "./GRS/mtCN_PRS/ms_adj_mtcn.tsv", row.names = F, quote = F, sep = "\t")
write.table(ms_raw, "./GRS/mtCN_PRS/ms_raw_mtcn.tsv", row.names = F, quote = F, sep = "\t")

ms_adj$CHR <- sub("...","",unlist(strsplit(ms_adj$variant_b38,":"))[4*(0:(nrow(ms_adj)-1))+1])
ms_adj$BP1 <- as.numeric(unlist(strsplit(ms_adj$variant_b38,":"))[4*(0:(nrow(ms_adj)-1))+2]) - 10
ms_adj$BP2 <- as.numeric(unlist(strsplit(ms_adj$variant_b38,":"))[4*(0:(nrow(ms_adj)-1))+2]) + 10
ms_adj$LABEL <- paste(unlist(strsplit(ms_adj$variant_b38,":"))[4*(0:(nrow(ms_adj)-1))+3],
                       unlist(strsplit(ms_adj$variant_b38,":"))[4*(0:(nrow(ms_adj)-1))+4], sep = "_")
write.table(ms_adj[,13:16], "./GRS/mtCN_PRS/missing_adj_mtcn.txt", quote = F, row.names = F, col.names = F, sep = "\t")

ms_raw$CHR <- sub("...","",unlist(strsplit(ms_raw$variant_b38,":"))[4*(0:(nrow(ms_raw)-1))+1])
ms_raw$BP1 <- as.numeric(unlist(strsplit(ms_raw$variant_b38,":"))[4*(0:(nrow(ms_raw)-1))+2]) - 10
ms_raw$BP2 <- as.numeric(unlist(strsplit(ms_raw$variant_b38,":"))[4*(0:(nrow(ms_raw)-1))+2]) + 10
ms_raw$LABEL <- paste(unlist(strsplit(ms_raw$variant_b38,":"))[4*(0:(nrow(ms_raw)-1))+3],
                       unlist(strsplit(ms_raw$variant_b38,":"))[4*(0:(nrow(ms_raw)-1))+4], sep = "_")
write.table(ms_raw[,13:16], "./GRS/mtCN_PRS/missing_raw_mtcn.txt", quote = F, row.names = F, col.names = F, sep = "\t")

# extract each SNP from plink file
# if not exist, use the LDProxy tool at https://ldlink.nih.gov/?tab=ldproxy
# select the SNP with D'=1 and check if it passed the filter
# if not, check the next SNP with D'=1

### improve using https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html