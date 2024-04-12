###### data processing ######
rm(list = ls())
library(data.table)
library(TwoSampleMR)
library(MRInstruments)
library(readxl)

GRS_all <- read_excel("/Users/wang.lo/Downloads/Nature_Supp.xlsx", sheet = "Table 2", skip = 1)
adj <- GRS_all[GRS_all$phenotype_id == "mtcn",]
raw <- GRS_all[GRS_all$phenotype_id == "raw_mtcn",]

##### adj mtcn exposure
adj$chr <- unlist(strsplit(adj$variant_b37,":"))[1+(0:91)*4]
adj$position <- unlist(strsplit(adj$variant_b37,":"))[2+(0:91)*4]
adj$other_allele <- unlist(strsplit(adj$variant_b37,":"))[3+(0:91)*4]
adj$effect_allele <- unlist(strsplit(adj$variant_b37,":"))[4+(0:91)*4]
adj$eaf <- 1-adj$AF_Allele2
adj$beta <- (-1) * adj$beta
adj_mtCN_GRS <- adj[,c("rsid","beta","se","effect_allele","other_allele","eaf","p","chr","position","N")]
colnames(adj_mtCN_GRS) <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval","chr","position","samplesize")
adj_mtCN_GRS$Phenotype <- "mtcn_adj"
write.table(adj_mtCN_GRS, "./Documents/adj_mtCN_exposure.txt", row.names = F, quote = F)

##### raw mtcn exposure
raw$chr <- unlist(strsplit(raw$variant_b37,":"))[1+(0:140)*4]
raw$position <- unlist(strsplit(raw$variant_b37,":"))[2+(0:140)*4]
raw$other_allele <- unlist(strsplit(raw$variant_b37,":"))[3+(0:140)*4]
raw$effect_allele <- unlist(strsplit(raw$variant_b37,":"))[4+(0:140)*4]
raw$eaf <- 1-raw$AF_Allele2
raw$beta <- (-1) * raw$beta
raw_mtCN_GRS <- raw[,c("rsid","beta","se","effect_allele","other_allele","eaf","p","chr","position","N")]
colnames(raw_mtCN_GRS) <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval","chr","position","samplesize")
raw_mtCN_GRS$Phenotype <- "mtcn_raw"
write.table(raw_mtCN_GRS, "./Documents/raw_mtCN_exposure.txt", row.names = F, quote = F)

exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
exposure_dat$snp_name <- paste(exposure_dat$chr.exposure,exposure_dat$pos.exposure, sep = "_")
##### adj mtcn outcome
adj_GRS_file <- "/Users/wang.lo/Downloads/GCST90268497.tsv.gz" #hg19
adj0 <- fread(adj_GRS_file)
adj0$snp_name <- paste(adj0$chromosome,adj0$base_pair_location, sep = "_")
adj1 <- adj0[adj0$snp_name %in% exposure_dat$snp_name,]
adj <- merge(exposure_dat[,c("SNP","snp_name")], adj1, by = "snp_name")
adj$eaf <- 1 - adj$effect_allele_frequency
adj$beta <- (-1) * adj$beta
adj_outcome <- adj[,c("SNP","beta","standard_error","effect_allele","other_allele","eaf","p_value","chromosome","base_pair_location")]
colnames(adj_outcome) <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval","chr","position")
write.table(adj_outcome, "./Documents/adj_mtcn_outcome.txt", row.names = F, quote = F)

##### raw mtcn outcome
raw_GRS_file <- "/Users/wang.lo/Downloads/GCST90268498.tsv.gz" #hg19
raw0 <- fread(raw_GRS_file)
raw0$snp_name <- paste(raw0$chromosome,raw0$base_pair_location, sep = "_")
raw1 <- raw0[raw0$snp_name %in% exposure_dat$snp_name,]
raw <- merge(exposure_dat[,c("SNP","snp_name")], raw1, by = "snp_name")
raw$eaf <- 1 - raw$effect_allele_frequency
raw$beta <- (-1) * raw$beta
raw_outcome <- raw[,c("SNP","beta","standard_error","effect_allele","other_allele","eaf","p_value","chromosome","base_pair_location")]
colnames(raw_outcome) <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval","chr","position")
write.table(raw_outcome, "./Documents/raw_mtcn_outcome.txt", row.names = F, quote = F)

