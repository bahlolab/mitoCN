###### Nature mtDNA-CN GWAS ######
rm(list = ls())
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(ggpubr) #for multiple plots

##### adj mtcn -> PD
mtcn_adj_file <- "./Documents/mtCN_adj_exposure_Gupta.txt"
mtcn_adj_dat <- read_exposure_data(filename = mtcn_adj_file)
mtcn_adj_dat_LD <- clump_data(mtcn_adj_dat)
outcome_dat <- extract_outcome_data(snps=mtcn_adj_dat_LD$SNP, outcomes = "ieu-b-7")
dat <- harmonise_data(mtcn_adj_dat_LD, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/mtcn_adj_PD_Gupta.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + ylab("SNP effect on Parkinson's disease") + xlab("SNP effect on mtDNA-CN (WGS,adj)")
ggsave(p.mr, file = "./Documents/mtcn_adj_PD_Gupta.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: mtDNA-CN (WGS,adj) on PD")
ggsave(p.snp, file = "./Documents/mtcn_adj_PD_singleSNP_Gupta.png", width = 7, height = 15)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out: mtDNA-CN (WGS,adj) on PD")
ggsave(p.loo, file = "./Documents/mtcn_adj_PD_LOO_Gupta.png", width = 7, height = 15)
# Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/mtcn_adj_PD_Funnel_Gupta.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp2.png", plot = p_both, dpi = 500, width = 10, height = 10)

##### raw mtcn -> PD
mtcn_raw_file <- "./Documents/mtCN_raw_exposure_Gupta.txt"
mtcn_raw_dat <- read_exposure_data(filename = mtcn_raw_file)
mtcn_raw_dat_LD <- clump_data(mtcn_raw_dat)
outcome_dat <- extract_outcome_data(snps=mtcn_raw_dat_LD$SNP, outcomes = "ieu-b-7")
dat <- harmonise_data(mtcn_raw_dat_LD, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/mtcn_raw_PD_Gupta.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + ylab("SNP effect on Parkinson's disease") + xlab("SNP effect on mtDNA-CN (WGS,raw)")
ggsave(p.mr, file = "./Documents/mtcn_raw_PD_Gupta.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: mtDNA-CN (WGS,raw) on PD")
ggsave(p.snp, file = "./Documents/mtcn_raw_PD_singleSNP_Gupta.png", width = 7, height = 15)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out: mtDNA-CN (WGS,raw) on PD")
ggsave(p.loo, file = "./Documents/mtcn_raw_PD_LOO_Gupta.png", width = 7, height = 15)
# Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/mtcn_raw_PD_Funnel_Gupta.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp1.png", plot = p_both, dpi = 500, width = 10, height = 10)

##### PD -> mtcn adj
exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
outcome_dat <- read_outcome_data(filename = "./Documents/adj_mtcn_outcome.txt")
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/PD_mtcn_adj_MR_Gupta.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + xlab("SNP effect on Parkinson's disease") + ylab("SNP effect on mtDNA-CN (WGS,adj)")
ggsave(p.mr, file = "./Documents/PD_mtcn_adj_MR_Gupta.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: PD on mtDNA-CN (WGS,adj)")
ggsave(p.snp, file = "./Documents/PD_mtcn_adj_singleSNP_Gupta.png", width = 7, height = 7)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out:  PD on mtDNA-CN (WGS,adj)")
ggsave(p.loo, file = "./Documents/PD_mtcn_adj_LOO_Gupta.png", width = 7, height = 7)
# Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/PD_mtcn_Funnel_Gupta.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp6.png", plot = p_both, dpi = 500, width = 10, height = 10)

##### PD -> mtcn raw
exposure_dat <- extract_instruments(outcomes = "ieu-b-7")
exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
outcome_dat <- read_outcome_data(filename = "./Documents/raw_mtcn_outcome.txt")
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/PD_mtcn_raw_MR_Gupta.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + xlab("SNP effect on Parkinson's disease") + ylab("SNP effect on mtDNA-CN (WGS,raw)")
ggsave(p.mr, file = "./Documents/PD_mtcn_raw_MR_Gupta.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: PD on mtDNA-CN (WGS,raw)")
ggsave(p.snp, file = "./Documents/PD_mtcn_raw_singleSNP_Gupta.png", width = 7, height = 7)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out:  PD on mtDNA-CN (WGS,raw)")
ggsave(p.loo, file = "./Documents/PD_mtcn_raw_LOO_Gupta.png", width = 7, height = 7)
# Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/PD_mtcn_Funnel_Gupta.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp5.png", plot = p_both, dpi = 500, width = 10, height = 10)
