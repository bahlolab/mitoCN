rm(list = ls())
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github('MRCIEU/MRInstruments') 
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)

### PD -> mtDNA-CN
# Get instruments
exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "ebi-a-GCST90026372")
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
# remove duplicates
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/PD_mtcn_MR_Chong.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + xlab("SNP effect on Parkinson's disease") + ylab("SNP effect on mtDNA-CN (genotyping)")
ggsave(p.mr, file = "./Documents/PD_mtcn_MR_Chong.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: PD on mtDNA-CN (genotyping)")
ggsave(p.snp, file = "./Documents/PD_mtcn_singleSNP_Chong.png", width = 7, height = 7)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out:  PD on mtDNA-CN (genotyping)")
ggsave(p.loo, file = "./Documents/PD_mtcn_LOO_Chong.png", width = 7, height = 7)
# Funnel plot
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/PD_mtcn_Funnel_Chong.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp7.png", plot = p_both, dpi = 500, width = 10, height = 10)

### mtDNA-CN -> PD
# Get instruments
exposure_dat <- extract_instruments(outcomes = "ebi-a-GCST90026372")
# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "ieu-b-7")
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/mtcn_PD_Chong.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + ylab("SNP effect on Parkinson's disease") + xlab("SNP effect on mtDNA-CN (genotyping)")
ggsave(p.mr, file = "./Documents/mtcn_PD_Chong.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: mtDNA-CN (genotyping) on PD")
ggsave(p.snp, file = "./Documents/mtcn_PD_singleSNP_Chong.png", width = 7, height = 15)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out: mtDNA-CN (genotyping) on PD")
ggsave(p.loo, file = "./Documents/mtcn_PD_LOO_Chong.png", width = 7, height = 15)
# Funnel plot
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/mtcn_PD_Funnel_Chong.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp3.png", plot = p_both, dpi = 500, width = 10, height = 10)
