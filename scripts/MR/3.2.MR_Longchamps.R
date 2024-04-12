###### Longchamps mtDNA-CN GWAS ######
rm(list = ls())
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)

##### mtcn -> PD
mtcn_file <- "./Documents/mtCN_WES_exposure.txt"
mtcn_dat <- read_exposure_data(filename = mtcn_file)
mtcn_dat_LD <- clump_data(mtcn_dat)
outcome_dat <- extract_outcome_data(snps=mtcn_dat_LD$SNP, outcomes = "ieu-b-7")
dat <- harmonise_data(mtcn_dat_LD, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/mtcn_PD_Longchamps.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + ylab("SNP effect on Parkinson's disease") + xlab("SNP effect on mtDNA-CN (WES+genotyping)")
ggsave(p.mr, file = "./Documents/mtcn_PD_Longchamps.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: mtDNA-CN (WES+genotyping) on PD")
ggsave(p.snp, file = "./Documents/mtcn_PD_singleSNP_Longchamps.png", width = 7, height = 15)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out: mtDNA-CN (WES+genotyping) on PD")
ggsave(p.loo, file = "./Documents/mtcn_PD_LOO_Longchamps.png", width = 7, height = 15)
# Funnel plot
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/mtcn_PD_Funnel_Longchamps.png", width = 7, height = 7)

ggarrange(p.mr, p.funnel, ncol = 2, labels = c("A", "B"))
ggarrange(p.snp, p.loo, ncol = 2, labels = c("C", "D"))

p_both <- ggarrange(ggarrange(p.mr, p.funnel, ncol = 2, labels = c("A", "B")),
                    ggarrange(p.snp, p.loo, ncol = 2, labels = c("C", "D")),
                    nrow = 2, heights = c(1, 2))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp4.png", plot = p_both, dpi = 500, width = 10, height = 12)

##### PD -> mtcn
exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
outcome_dat <- read_outcome_data(filename = "./Documents/mtcn_WES_outcome.txt")
# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
dat <- dat[!duplicated(dat$SNP),]
# Perform MR
res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
write.table(res,"./Documents/PD_mtcn_MR_Longchamps.tsv", row.names = F, quote = F, sep = "\t")
p1 <- mr_scatter_plot(res, dat)
p.mr <- p1[[1]] + xlab("SNP effect on Parkinson's disease") + ylab("SNP effect on mtDNA-CN (WES+genotyping)")
ggsave(p.mr, file = "./Documents/PD_mtcn_MR_Longchamps.png", width = 7, height = 7)
## Sensitivity analyses
# Heterogeneity statistics
mr_heterogeneity(dat)
# Horizontal pleiotropy
mr_pleiotropy_test(dat)
# Single SNP analysis
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p.snp <- p2[[1]] + xlab("MR effect size: PD on mtDNA-CN (WES+genotyping)")
ggsave(p.snp, file = "./Documents/PD_mtcn_singleSNP_Longchamps.png", width = 7, height = 7)
# Leave-one-out analysis
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p.loo <- p3[[1]] + xlab("MR leave-one-out:  PD on mtDNA-CN (WES+genotyping)")
ggsave(p.loo, file = "./Documents/PD_mtcn_LOO_Longchamps.png", width = 7, height = 7)
# Funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p.funnel <- p4[[1]]
ggsave(p.funnel, file = "./Documents/PD_mtcn_Funnel_Longchamps.png", width = 7, height = 7)

p_both <- ggarrange(p.mr, p.funnel, p.snp, p.loo, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(filename = "/Users/wang.lo/Documents/MR.supp8.png", plot = p_both, dpi = 500, width = 10, height = 10)
