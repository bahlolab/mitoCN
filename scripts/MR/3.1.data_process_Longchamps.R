###### data processing ######
rm(list = ls())
library(data.table)
library(TwoSampleMR)
library(MRInstruments)

### WGS + SNPchip
WES_file <- "/Users/wang.lo/Downloads/mtDNA_CN.ALLm2.bgen.stats.gz"
WES <- fread(WES_file)

##### mtcn exposure
WES <- WES[,c("SNP","BETA","SE","ALLELE1","ALLELE0","A1FREQ","P_LINREG","CHR","BP")]
colnames(WES) <- c("SNP","beta","se","effect_allele","other_allele","eaf","pval","chr","position")
WES$Phenotype <- "mtDNA-CN"
WES.p <- WES[WES$pval < 5e-8,]
write.table(WES.p, "./Documents/mtcn_WES_exposure.txt", row.names = F, quote = F)

exposure_dat <- extract_instruments(outcomes = "ieu-b-7", p1 = 5e-06, p2 = 5e-06)
exposure_dat$snp_name <- paste(exposure_dat$chr.exposure,exposure_dat$pos.exposure, sep = "_")
##### mtcn outcome
WES$snp_name <- paste(WES$chr,WES$position, sep = "_")
WES1 <- WES[WES$snp_name %in% exposure_dat$snp_name,]
write.table(WES1, "./Documents/mtcn_WES_outcome.txt", row.names = F, quote = F)
