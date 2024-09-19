rm(list=ls())
library(magrittr)
library(sfsmisc)
library(MASS)
library(ggplot2)
library(ggstatsplot)

setwd("/vast/projects/bahlo_amppd/PPMI/analysis/")
input_path <- "/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/"

### genetically defined ancestry
ancestry <- read.table(paste0(input_path,"ancestry/version2/PC_table_with_estimated_superpop_PC123.txt"), sep = "\t",header = T)
AFR <- ancestry[ancestry$AFR == "AFR",c("sample",paste0("PC",1:5))]

### load cell composition estimates
BF <- read.table(paste0(input_path,"meta_analysis/cell_composition/results/CIBERSORTx_BF_CPM.txt"), header = T, sep = "\t")
PD <- read.table(paste0(input_path,"meta_analysis/cell_composition/results/CIBERSORTx_PD_CPM.txt"), header = T, sep = "\t")
PP <- read.table(paste0(input_path,"meta_analysis/cell_composition/results/CIBERSORTx_PP_CPM.txt"), header = T, sep = "\t")
res <- rbind(BF, PD, PP)
convertID <- function(x) y <- paste(x[1:2], collapse = "-")
res$ID <- unlist(lapply(strsplit(res$Mixture,"[.]"), convertID))
### function: data pre-processing
process <- function(data, percent = TRUE){
  data$B.cells <- rowSums(data[,startsWith(colnames(data),"B.cells")])
  data$NK.cells <- rowSums(data[,startsWith(colnames(data),"NK.cells")])
  data$T.cells <- rowSums(data[,startsWith(colnames(data),"T.cells")])
  data$Lymphocytes <- rowSums(data[,c("B.cells","T.cells","NK.cells")])
  data$Lymphocytes.p <- data$Lymphocyte/data$Absolute.score..sig.score.
  data$Neutrophils.p <- data$Neutrophils/data$Absolute.score..sig.score.
  data$NLR <- data$Neutrophils/data$Lymphocytes
  
  data1 <- data[,c("ID","Neutrophils.p","Lymphocytes.p","NLR")]
  return(data1)
}
amp <- process(res)
amp.afr <- amp[amp$ID %in% AFR$sample,]

metadata <- read.csv(paste0(input_path,"metadata/v3/metadata_20230525.csv"))
meta <- metadata[,c("participant_id","case_control_other_at_baseline","age_at_baseline","sex","study")]
data <- merge(amp.afr, meta, by.x = "ID", by.y = "participant_id", all.x = T)
data <- data[data$case_control_other_at_baseline != "Other",]
#! 48 individuals of African ancestry with WB RNA-seq data
columns_to_factor <- c("sex", "study")
data[columns_to_factor] <- lapply(data[columns_to_factor], factor)
data$case_control_other_at_baseline <- factor(data$case_control_other_at_baseline, levels = c("Control","Case"))
data <- merge(data, AFR, by.x = "ID", by.y = "sample")
### check if batch effect between studies: no difference
ggbetweenstats(
  data  = data[data$study != "BioFIND" & data$case_control_other_at_baseline == "Control",],
  x     = study,
  y     = NLR,
)

### association tests: rlm(blood marker) ~ diagnosis + age + sex
res <- summary(m1 <- rlm(formula = Neutrophils.p ~ case_control_other_at_baseline + age_at_baseline + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(m1, var = 2)$p.value

res1 <- summary(m1 <- rlm(formula = Lymphocytes.p ~ case_control_other_at_baseline + age_at_baseline + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res1$p.value <- f.robftest(m1, var = 2)$p.value
res <- rbind(res, res1)

res1 <- summary(m1 <- rlm(formula = NLR ~ case_control_other_at_baseline + age_at_baseline + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res1$p.value <- f.robftest(m1, var = 2)$p.value
res <- rbind(res, res1)

colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$blood.vars <- colnames(data)[2:4]
res$p.adj <- p.adjust(res$`Pr(>|t|)`, method = "fdr")
write.table(res, "./AMPPD.AFR.tsv", quote = F, row.names = F, sep = "\t")
