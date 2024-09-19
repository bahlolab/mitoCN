rm(list=ls())
library(magrittr)
library(sfsmisc)
library(MASS)
library(ggplot2)
library(ggstatsplot)

setwd("/vast/projects/bahlo_amppd/metadata/")
input_path <- "/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/"
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
metadata <- read.csv(paste0(input_path,"metadata/v3/metadata_20230525.csv"))
meta <- metadata[,c("participant_id","case_control_other_at_baseline","age_at_baseline","sex")]
PD.meta <- metadata[metadata$study == "PDBP",]

### load cell composition estimates
PD <- read.table(paste0(input_path,"meta_analysis/cell_composition/results/CIBERSORTx_PD_CPM.txt"), header = T, sep = "\t")
convertID <- function(x) y <- paste(x[1:2], collapse = "-")
PD$ID <- unlist(lapply(strsplit(PD$Mixture,"[.]"), convertID))

### load ancestry
ancestry <- read.table(paste0(input_path,"ancestry/version2/PC_table_with_estimated_superpop_PC123.txt"), sep = "\t",header = T)[,c("sample",paste0("PC",1:5))]

### load medication history
med <- read.csv(paste0(input_path,"metadata/v3/clinical/PD_Medical_History.csv"))
med.dic <- read.csv(paste0(input_path,"metadata/v3/clinical/PD_Medical_History_dictionary.csv"))
med <- merge(med, ancestry, by.x = "participant_id", by.y = "sample")

PD.med <- med[startsWith(med$participant_id,"PD") & med$visit_name == "M0" & med$diagnosis %in% c("Parkinson's Disease","No PD Nor Other Neurological Disorder"),
              c("participant_id","diagnosis","age_at_diagnosis","on_levodopa","on_dopamine_agonist","on_other_pd_medications", paste0("PC",1:5))]
PD.med <- PD.med[PD.med$diagnosis == "Parkinson's Disease",]
PD.med <- PD.med[PD.med$on_levodopa != "" & PD.med$on_dopamine_agonist != "" & PD.med$on_other_pd_medications != "",]
columns_to_factor <- c("diagnosis","on_levodopa","on_dopamine_agonist","on_other_pd_medications")
PD.med[columns_to_factor] <- lapply(PD.med[columns_to_factor], factor)

summary(PD.med)

### PDBP
PD.blood <- process(PD)
PD1 <- merge(PD.med, meta, by = "participant_id", all.x = T)
PD2 <- merge(PD1, PD.blood, by.x = "participant_id", by.y = "ID")
columns_to_factor <- c("case_control_other_at_baseline","sex")
PD2[columns_to_factor] <- lapply(PD2[columns_to_factor], factor)

### clinic
load(file = "/vast/scratch/users/wang.lo/PDBP/data.RData")
clinic <- data[,-(2:9)]

data <- merge(PD2, clinic, by = "participant_id")
summary(data)

### association tests: rlm, blood marker ~ age + sex + PC1-5 + medication
# NLR
res1 <- summary(m1 <- rlm(formula = NLR ~ age_at_baseline + sex + on_levodopa + on_dopamine_agonist + on_other_pd_medications + 
                            PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2:7,] %>% as.data.frame
p.values <- numeric()  
for (i in 2:7) {  # Adjust the range as needed based on your model
  test_result <- f.robftest(m1, var = i)
  p.values <- c(p.values, test_result$p.value)
}
res1$p.value <- p.values

# Neutrophils
res2 <- summary(m1 <- rlm(formula = Neutrophils.p ~ age_at_baseline + sex + on_levodopa + on_dopamine_agonist + on_other_pd_medications + 
                            PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2:7,] %>% as.data.frame
p.values <- numeric()  
for (i in 2:7) {  # Adjust the range as needed based on your model
  test_result <- f.robftest(m1, var = i)
  p.values <- c(p.values, test_result$p.value)
}
res2$p.value <- p.values

# Lymphocytes
res3 <- summary(m1 <- rlm(formula = Lymphocytes.p ~ age_at_baseline + sex + on_levodopa + on_dopamine_agonist + on_other_pd_medications +
                            PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2:7,] %>% as.data.frame
p.values <- numeric()  
for (i in 2:7) {  # Adjust the range as needed based on your model
  test_result <- f.robftest(m1, var = i)
  p.values <- c(p.values, test_result$p.value)
}
res3$p.value <- p.values

res <- cbind(res1, res2, res3)
write.table(res, "blood_meds.tsv", quote = F, sep = "\t")

##### test if medication influence on the observed associations
### 1. association tests: rlm, NLR ~ PD severity + age + sex + PC1-5
### 2. association tests: rlm, NLR ~ PD severity + age + sex + PC1-5 + medication

res0 <- summary(m0 <- rlm(formula = NLR ~ age_at_baseline + sex + mds_updrs_part_iii_summary_score + 
                            PC1 + PC2 + PC3 + PC4 + PC5, data = data))$coefficients[2:4,] %>% as.data.frame
p.values <- numeric()  
for (i in 2:4) {  # Adjust the range as needed based on your model
  test_result <- f.robftest(m0, var = i)
  p.values <- c(p.values, test_result$p.value)
}
res0$p.value <- p.values

### without medication
test <- data[,c(7:11,13:14,17,18)]
res <- summary(m1 <- rlm(formula = NLR ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(m1, var = 9)$p.value

for(i in 19:25){
  test <- data[,c(7:11,13:14,17,i)]
  res_new <- summary(m1 <- rlm(formula = NLR ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(m1, var = 9)$p.value
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$PD.vars <- colnames(data)[c(18:25)]
res$p.adj <- p.adjust(res$`Pr(>|t|)`, method="fdr")

write.table(res, "/vast/scratch/users/wang.lo/PDBP/PDBP.tsv", quote = F, row.names = F, sep = "\t")

### adjusted mtDNA-CN
test <- data[,c(4:6,7:11,13:14,17,18)]
res <- summary(m1 <- rlm(formula = NLR ~ ., data = test))$coefficients[12,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(m1, var = 12)$p.value

for(i in 19:25){
  test <- data[,c(4:6,7:11,13:14,17,i)]
  res_new <- summary(m1 <- rlm(formula = NLR ~ ., data = test))$coefficients[12,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(m1, var = 12)$p.value
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$PD.vars <- colnames(data)[c(18:25)]
res$p.adj <- p.adjust(res$`Pr(>|t|)`, method="fdr")
res$n.samp <- apply(data[,c(18:25)],2, function(x) sum(!(is.na(x))))
write.table(res, "/vast/scratch/users/wang.lo/PDBP/PDBP.med.tsv", quote = F, row.names = F, sep = "\t")
