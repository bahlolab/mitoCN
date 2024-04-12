rm(list=ls())
library(sfsmisc)
library(MASS)
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")

# load mtDNA-CN results
load(file = "./meta_analysis/mt.all.RData")
# exclude LBD cohort
mt1 <- mt[mt$study != "LBD",]
# remove outliers from HBS
HBS <- mt1[mt1$study == "HBS",]
c2 <- mt1$logmt[mt1$cluster == 2 & mt1$study.cluster != "HBS"]
HBS$z <- (HBS$logmt - mean(c2)) / sd(c2)
outlier <- HBS$ID[HBS$z > 3 | HBS$z < (-3)]
mt2 <- mt1[!(mt1$ID %in% outlier),c("ID","mt","logmt","study.cluster","study")]

# load metadata
metadata <- read.csv("./metadata/v3/metadata_20230320.csv")[,c("participant_id","sex","age_at_baseline")]
eigenVectors <- read.table("./meta_analysis/PCA/results/plinkPCA_all.eigenvec")[,c(1,3:7)]
colnames(eigenVectors) <- c("ID",paste0("PC",1:5))
meta <- merge(metadata, eigenVectors, by.x = "participant_id", by.y = "ID")
# merge data
data <- merge(meta, mt2, by.x = "participant_id", by.y = "ID")

# association test
asso <- list()
asso$PPMI1 <- data[data$study.cluster == "PPMI.1",]
asso$PDBP1 <- data[data$study.cluster == "PDBP.1",]
asso$BioFIND1 <- data[data$study.cluster == "BioFIND.1",]
asso$PPMI2 <- data[data$study.cluster == "PPMI.2",]
asso$PDBP2 <- data[data$study.cluster == "PDBP.2",]
asso$HBS2 <- data[data$study.cluster == "HBS.2",]
asso$SURE2 <- data[data$study.cluster == "SURE-PD.2",]
asso$STEADY2 <- data[data$study.cluster == "STEADY-PD3.2",]

##### sex
test = asso[[1]]
res <- summary(lm.sex <- rlm(formula = logmt ~ sex + PC1 + PC2 + PC3 + PC4 + PC5, data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.sex, var = "sexMale")$p.value
res$n.samp <- nrow(test)
res$female <- as.numeric(table(test$sex)[1])
res$male <- as.numeric(table(test$sex)[2])

for(i in 2:8){
  test = asso[[i]]
  res_new <- summary(lm.sex <- rlm(formula = logmt ~ sex + PC1 + PC2 + PC3 + PC4 + PC5, data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(lm.sex, var = "sexMale")$p.value
  res_new$n.samp <- nrow(test)
  res_new$female <- as.numeric(table(test$sex)[1])
  res_new$male <- as.numeric(table(test$sex)[2])
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")
res$prop <- exp(as.numeric(res$Estimate))

write.table(res, "./meta_analysis/res/sex.rlm.tsv", quote = F, row.names = F, sep = "\t")

##### age
test = asso[[1]]
res <- summary(rlm.age<- rlm(formula = logmt ~ age_at_baseline + PC1 + PC2 + PC3 + PC4 + PC5, data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(rlm.age, var = "age_at_baseline")$p.value
res$n.samp <- nrow(test)
res$age_mean <- mean(test$age_at_baseline)
res$age_sd <- sd(test$age_at_baseline)

for(i in 2:8){
  test = asso[[i]]
  res_new <- summary(rlm.age<- rlm(formula = logmt ~ age_at_baseline + PC1 + PC2 + PC3 + PC4 + PC5, data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(rlm.age, var = "age_at_baseline")$p.value
  res_new$n.samp <- nrow(test)
  res_new$age_mean <- mean(test$age_at_baseline)
  res_new$age_sd <- sd(test$age_at_baseline)
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")

write.table(res, "./meta_analysis/res/age.rlm.tsv", quote = F, row.names = F, sep = "\t")


