rm(list=ls())
library(magrittr)
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
metadata <- read.table("./meta_analysis/res/PD.vars.tsv", header = T, sep = "\t")
eigenVectors <- read.table("./meta_analysis/PCA/results/plinkPCA_all.eigenvec")[,c(1,3:7)]
colnames(eigenVectors) <- c("ID",paste0("PC",1:5))
meta <- merge(metadata, eigenVectors, by = "ID")
# merge data
meta.mt <- merge(meta, mt2, by = "ID")

# association test
data <- meta.mt[meta.mt$diagnosis != "Other",]
data$diagnosis <- factor(data$diagnosis, level = c("Control", "Case"))
summary(data)
asso <- list()
asso$PPMI1 <- data[data$study.cluster == "PPMI.1",]
asso$PDBP1 <- data[data$study.cluster == "PDBP.1",]
asso$BioFIND1 <- data[data$study.cluster == "BioFIND.1",]
asso$PPMI2 <- data[data$study.cluster == "PPMI.2",]
asso$PDBP2 <- data[data$study.cluster == "PDBP.2",]
asso$HBS2 <- data[data$study.cluster == "HBS.2",]
# asso$SURE2 <- data[data$study.cluster == "SURE-PD.2",]
# asso$STEADY2 <- data[data$study.cluster == "STEADY-PD3.2",]

##### diagnosis
test = asso[[1]]
res <- summary(lm.dx <- rlm(formula = logmt ~ diagnosis + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                            data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.dx, var = "diagnosisCase")$p.value
res$n.samp <- nrow(test)
res$Control <- paste0(round(mean(test$mt[test$diagnosis == "Control"]),2),
                      "(",round(sd(test$mt[test$diagnosis == "Control"]),2),")")
res$Case <- paste0(round(mean(test$mt[test$diagnosis == "Case"]),2),
                   "(",round(sd(test$mt[test$diagnosis == "Case"]),2),")")

for(i in 2:6){
  test = asso[[i]]
  res_new <- summary(rlm.pd <- rlm(formula = logmt ~ diagnosis + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                                   data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(rlm.pd, var = "diagnosisCase")$p.value
  res_new$n.samp <- nrow(test)
  res_new$Control <- paste0(round(mean(test$mt[test$diagnosis == "Control"]),2),
                        "(",round(sd(test$mt[test$diagnosis == "Control"]),2),")")
  res_new$Case <- paste0(round(mean(test$mt[test$diagnosis == "Case"]),2),
                     "(",round(sd(test$mt[test$diagnosis == "Case"]),2),")")
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS")

write.table(res, "./meta_analysis/res/PD.rlm.tsv", quote = F, row.names = F, sep = "\t")

#################################################################################
##### clinical assessments
data <- meta.mt
summary(data)
asso <- list()
asso$PPMI1 <- data[data$study.cluster == "PPMI.1",]
asso$PDBP1 <- data[data$study.cluster == "PDBP.1",]
asso$BioFIND1 <- data[data$study.cluster == "BioFIND.1",]
asso$PPMI2 <- data[data$study.cluster == "PPMI.2",]
asso$PDBP2 <- data[data$study.cluster == "PDBP.2",]
asso$HBS2 <- data[data$study.cluster == "HBS.2",]
asso$SURE2 <- data[data$study.cluster == "SURE-PD.2",]
asso$STEADY2 <- data[data$study.cluster == "STEADY-PD3.2",]

##### MDS_UPDRS_III
test = asso[[1]]
res <- summary(lm.dx <- rlm(formula = logmt ~ MDS_UPDRS_III + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                            data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.dx, var = "MDS_UPDRS_III")$p.value
res$n.samp <- sum(!(is.na(test$MDS_UPDRS_III)))

for(i in 2:8){
  test = asso[[i]]
  if(all(is.na(test$MDS_UPDRS_III))){
    res_new <- matrix(rep(NA, 3), nrow=1) %>% as.data.frame
    res_new$p.value <- NA
  }else{
    res_new <- summary(rlm.pd <- rlm(formula = logmt ~ MDS_UPDRS_III + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                                     data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
    res_new$p.value <- f.robftest(rlm.pd, var = "MDS_UPDRS_III")$p.value
  }
  res_new$n.samp <- sum(!(is.na(test$MDS_UPDRS_III)))
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")
res <- res[res$n.samp>0,]
write.table(res, "./meta_analysis/res/MDS_UPDRS_III.rlm.tsv", quote = F, row.names = F, sep = "\t")

#################################################################################
##### ADL
test = asso[[1]]
res <- summary(lm.dx <- rlm(formula = logmt ~ ADL + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                            data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.dx, var = "ADL")$p.value
res$n.samp <- sum(!(is.na(test$ADL)))

for(i in 2:8){
  test = asso[[i]]
  if(all(is.na(test$ADL))){
    res_new <- matrix(rep(NA, 3), nrow=1) %>% as.data.frame
    res_new$p.value <- NA
  }else{
    res_new <- summary(rlm.pd <- rlm(formula = logmt ~ ADL + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                                     data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
    res_new$p.value <- f.robftest(rlm.pd, var = "ADL")$p.value
  }
  res_new$n.samp <- sum(!(is.na(test$ADL)))
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")
res <- res[res$n.samp>0,]
write.table(res, "./meta_analysis/res/ADL.rlm.tsv", quote = F, row.names = F, sep = "\t")

#################################################################################
##### MoCA
test = asso[[1]]
res <- summary(lm.dx <- rlm(formula = logmt ~ MOCA + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                            data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.dx, var = "MOCA")$p.value
res$n.samp <- sum(!(is.na(test$MOCA)))

for(i in 2:8){
  test = asso[[i]]
  if(all(is.na(test$MOCA))){
    res_new <- matrix(rep(NA, 3), nrow=1) %>% as.data.frame
    res_new$p.value <- NA
  }else{
    res_new <- summary(rlm.pd <- rlm(formula = logmt ~ MOCA + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                                     data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
    res_new$p.value <- f.robftest(rlm.pd, var = "MOCA")$p.value
  }
  res_new$n.samp <- sum(!(is.na(test$MOCA)))
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")
res <- res[res$n.samp>0,]
write.table(res, "./meta_analysis/res/MOCA.rlm.tsv", quote = F, row.names = F, sep = "\t")

##### UPSIT
test = asso[[1]]
res <- summary(lm.dx <- rlm(formula = logmt ~ UPSIT + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                            data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(lm.dx, var = "UPSIT")$p.value
res$n.samp <- sum(!(is.na(test$UPSIT)))

for(i in 2:8){
  test = asso[[i]]
  if(all(is.na(test$UPSIT))){
    res_new <- matrix(rep(NA, 3), nrow=1) %>% as.data.frame
    res_new$p.value <- NA
  }else{
    res_new <- summary(rlm.pd <- rlm(formula = logmt ~ UPSIT + sex + age + PC1 + PC2 + PC3 + PC4 + PC5, 
                                     data = test))$coefficients[2,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
    res_new$p.value <- f.robftest(rlm.pd, var = "UPSIT")$p.value
  }
  res_new$n.samp <- sum(!(is.na(test$UPSIT)))
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$cluster <- c(1,1,1,2,2,2,2,2)
res$study <- c("PPMI","PDBP","BioFIND","PPMI","PDBP","HBS","SURE","STEADY")
res <- res[res$n.samp>0,]
write.table(res, "./meta_analysis/res/UPSIT.rlm.tsv", quote = F, row.names = F, sep = "\t")

