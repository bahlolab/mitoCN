rm(list=ls())
library(magrittr)
library(sfsmisc)
library(MASS)
library(ggplot2)
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")

# load mtDNA-CN results
load(file = "./PPMI/analysis/mtDNA_CN/results/data.RData")
colnames(data)[6:7] <- c("age","diagnosis")
summary(data)

### raw mtDNA-CN
test <- data[,c(3,5,6,8:12,7)]
res <- summary(m1 <- rlm(formula = logmt ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(m1, var = 9)$p.value

for(i in 13:25){
  test <- data[,c(3,5,6,8:12,i)]
  res_new <- summary(m1 <- rlm(formula = logmt ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(m1, var = 9)$p.value
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$PD.vars <- colnames(data)[c(7,13:25)]

write.table(res, "./PPMI/analysis/mtDNA_CN/results/PD.raw.mt.tsv", quote = F, row.names = F, sep = "\t")

### adjusted mtDNA-CN
test <- data[,c(4,5,6,8:12,7)]
res <- summary(m1 <- rlm(formula = mt.adj ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
res$p.value <- f.robftest(m1, var = 9)$p.value

for(i in 13:25){
  test <- data[,c(4,5,6,8:12,i)]
  res_new <- summary(m1 <- rlm(formula = mt.adj ~ ., data = test))$coefficients[9,] %>% unlist %>% matrix(nrow=1) %>% as.data.frame
  res_new$p.value <- f.robftest(m1, var = 9)$p.value
  res <- rbind(res, res_new)
}
colnames(res)[1:4] <- c("Estimate","Std. Error","t value","Pr(>|t|)")
res$PD.vars <- colnames(data)[c(7,13:25)]
res$n.samp <- apply(data[,c(7,13:25)],2, function(x) sum(!(is.na(x))))
write.table(res, "./PPMI/analysis/mtDNA_CN/results/PD.adj.mt.tsv", quote = F, row.names = F, sep = "\t")

