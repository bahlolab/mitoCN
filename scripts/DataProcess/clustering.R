rm(list=ls())
library(ggplot2)
library(mclust)

study <- "PPMI"
setwd(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/",study,"/analysis/mtDNA_CN/"))

# load mtDNA-CN results
data <- read.table("mtDNA_CN.txt", header = T, sep = "\t")[,1:10]
data$logmt <- log(data$mt)
ggplot(data, aes(x=logmt)) +  geom_histogram()

# mixture model
df <- data$logmt
BIC <- mclustBIC(df); plot(BIC); summary(BIC)
ICL <- mclustICL(df); plot(ICL); summary(ICL)
gmm2 <- Mclust(df, 2) 

gg2 <- cbind(data, cluster = gmm2$classification)
gg2$cluster <-  as.factor(gg2$cluster)
summary(gg2)

# histogram 
ggplot(gg2, aes(x=logmt, color=cluster, fill=cluster)) +
  geom_histogram(position="identity", alpha=0.5)

write.table(gg2, "mtDNA_CN.txt", row.names = F, quote = F, sep = "\t")

