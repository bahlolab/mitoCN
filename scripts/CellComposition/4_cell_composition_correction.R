rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")

### cell composition data
# load cell composition estimates
cell.comp <- read.table("./meta_analysis/cell_composition/results/CIBERSORTx_PP_CPM.txt", header = T, sep = "\t")
# QC: select cell types with proportion >1%
vars <- c(colnames(cell.comp[,2:23])[colMeans(cell.comp[,2:23]) > 0.01],"Absolute.score..sig.score.")
cells <- cell.comp[,c(vars)]
convertID <- function(x) y <- paste(x[1:2], collapse = "-")
cells$ID <- unlist(lapply(strsplit(cell.comp$Mixture,"[.]"), convertID))
# absolute proportions for cluster 2
abs <- cells[,c(10,1:9)]
colnames(abs)[10] <- "absolute.score"
write.table(abs, "./PPMI/analysis/mtDNA_CN/results/abs.cell.composition.txt", row.names = F, quote = F, sep = "\t")
# relative proportions for cluster 1
# rlt <- cbind(abs$ID,abs[,10], abs[,2:9]/abs[,10])
# colnames(rlt)[1:2] <- c("ID","abs.score")
# write.table(rlt,"./PPMI/analysis/mtDNA_CN/results/relative.cell.composition.txt", row.names = F, quote = F, sep = "\t")

### cell composition correction
# load mtDNA-CN estimates
load(file = "./meta_analysis/mt.all.RData")
mt <- mt[mt$study.cluster == "PPMI.2",c("ID","mt","logmt","cluster")]

# cluster 2
mt.cell <- merge(mt, abs, by = "ID", all.x = T)
mt.cell <- na.omit(mt.cell)
### model selection
# full model
m0 <- lm(logmt ~ ., data = mt.cell[,-c(1,2,4)])
summary(m0)
### model selection
# Step-wise regression
m1 <- step(lm0, direction = "both")
summary(m1) #R2 = 16%

# cell composition correction
mt.adj <- round(residuals(m1),4)
adj.mt <- cbind(mt.cell, mt.adj)
result <- adj.mt[,c("ID","mt","logmt","mt.adj")]
write.table(result, "./PPMI/analysis/mtDNA_CN/results/adj.mt.txt", row.names = F, quote = F, sep = "\t")
