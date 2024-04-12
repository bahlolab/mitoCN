rm(list = ls())
setwd("~/OneDrive - wehi.edu.au/Projects/mtDNA-CN/MJFF_mtDNA-CN/analysis/CellComposition/")

### overlap genes: 75.5%
LM22 <- read.table("./data/LM22.txt", header = T, sep = "\t")
data <- read.table("./data/cpm.txt", header = T, sep = "\t")
sum(LM22$Gene.symbol %in% data$Gene_symbol)/nrow(LM22)
samplenames <- colnames(data)

### BioFIND
BF <- samplenames[startsWith(samplenames, "BF")]
data.BF <- cbind(data$Gene_symbol, data[,colnames(data) %in% BF])
write.table(data.BF, "BF.CPM.txt", row.names = F, quote = F, sep = "\t")

### PPMI
PPMI <- samplenames[startsWith(samplenames, "PP")]
data.PPMI <- cbind(data$Gene_symbol, data[,colnames(data) %in% PPMI])
write.table(data.PPMI, "PPMI.CPM.txt", row.names = F, quote = F, sep = "\t")
PPMI.BL <- PPMI[endsWith(PPMI, "BLM0T1")]
data.PPMI.BL <- cbind(data$Gene_symbol, data[,colnames(data) %in% PPMI.BL])
write.table(data.PPMI.BL, "PPMI.BL.CPM.txt", row.names = F, quote = F, sep = "\t")

PPMI.M6 <- PPMI[endsWith(PPMI, "SVM6T1")]
data.PPMI.M6 <- cbind(data.PPMI$data.Gene_symbol, data.PPMI[,colnames(data.PPMI) %in% PPMI.M6])
write.table(data.PPMI.M6, "PPMI.M6.CPM.txt", row.names = F, quote = F, sep = "\t")
PPMI.M12 <- PPMI[endsWith(PPMI, "SVM12T1")]
data.PPMI.M12 <- cbind(data.PPMI$data.Gene_symbol, data.PPMI[,colnames(data.PPMI) %in% PPMI.M12])
write.table(data.PPMI.M12, "PPMI.M12.CPM.txt", row.names = F, quote = F, sep = "\t")
PPMI.M24 <- PPMI[endsWith(PPMI, "SVM24T1")]
data.PPMI.M24 <- cbind(data.PPMI$data.Gene_symbol, data.PPMI[,colnames(data.PPMI) %in% PPMI.M24])
write.table(data.PPMI.M24, "PPMI.M24.CPM.txt", row.names = F, quote = F, sep = "\t")
PPMI.M36 <- PPMI[endsWith(PPMI, "SVM36T1")]
data.PPMI.M36 <- cbind(data.PPMI$data.Gene_symbol, data.PPMI[,colnames(data.PPMI) %in% PPMI.M36])
write.table(data.PPMI.M36, "PPMI.M36.CPM.txt", row.names = F, quote = F, sep = "\t")

### PPMI
PDBP <- samplenames[startsWith(samplenames, "PD")]
data.PDBP <- cbind(data$Gene_symbol, data[,colnames(data) %in% PDBP])
write.table(data.PDBP, "PDBP.CPM.txt", row.names = F, quote = F, sep = "\t")
PDBP.BL <- PDBP[endsWith(PDBP, "BLM0T1")]
data.PDBP.BL <- cbind(data$Gene_symbol, data[,colnames(data) %in% PDBP.BL])
write.table(data.PDBP.BL, "PDBP.BL.CPM.txt", row.names = F, quote = F, sep = "\t")

PDBP.M6 <- PDBP[endsWith(PDBP, "SVM6T1")]
data.PDBP.M6 <- cbind(data$data.Gene_symbol, data[,colnames(data) %in% PDBP.M6])
write.table(data.PDBP.M6, "PDBP.M6.CPM.txt", row.names = F, quote = F, sep = "\t")
PDBP.M12 <- PDBP[endsWith(PDBP, "SVM12T1")]
data.PDBP.M12 <- cbind(data$data.Gene_symbol, data[,colnames(data) %in% PDBP.M12])
write.table(data.PDBP.M12, "PDBP.M12.CPM.txt", row.names = F, quote = F, sep = "\t")
PDBP.M18 <- PDBP[endsWith(PDBP, "SVM18T1")]
data.PDBP.M18 <- cbind(data$data.Gene_symbol, data[,colnames(data) %in% PDBP.M18])
write.table(data.PDBP.M18, "PDBP.M18.CPM.txt", row.names = F, quote = F, sep = "\t")
PDBP.M24 <- PDBP[endsWith(PDBP, "SVM24T1")]
data.PDBP.M24 <- cbind(data$data.Gene_symbol, data[,colnames(data) %in% PDBP.M24])
write.table(data.PDBP.M24, "PDBP.M24.CPM.txt", row.names = F, quote = F, sep = "\t")
