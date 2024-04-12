# make sample list files by cohort, make sure all samples with WGS included
load(file = "/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/mt.all.RData")
table(mt$study)
BioFIND <- data.frame(cbind(mt$ID[mt$study == "BioFIND"],mt$ID[mt$study == "BioFIND"]))
HBS <- data.frame(cbind(mt$ID[mt$study == "HBS"],mt$ID[mt$study == "HBS"]))
LBD <- data.frame(cbind(mt$ID[mt$study == "LBD"],mt$ID[mt$study == "LBD"]))
PDBP <- data.frame(cbind(mt$ID[mt$study == "PDBP"],mt$ID[mt$study == "PDBP"]))
PPMI <- data.frame(cbind(mt$ID[mt$study == "PPMI"],mt$ID[mt$study == "PPMI"]))
STEADY <- data.frame(cbind(mt$ID[mt$study == "STEADY-PD3"],mt$ID[mt$study == "STEADY-PD3"]))
SURE <- data.frame(cbind(mt$ID[mt$study == "SURE-PD"],mt$ID[mt$study == "SURE-PD"]))

colnames(BioFIND) <- colnames(HBS) <- colnames(LBD) <- colnames(PDBP) <-
  colnames(PPMI) <- colnames(STEADY) <- colnames(SURE) <- c("#FID","IID")

write.table(BioFIND, "BioFIND_ID.txt", row.names = F, quote = F)
write.table(HBS, "HBS_ID.txt", row.names = F, quote = F)
write.table(LBD, "LBD_ID.txt", row.names = F, quote = F)
write.table(PDBP, "PDBP_ID.txt", row.names = F, quote = F)
write.table(PPMI, "PPMI_ID.txt", row.names = F, quote = F)
write.table(STEADY, "STEADY_ID.txt", row.names = F, quote = F)
write.table(SURE, "SURE_ID.txt", row.names = F, quote = F)

