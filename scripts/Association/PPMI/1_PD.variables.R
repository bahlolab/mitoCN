rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")

mt <- read.table("./PPMI/analysis/mtDNA_CN/results/adj.mt.txt", header = T, sep = "\t")
IDs <- mt$ID
# load metadata
metadata <- read.csv("./metadata/v3/metadata_20230320.csv")[,c("participant_id","sex","age_at_baseline","diagnosis_at_baseline")]
metadata <- metadata[metadata$participant_id %in% IDs,]
table(metadata$diagnosis_at_baseline)
metadata$diagnosis_at_baseline[metadata$diagnosis_at_baseline %in% c("Prodromal motor PD","Prodromal non-motor PD")] <- "Prodromal PD"
metadata$diagnosis_at_baseline[metadata$diagnosis_at_baseline %in% c("Essential Tremor","Other Neurological Disorder(s)")] <- NA
metadata$diagnosis_at_baseline[metadata$diagnosis_at_baseline == "Idiopathic PD"] <- "PD"
metadata$diagnosis_at_baseline[metadata$diagnosis_at_baseline == "No PD Nor Other Neurological Disorder"] <- "HC"
metadata$diagnosis_at_baseline <- as.factor(metadata$diagnosis_at_baseline)
eigenVectors <- read.table("./meta_analysis/PCA/results/plinkPCA_all.eigenvec")[,c(1,3:7)]
colnames(eigenVectors) <- c("ID",paste0("PC",1:5))
meta <- merge(metadata, eigenVectors, by.x = "participant_id", by.y = "ID", all.x = T)

dat0 <- merge(mt, meta, by.x = "ID", by.y = "participant_id")

### MDS_UPDRS_I
MDS_UPDRS_I <- read.csv("./metadata/v3/clinical/MDS_UPDRS_Part_I.csv")[,c("participant_id","visit_name","mds_updrs_part_i_summary_score")]
MDS_UPDRS_I <- MDS_UPDRS_I[MDS_UPDRS_I$participant_id %in% IDs & 
                                 !is.na(MDS_UPDRS_I$mds_updrs_part_i_summary_score) & 
                                 (MDS_UPDRS_I$visit_name == "M0" | MDS_UPDRS_I$visit_name == "SC"),]
dup.ID <- MDS_UPDRS_I$participant_id[duplicated(MDS_UPDRS_I$participant_id)]
rm <- which((MDS_UPDRS_I$participant_id %in% dup.ID) & MDS_UPDRS_I$visit_name == "SC")
MDS_UPDRS_I <- MDS_UPDRS_I[-rm,]
MDS_UPDRS_I <- MDS_UPDRS_I[!(duplicated(MDS_UPDRS_I$participant_id)),-2]

### MDS_UPDRS_II
MDS_UPDRS_II <- read.csv("./metadata/v3/clinical/MDS_UPDRS_Part_II.csv")[,c("participant_id","visit_name","mds_updrs_part_ii_summary_score")]
MDS_UPDRS_II <- MDS_UPDRS_II[MDS_UPDRS_II$participant_id %in% IDs & 
                             !is.na(MDS_UPDRS_II$mds_updrs_part_ii_summary_score) & 
                             (MDS_UPDRS_II$visit_name == "M0" | MDS_UPDRS_II$visit_name == "SC"),]
dup.ID <- MDS_UPDRS_II$participant_id[duplicated(MDS_UPDRS_II$participant_id)]
rm <- which((MDS_UPDRS_II$participant_id %in% dup.ID) & MDS_UPDRS_II$visit_name == "SC")
MDS_UPDRS_II <- MDS_UPDRS_II[-rm,]
MDS_UPDRS_II <- MDS_UPDRS_II[!(duplicated(MDS_UPDRS_II$participant_id)),-2]

### MDS_UPDRS_III
MDS_UPDRS_III <- read.csv("./metadata/v3/clinical/MDS_UPDRS_Part_III.csv")[,c("participant_id","visit_name","mds_updrs_part_iii_summary_score")]
MDS_UPDRS_III <- MDS_UPDRS_III[MDS_UPDRS_III$participant_id %in% IDs & 
                                 !is.na(MDS_UPDRS_III$mds_updrs_part_iii_summary_score) & 
                                 (MDS_UPDRS_III$visit_name == "M0" | MDS_UPDRS_III$visit_name == "SC"),]
dup.ID <- MDS_UPDRS_III$participant_id[duplicated(MDS_UPDRS_III$participant_id)]
rm <- which((MDS_UPDRS_III$participant_id %in% dup.ID) & MDS_UPDRS_III$visit_name == "SC")
MDS_UPDRS_III <- MDS_UPDRS_III[-rm,]
MDS_UPDRS_III <- MDS_UPDRS_III[!(duplicated(MDS_UPDRS_III$participant_id)), -2]

### MDS_UPDRS_IV
MDS_UPDRS_IV <- read.csv("./metadata/v3/clinical/MDS_UPDRS_Part_IV.csv")[,c("participant_id","visit_name","mds_updrs_part_iv_summary_score")]
MDS_UPDRS_IV <- MDS_UPDRS_IV[MDS_UPDRS_IV$participant_id %in% IDs & 
                                 !is.na(MDS_UPDRS_IV$mds_updrs_part_iv_summary_score) & 
                                 (MDS_UPDRS_IV$visit_name == "M0" | MDS_UPDRS_IV$visit_name == "SC"),]
dup.ID <- MDS_UPDRS_IV$participant_id[duplicated(MDS_UPDRS_IV$participant_id)]
rm <- which((MDS_UPDRS_IV$participant_id %in% dup.ID) & MDS_UPDRS_IV$visit_name == "SC")
MDS_UPDRS_IV <- MDS_UPDRS_IV[-rm,]
MDS_UPDRS_IV <- MDS_UPDRS_IV[!(duplicated(MDS_UPDRS_IV$participant_id)), -2]

### MDS_UPDRS
MDS_UPDRS <- merge(merge(merge(MDS_UPDRS_I, MDS_UPDRS_II, by = "participant_id", all = T), 
                         MDS_UPDRS_III, by = "participant_id", all = T), 
                   MDS_UPDRS_IV, by = "participant_id", all = T)
colnames(MDS_UPDRS) <- c("ID", "MDS_UPDRS_I","MDS_UPDRS_II","MDS_UPDRS_III","MDS_UPDRS_IV")
write.table(MDS_UPDRS, "./PPMI/analysis/mtDNA_CN/results/MDS_UPDRS.txt", row.names = F, sep = "\t", quote = F)

##### MoCA #####
MOCA <- read.csv("./metadata/v3/clinical/MOCA.csv")[,c("participant_id","visit_name","moca_total_score")]
MOCA <- MOCA[MOCA$participant_id %in% IDs & 
                                 !is.na(MOCA$moca_total_score) & 
                                 (MOCA$visit_name == "M0" | MOCA$visit_name == "SC"),]
dup.ID <- MOCA$participant_id[duplicated(MOCA$participant_id)]
rm <- which((MOCA$participant_id %in% dup.ID) & MOCA$visit_name == "SC")
MOCA <- MOCA[-rm,]
MOCA <- MOCA[!(duplicated(MOCA$participant_id)), -2]
colnames(MOCA) <- c("ID","MoCA")

##### RBD #####
RBD <- read.csv("./metadata/v3/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv")[,c("participant_id","visit_name","rbd_summary_score")]
RBD <- RBD[RBD$participant_id %in% IDs & 
               !is.na(RBD$rbd_summary_score) & 
               (RBD$visit_name == "M0" | RBD$visit_name == "SC"),]
dup.ID <- RBD$participant_id[duplicated(RBD$participant_id)]
RBD <- RBD[,-2]
colnames(RBD) <- c("ID","RBD")

##### ESS #####
ESS <- read.csv("./metadata/v3/clinical/Epworth_Sleepiness_Scale.csv")[,c("participant_id","visit_name","ess_summary_score")]
ESS <- ESS[ESS$participant_id %in% IDs & 
             !is.na(ESS$ess_summary_score) & 
             (ESS$visit_name == "M0" | ESS$visit_name == "SC"),]
dup.ID <- ESS$participant_id[duplicated(ESS$participant_id)]
ESS <- ESS[,-2]
colnames(ESS) <- c("ID","ESS")

##### UPSIT #####
UPSIT <- read.csv("./metadata/v3/clinical/UPSIT.csv")[,c("participant_id","visit_name","upsit_total_score")]
UPSIT <- UPSIT[UPSIT$participant_id %in% IDs & 
             !is.na(UPSIT$upsit_total_score) & 
             (UPSIT$visit_name == "M0" | UPSIT$visit_name == "SC"),]
dup.ID <- UPSIT$participant_id[duplicated(UPSIT$participant_id)]
UPSIT <- UPSIT[,-2]
colnames(UPSIT) <- c("ID","UPSIT")

##### ADL #####
ADL <- read.csv("./metadata/v3/clinical/Modified_Schwab___England_ADL.csv")[,c("participant_id","visit_name","mod_schwab_england_pct_adl_score")]
ADL <- ADL[ADL$participant_id %in% IDs & 
                 !is.na(ADL$mod_schwab_england_pct_adl_score) & 
                 (ADL$visit_name == "M0" | ADL$visit_name == "SC"),]
dup.ID <- ADL$participant_id[duplicated(ADL$participant_id)]
rm <- which((ADL$participant_id %in% dup.ID) & ADL$visit_name == "SC")
ADL <- ADL[-rm,]
ADL <- ADL[!(duplicated(ADL$participant_id)), -2]
colnames(ADL) <- c("ID","ADL")

##### MRI #####
MRI <- read.csv("./metadata/v3/clinical/MRI.csv")[,c("participant_id","visit_name","mri_results")]
MRI <- MRI[MRI$participant_id %in% IDs & 
             !is.na(MRI$mri_results) & 
             (MRI$visit_name == "M0" | MRI$visit_name == "SC"),]
dup.ID <- MRI$participant_id[duplicated(MRI$participant_id)]
MRI <- MRI[,-2]
colnames(MRI) <- c("ID","MRI")
MRI$MRI[MRI$MRI == ""] <- NA

##### CSF #####
CSF <- read.csv("./metadata/v3/clinical/Biospecimen_analyses_CSF_abeta_tau_ptau.csv")[,c("participant_id","visit_name","test_name","test_value")]
CSF <- CSF[CSF$participant_id %in% IDs & 
             (CSF$visit_name == "M0" | CSF$visit_name == "SC"),]
Abeta <- CSF[CSF$test_name == "Abeta",c("participant_id","visit_name","test_value")]
colnames(Abeta)[3] <- "Abeta"
p.tau <- CSF[CSF$test_name == "p-Tau",c("participant_id","visit_name","test_value")]
colnames(p.tau)[3] <- "p_Tau"
tau <- CSF[CSF$test_name == "Tau",c("participant_id","visit_name","test_value")]
colnames(tau)[3] <- "Tau"
CSF <- merge(Abeta, p.tau, by = c("participant_id","visit_name"), all = T)
CSF <- merge(CSF, tau, by = c("participant_id","visit_name"), all = T)[,-2]
dup.ID <- CSF$participant_id[duplicated(CSF$participant_id)]

##### DaTSCAN #####
DaTSCAN <- read.csv("./metadata/v3/clinical/DaTSCAN_visual_interpretation.csv")[,c("participant_id","visit_name","scan_months_after_baseline","datscan_visual_interpretation")]
DaTSCAN <- DaTSCAN[DaTSCAN$participant_id %in% IDs & DaTSCAN$visit_name == "LOG",]
dup.ID <- DaTSCAN$participant_id[duplicated(DaTSCAN$participant_id)]
DaTSCAN <-DaTSCAN[,c(1,4)]
colnames(DaTSCAN) <- c("ID","DaTSCAN")

#### merge all PD variables
dat1 <- merge(dat0, MDS_UPDRS, by = "ID", all = T)
dat2 <- merge(dat1, MOCA, by = "ID", all = T)
dat3 <- merge(dat2, RBD, by = "ID", all = T)
dat4 <- merge(dat3, ESS, by = "ID", all = T)
dat5 <- merge(dat4, UPSIT, by = "ID", all = T)
dat6 <- merge(dat5, ADL, by = "ID", all = T)
dat7 <- merge(dat6, MRI, by = "ID", all = T)
dat8 <- merge(dat7, CSF, by.x = "ID", by.y = "participant_id", all = T)
data <- merge(dat8, DaTSCAN, by = "ID", all = T)
data$sex <- as.factor(data$sex)
data$MRI <- factor(data$MRI, levels = c("Normal","Abnormal, not clinically significant","Abnormal, clinically significant"))
data$DaTSCAN <- as.factor(data$DaTSCAN)
summary(data)
save(data, file = "./PPMI/analysis/mtDNA_CN/results/data.RData")
