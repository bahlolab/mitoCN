rm(list=ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")

metadata <- read.csv("./metadata/v3/metadata_20230320.csv")
IDs <- metadata$participant_id[!(metadata$study %in% c("LBD","LCC"))]
meta <- metadata[!(metadata$study %in% c("LBD","LCC")) & metadata$wgs_jnt_sample==1,c(1,25,28,32:46)]
load(file = "./meta_analysis/mt.all.RData")
cluster <- mt[,c("ID","study.cluster")]
meta <- merge(meta, cluster, by.x = "participant_id", by.y = "ID")
data.clinic <- aggregate(. ~ study.cluster, meta[,-1], sum)
write.table(data.clinic, "./meta_analysis/res/clinical.summary.tsv", quote = F, row.names = F, sep = "\t")
# based on the summary, select MDS_UPDRS_III, ADL, MOCA for association tests

### MDS_UPDRS_III
MDS_UPDRS_III <- read.csv("./metadata/v3/clinical/MDS_UPDRS_Part_III.csv")[,c("participant_id","visit_name","mds_updrs_part_iii_summary_score")]
MDS_UPDRS_III <- MDS_UPDRS_III[MDS_UPDRS_III$participant_id %in% IDs & 
                                 !is.na(MDS_UPDRS_III$mds_updrs_part_iii_summary_score) & 
                                 (MDS_UPDRS_III$visit_name == "M0" | MDS_UPDRS_III$visit_name == "SC"),]
dup.ID <- MDS_UPDRS_III$participant_id[duplicated(MDS_UPDRS_III$participant_id)]
rm <- which((MDS_UPDRS_III$participant_id %in% dup.ID) & MDS_UPDRS_III$visit_name == "SC")
MDS_UPDRS_III <- MDS_UPDRS_III[-rm,]
MDS_UPDRS_III <- MDS_UPDRS_III[!(duplicated(MDS_UPDRS_III$participant_id)),]

### ADL
ADL <- read.csv("./metadata/v3/clinical/Modified_Schwab___England_ADL.csv")[,c("participant_id","visit_name","mod_schwab_england_pct_adl_score")]
ADL <- ADL[ADL$participant_id %in% IDs & !is.na(ADL$mod_schwab_england_pct_adl_score) & (ADL$visit_name == "M0" | ADL$visit_name == "SC"),]
dup.ID <- ADL$participant_id[duplicated(ADL$participant_id)]
rm <- which((ADL$participant_id %in% dup.ID) & ADL$visit_name == "SC")
ADL <- ADL[-rm,]
ADL <- ADL[!(duplicated(ADL$participant_id)),]

### MOCA
MOCA <- read.csv("./metadata/v3/clinical/MOCA.csv")[,c("participant_id","visit_name","code_education_12years_complete","moca_total_score")]
MOCA <- MOCA[(MOCA$participant_id %in% IDs) & (!is.na(MOCA$moca_total_score)) & (MOCA$visit_name == "M0" | MOCA$visit_name == "SC"),]
dup.ID <- MOCA$participant_id[duplicated(MOCA $participant_id)]
rm <- which((MOCA$participant_id %in% dup.ID) & MOCA$visit_name == "SC")
MOCA  <- MOCA[-rm,]
MOCA  <- MOCA[!(duplicated(MOCA$participant_id)),]

### RBD
RBD <- read.csv("./metadata/v3/clinical/REM_Sleep_Behavior_Disorder_Questionnaire_Stiasny_Kolster.csv")[,c("participant_id","visit_name","rbd_summary_score")]
RBD <- RBD[RBD$participant_id %in% IDs & !is.na(RBD$rbd_summary_score) & (RBD$visit_name == "M0" | RBD$visit_name == "SC"),]
dup.ID <- RBD$participant_id[duplicated(RBD$participant_id)]
RBD  <- RBD[!(duplicated(RBD$participant_id)),]

### UPSIT
UPSIT <- read.csv("./metadata/v3/clinical/UPSIT.csv")[,c("participant_id","visit_name","upsit_total_score")]
UPSIT <- UPSIT[UPSIT$participant_id %in% IDs & !is.na(UPSIT$upsit_total_score) & (UPSIT$visit_name == "M0" | UPSIT$visit_name == "SC"),]
dup.ID <- UPSIT$participant_id[duplicated(UPSIT$participant_id)]
UPSIT <- UPSIT[!(duplicated(UPSIT$participant_id)),]

metadata <- read.csv("./metadata/v3/metadata_20230320.csv")[,c("participant_id","sex","age_at_baseline","case_control_other_at_baseline")]
data0 <- metadata[metadata$participant_id %in% IDs,]
data1 <- merge(data0, MDS_UPDRS_III[,-2], by = "participant_id", all.x = T)
data2 <- merge(data1, ADL[,-2], by = "participant_id", all.x = T)
data3 <- merge(data2, MOCA[,-2], by = "participant_id", all.x = T)
data4 <- merge(data3, RBD[,-2], by = "participant_id", all.x = T)
data5 <- merge(data4, UPSIT[,-2], by = "participant_id", all.x = T)
colnames(data5) <- c("ID","sex","age","diagnosis","MDS_UPDRS_III","ADL","education","MOCA","RBD","UPSIT")
write.table(data5, "./meta_analysis/res/PD.vars.tsv", quote = F, row.names = F, sep = "\t")

