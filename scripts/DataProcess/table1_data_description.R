#######################################################################################################
#################################### AMP-PD data description ##########################################
#######################################################################################################
rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")


### load metadata
metadata <- read.csv("./metadata/v3/metadata_20230525.csv")
metadata$wgs_jnt_sample <- as.factor(metadata$wgs_jnt_sampl)
metadata$study <- as.factor(metadata$study)
summary(metadata)

# participants with WGS data
# exclude LBD & LCC
dat <- metadata[metadata$wgs_jnt_sample==1 & !(metadata$study %in% c("LBD","LCC")) ,-1]

# combine with cluster info
load(file = "./meta_analysis/mt.all.RData")
mt <- mt[,c("ID","study.cluster","cluster","mt","logmt")]
mt.data <- merge(mt, dat, by.x = "ID", by.y = "participant_id", )


### generate table
# sample size
mt.data$study.cluster <- factor(mt.data$study.cluster, exclude = c("LBD.2","LBD.3","LBD.4"))
samp.size <- table(mt.data$study.cluster)

# sex
mt.data$sex <- as.factor(mt.data$sex)
tab.sex <- table(mt.data$study.cluster, mt.data$sex)
tab1 <- cbind(samp.size, tab.sex)

# case/control/other
mt.data$case_control_other_at_baseline <- as.factor(mt.data$case_control_other_at_baseline)
tab.dx <- table(mt.data$study.cluster,mt.data$case_control_other_at_baseline)
tab2 <- cbind(tab1, tab.dx)

# has known PD mutations
mt.data$has_known_PD_Mutation_in_WGS <- as.factor(mt.data$has_known_PD_Mutation_in_WGS)
tab.mut <- table(mt.data$study.cluster, mt.data$has_known_PD_Mutation_in_WGS)
tab3 <- cbind(tab2, tab.mut)

# RNAseq data
rna <- read.csv("./metadata/v3/rnaseq_WB-RWTS_sample_inventory.csv")
ran.ids <- unique(rna$participant_id)
mt.data$rnaseq <- 0
mt.data$rnaseq[mt.data$ID %in% ran.ids] <- 1
mt.data$rnaseq <- as.factor(mt.data$rnaseq)
tab.rna <- table(mt.data$study.cluster, mt.data$rnaseq)
tab4 <- cbind(tab3, tab.rna)

# age: mean, SD
age.mean <- aggregate(mt.data$age_at_baseline, list(mt.data$study.cluste), FUN=mean)
age.sd <- aggregate(mt.data$age_at_baseline, list(mt.data$study.cluste), FUN=sd)
age <- merge(age.mean, age.sd, by = "Group.1")
colnames(age) <- c("study","age.mean","age.sd")
rownames(age) <- age$study
age <- age[,-1]
tab5 <- cbind(tab4, age)
mean(mt.data$age_at_baseline); sd(mt.data$age_at_baseline)
colnames(tab5)[7:10] <- c("mut.no","mut.yes","rna.0","rna.1")
write.table(tab5, "./meta_analysis/description.table.tsv", quote = F, sep = "\t")
