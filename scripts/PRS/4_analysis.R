##### mtCN PRS validation #####
rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/")
library(ggplot2)
library(ggstatsplot) # for boxplot with stats
library(ggpubr) #for multiple plots

# combine two PRSs
adj_mtCN_PRS <- read.table("./PRS/results/mtCN/adj_mtcn_PRS.sscore", header = T)[,c(2,5)]
raw_mtCN_PRS <- read.table("./PRS/results/mtCN/raw_mtCN_PRS.sscore", header = T)[,c(2,5)]
PRS <- merge(adj_mtCN_PRS, raw_mtCN_PRS, by = "IID")
colnames(PRS)[2:3] <- c("adj_PRS","raw_PRS")

# load metadata
metadata <- read.csv("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/metadata/v3/metadata_20230320.csv")
# define true healthy control by removing the known mutation carriers and individuals with family history
HC <- metadata$participant_id[metadata$case_control_other_latest=="Control" & 
                                (metadata$has_known_PD_Mutation_in_WGS != "Yes" | is.na(metadata$has_known_PD_Mutation_in_WGS)) & 
                                (metadata$biological_father_with_pd != "Yes" | is.na(metadata$biological_father_with_pd)) & 
                                (metadata$biological_mother_with_pd != "Yes" | is.na(metadata$biological_mother_with_pd))]
PD <- metadata$participant_id[metadata$case_control_other_latest=="Case"]
Control <- metadata$participant_id[metadata$case_control_other_latest=="Control"]
LBD <- metadata$participant_id[metadata$diagnosis_latest == "LBD"]
PRS$DX <- NA
PRS$DX[PRS$IID %in% Control] <- "Control"
PRS$DX[PRS$IID %in% HC] <- "HC"
PRS$DX[PRS$IID %in% PD] <- "PD"
PRS$DX[PRS$IID %in% LBD] <- "LBD"

PRS_QC <- PRS[!(is.na(PRS$DX)),]

# re-scale PRS
adj_mean <- mean(PRS_QC$adj_PRS[PRS_QC$DX=="HC"])
adj_sd <- sd(PRS_QC$adj_PRS[PRS_QC$DX=="HC"])
PRS_QC$adj_PRS_scaled <- (PRS_QC$adj_PRS-adj_mean)/adj_sd

raw_mean <- mean(PRS_QC$raw_PRS[PRS_QC$DX=="HC"])
raw_sd <- sd(PRS_QC$raw_PRS[PRS_QC$DX=="HC"])
PRS_QC$raw_PRS_scaled <- (PRS_QC$raw_PRS-raw_mean)/raw_sd

load(file = "./mt.all.RData")
mtSwirl <- read.table("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/mtSwirl/mtcn/results/mtSwirl_mtcn.txt", header = T)[,3:4]
metadata <- read.csv("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/metadata/v3/metadata_20230320.csv")[,c(1,5,19,20,47,51)]

mt2 <- merge(mt, mtSwirl, by = "ID")
mt.meta <- merge(mt2, metadata, by.x = "ID", by.y = "participant_id")
mt_PRS <- merge(PRS_QC, mt.meta, by.x = "IID", by.y = "ID")

### PRS density plot by diagnosis
raw_ds <- ggplot(data=mt_PRS[mt_PRS$DX %in% c("HC","PD"),], aes(x=raw_PRS_scaled, group=DX, fill=DX)) +
  geom_density(adjust=1.5, alpha=.4) + xlab("scaled raw mtDNA-CN PRS")
adj_ds <- ggplot(data=mt_PRS[mt_PRS$DX %in% c("HC","PD"),], aes(x=adj_PRS_scaled, group=DX, fill=DX)) +
  geom_density(adjust=1.5, alpha=.4)+ xlab("scaled adjusted mtDNA-CN PRS")
ds_both <- ggarrange(raw_ds, adj_ds, ncol = 2, labels = c("A", "B"))
ggsave("./PRS/img/density.png", plot = ds_both, width = 12, height = 6)

### t test
t.test(raw_PRS_scaled ~ DX, data = mt_PRS[mt_PRS$DX %in% c("HC","PD"),])
t.test(adj_PRS_scaled ~ DX, data = mt_PRS[mt_PRS$DX %in% c("HC","PD"),])

### PRS boxplot by diagnosis
raw_box <- ggbetweenstats(
  data  = mt_PRS[mt_PRS$DX %in% c("HC","PD"),],
  x     = DX,
  y     = raw_PRS_scaled,
  title = "Distribution of scaled raw mtcn PRS across diagnosis"
)

adj_box <- ggbetweenstats(
  data  = mt_PRS[mt_PRS$DX %in% c("HC","PD"),],
  x     = DX,
  y     = adj_PRS_scaled,
  title = "Distribution of scaled adjusted mtcn PRS across diagnosis"
)
ggsave("./PRS/img/boxplot_raw.png", plot = raw_box, width = 8, height = 5)
ggsave("./PRS/img/boxplot_adj.png", plot = adj_box, width = 8, height = 5)


