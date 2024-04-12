#######################################################################################################
####################################### PCA using Plink 2 #############################################
#######################################################################################################
rm(list = ls())
library(ggplot2)
cohort <- "LBD"
work_dir <- "/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/PCA/"
# read in result files
eigenValues <- read.table(paste0(work_dir,"results/plinkPCA_",cohort,".eigenval"))
eigenVectors <- read.table(paste0(work_dir,"results/plinkPCA_",cohort,".eigenvec"))
## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)
# get population information
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/",
            cohort,"/analysis/metadata/",cohort,".v3.with.control.info.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/SURE-PD/analysis/metadata2/Sure.v3.with.control.info.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/STEADY-PD3/analysis/metadata2/Steady.v3.with.control.info.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/LBD/analysis/metadata/LBD.v3.RData"))

meta <- LBD.mt
table(meta$ethnicity, meta$race)
colnames(eigenVectors) <- c("FID","IID",paste0("PC",1:10))
PCs <- merge(meta[,c("participant_id","sex","race","ethnicity")], eigenVectors, by.x = "participant_id", by.y = "IID")

# PCA plot
PC12 <- ggplot(data = PCs) +
  geom_point(mapping = aes(x = PC1, y = PC2, colour = race, shape = ethnicity), size = 3, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = paste0("PCA of ",cohort),
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
       colour = "Ancestry", shape = "Ethnicity") +
  theme_minimal()
ggsave(filename = paste0(work_dir,"img/",cohort,"_PC12.png"),
      plot = PC12, width = 10, height = 5,)

PC23 <-ggplot(data = PCs) +
  geom_point(mapping = aes(x = PC2, y = PC3, colour = race, shape = ethnicity), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = paste0("PCA of ",cohort),
       x = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
       y = paste0("Principal component 3 (",eigen_percent[3,1]," %)"),
       colour = "Ancestry", shape = "Ethnicity") +
  theme_minimal()
ggsave(filename = paste0(work_dir,"img/",cohort,"_PC23.png"),
       plot = PC23, width = 10, height = 5,)

#######################################################################################################
####################################### PCA for all cohort ############################################
#######################################################################################################
rm(list = ls())
library(ggplot2)
work_dir <- "/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/PCA/"
# read in result files
eigenValues <- read.table(paste0(work_dir,"results/plinkPCA_all.eigenval"))
eigenVectors <- read.table(paste0(work_dir,"results/plinkPCA_all.eigenvec"))
## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)
# get population information
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/BioFIND/analysis/metadata/BioFIND.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/PPMI/analysis/metadata/ppmi.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/PDBP/analysis/metadata/PDBP.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/HBS/analysis/metadata/HBS.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/SURE-PD/analysis/metadata/Sure.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/STEADY-PD3/analysis/metadata/Steady.v3.RData"))
load(paste0("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/LBD/analysis/metadata/LBD.v3.RData"))

BioFIND <- BioFIND[,c("participant_id","race","ethnicity","study")]
PPMI <- ppmi[,c("participant_id","race","ethnicity","cluster")]
PDBP <- PDBP[,c("participant_id","race","ethnicity","cluster")]
HBS <- HBS[,c("participant_id","race","ethnicity","study")]
Sure <- Sure.mt[,c("participant_id","race","ethnicity","study","cluster")]
Steady <- Steady[,c("participant_id","race","ethnicity","study")]
LBD <- LBD.mt[,c("participant_id","race","cluster")]
LBD$ethnicity <- "Unknown"

BioFIND$cluster <- 1
PPMI$study <- "PPMI"
PDBP$study <- "PDBP"
HBS$cluster <- 2
Sure$study <- "Sure"
Steady$cluster <- 2
LBD$study <- "LBD"

meta <- rbind(BioFIND, PPMI, PDBP, HBS, Sure, Steady, LBD)
meta$study <- as.factor(meta$study)
meta$cluster <- as.factor(meta$cluster)
table(meta$ethnicity, meta$race)
colnames(eigenVectors) <- c("FID","IID",paste0("PC",1:10))
PCs <- merge(meta[,c("participant_id","race","ethnicity","study","cluster")], eigenVectors, by.x = "participant_id", by.y = "IID")

# PCA plot
PC12 <- ggplot(data = PCs[PCs,]) +
  geom_point(mapping = aes(x = PC1, y = PC2, colour = race, shape = ethnicity), size = 3, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = paste0("PCA of all"),
       x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
       colour = "Ancestry", shape = "Ethnicity") +
  theme_minimal()
ggsave(filename = paste0(work_dir,"img/all_PC12.png"),
       plot = PC12, width = 10, height = 5,)

PC23 <-ggplot(data = PCs) +
  geom_point(mapping = aes(x = PC2, y = PC3, colour = race, shape = ethnicity), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = paste0("PCA of all"),
       x = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
       y = paste0("Principal component 3 (",eigen_percent[3,1]," %)"),
       colour = "Ancestry", shape = "Ethnicity") +
  theme_minimal()
ggsave(filename = paste0(work_dir,"img/all_PC23.png"),
       plot = PC23, width = 10, height = 5,)


ggplot(data = PCs[PCs$study=="PDBP",]) +
  geom_point(mapping = aes(x = PC1, y = PC2, colour = race, shape = ethnicity), size = 3, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = paste0("PCA of all"))+
  theme_minimal()
