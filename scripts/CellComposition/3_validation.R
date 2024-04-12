rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/cell_composition/")
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr) # for adding correlation coef 
library(ggstatsplot)

### function: data pre-processing
process <- function(data, percent = TRUE){
  data$B.cells <- rowSums(data[,startsWith(colnames(data),"B.cells")])
  data$NK.cells <- rowSums(data[,startsWith(colnames(data),"NK.cells")])
  data$T.cells <- rowSums(data[,startsWith(colnames(data),"T.cells")])
  data$Lymphocytes <- rowSums(data[,c("B.cells","T.cells","NK.cells")])
  data$T.cells.CD4 <- rowSums(data[,startsWith(colnames(data),"T.cells.CD4")])
  data1 <- data[,c("Neutrophils","Lymphocytes","Monocytes","T.cells",
                   "T.cells.CD8","T.cells.CD4","B.cells","NK.cells")]
  if(percent){
    data2 <- data1/data$Absolute.score..sig.score.
    return(data2)
  }else{
    return(data1)
  }
}

###### whole blood data from CIBERSORTx
# load true cell type
true <- read.table("./data/validation/Fig2b_ground_truth_whole_blood.txt", header = T, sep = "\t")
est <- read.table("./data/validation/CIBERSORTx_Job1_Adjusted.txt", header = T, sep = "\t")
data <- process(est)
tab1 <- rbind(unlist(lapply(data, mean, 1)), unlist(lapply(data, sd, 1)))
ggcorrmat(data=true,cor.vars = Neutrophils:NK.cells)
ggcorrmat(data)
ggcorrmat(process(est, percent = F))

### boxplot of cell types
data1 <- data %>% gather(key="CellType", value="proportion") 
data1$CellType <- factor(data1$CellType, 
                             levels=c("Neutrophils","Lymphocytes","Monocytes","T.cells",
                                      "T.cells.CD8","T.cells.CD4","B.cells","NK.cells"))
box_type <- ggplot(data1, aes(x=CellType, y=proportion, fill=CellType)) +
  geom_boxplot() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=0.6))

### correlation across cell type
true1 <- true[,-1] %>% gather(key="CellType", value="proportion") 
cmp <- cbind(true1, data1[,2])
colnames(cmp)[2:3] <- c("true","estimate")
p <- ggplot(data = cmp, aes(x = true, y = estimate)) + geom_point() +
  stat_cor(method = "pearson", label.x = 0.15, label.y = 0.05) + xlab("true proportion")
cor_type <- p + facet_wrap(~CellType)

###### whole blood data from AMP PD
BF <- read.table("./results/CIBERSORTx_BF_CPM.txt", header = T, sep = "\t")
PD <- read.table("./results/CIBERSORTx_PD_CPM.txt", header = T, sep = "\t")
PP <- read.table("./results/CIBERSORTx_PP_CPM.txt", header = T, sep = "\t")
res <- rbind(BF, PD, PP)
convertID <- function(x) y <- paste(x[1:2], collapse = "-")
res$ID <- unlist(lapply(strsplit(res$Mixture,"[.]"), convertID))

meta <- read.csv("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/metadata/v3/amp_pd_case_control.csv")
ID.cntl <- meta$participant_id[meta$case_control_other_at_baseline == "Control"]

amp <- process(res[res$ID %in% ID.cntl,])
tab2 <- rbind(unlist(lapply(amp, mean, 1)), unlist(lapply(amp, sd, 1)))
tab <- rbind(tab1, tab2)

### boxplot of cell types
amp1 <- amp %>% gather(key="CellType", value="proportion") 
amp1$CellType <- factor(amp1$CellType, 
                         levels=c("Neutrophils","Lymphocytes","Monocytes","T.cells",
                                  "T.cells.CD8","T.cells.CD4","B.cells","NK.cells"))
amp_box_type <- ggplot(amp1, aes(x=CellType, y=proportion, fill=CellType)) +
  geom_boxplot() + 
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=0.6))

p1 <- ggarrange(box_type, amp_box_type, nrow = 2, labels = c("B", "C"))
fig <- ggarrange(cor_type, p1, ncol = 2, labels = c("A",""), widths = c(1.5,1))
ggsave(filename = "./results/fig.cibersortx.png", plot = fig, width = 12, height = 6)

