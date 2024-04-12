rm(list=ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")
library(ggridges)
library(ggplot2)
# load mtDNA-CN results
PPMI <- read.table("./PPMI/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt","cluster")]
PDBP <- read.table("./PDBP/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt","cluster")]
HBS <- read.table("./HBS/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt")]
BioFIND <- read.table("./BioFIND/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt")]
Sure <- read.table("./SURE-PD/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt")]
Steady<- read.table("./STEADY-PD3/analysis/mtDNA_CN/mtDNA_CN.txt", header = T, sep = "\t")[,c("ID","mt")]

BioFIND$cluster <- 1
HBS$cluster <- 2
Sure$cluster <- 2
Steady$cluster <- 2

# add study column
PPMI$study <- "PPMI"
PDBP$study <- "PDBP"
HBS$study <- "HBS"
BioFIND$study <- "BioFIND"
Sure$study <- "SURE-PD"
Steady$study <- "STEADY-PD3"

# merge all data
mt <- rbind(PPMI,PDBP,HBS,BioFIND,Sure,Steady)
mt$study <- as.factor(mt$study)
mt$logmt <- log(mt$mt)
mt$cluster <- as.factor(mt$cluster)

ridge <- ggplot(mt, aes(x = logmt, y = study, fill = study)) +
  geom_density_ridges() +
  theme(legend.position = "none") + 
  xlab("log(mtDNA-CN)") +
  ylab("Cohort")

cluster <- ggplot(mt, aes(x = logmt, y = study, fill = cluster)) +
  geom_density_ridges() +
  xlab("log(mtDNA-CN)") +
  ylab("Cohort") +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75)

p_both <- ggarrange(ridge, cluster, nrow = 2, labels = c("A", "B"))

ggsave(
  filename = "./meta_analysis/img/cluster.png",
  plot = p_both,
  width = 16,
  height = 14,
  device = "png"
)

