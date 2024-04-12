##### mtCN PRS validation #####
rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/")
library(ggplot2)
library(ggpubr) #for multiple plots
library(RNOmni)
library(plyr)
library(dplyr)
library(ggstatsplot) # for boxplot with stats

adj_mtCN_PRS <- read.table("./PRS/results/mtCN/adj_mtcn_PRS.sscore", header = T)[,c(2,5)]
raw_mtCN_PRS <- read.table("./PRS/results/mtCN/raw_mtCN_PRS.sscore", header = T)[,c(2,5)]
PRS <- merge(adj_mtCN_PRS, raw_mtCN_PRS, by = "IID")
colnames(PRS)[2:3] <- c("adj_PRS","raw_PRS")

## correlation of adj & raw PRS
comp <- ggplot(PRS, aes(x=adj_PRS, y=raw_PRS)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + xlim(-0.007,0.005) + ylim(-0.007,0.005)
comp_p <- comp + stat_cor(label.x = -0.005, r.digits = 2, p.digits = 1)
ggsave(filename = "./PRS/img/correlation2PRS.png", plot = comp_p, width = 10, height = 8)

df <- data.frame(
  PRS_type=factor(rep(c("raw", "adj"), each=nrow(PRS))),
  PRS_value=c(PRS$raw_PRS, PRS$adj_PRS))
head(df)
mu <- ddply(df, "PRS_type", summarise, grp.mean=mean(PRS_value))
p<-ggplot(df, aes(x=PRS_value, color=PRS_type)) +
  geom_density()+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=PRS_type),
             linetype="dashed")
ggsave(filename = "./PRS/img/density2PRS.png", plot = p, width = 10, height = 8)

### validate PRS using related individuals
relatedness <- read.table("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/relatedness/king/maf-prune-version1/king.kin0", header = T)[,c("ID1","ID2","Kinship","InfType")]
related1 <- merge(relatedness, PRS, by.x = "ID1", by.y = "IID")
related12 <- merge(related1, PRS, by.x = "ID2", by.y = "IID")
related12$adj_diff_PRS <- abs(related12$adj_PRS.x - related12$adj_PRS.y)
related12$raw_diff_PRS <- abs(related12$raw_PRS.x - related12$raw_PRS.y)
table(related12$InfType)
related12$Kinship <- as.numeric(related12$Kinship)
related12$InfType <- factor(related12$InfType, levels = c("FS", "PO", "2nd", "3rd", "4th", "UN"))

kin_boxplot_adj <- ggbetweenstats(data = related12, x = InfType, y = adj_diff_PRS,
                                  title = "Distribution of adjusted mtcn PRS distance across kinship type")
kin_boxplot_raw <- ggbetweenstats(data = related12, x = InfType, y = raw_diff_PRS,
                                  title = "Distribution of raw mtcn PRS distance across kinship type")
ggsave("./PRS/img/kin_boxplot_raw.png", plot = kin_boxplot_raw, width = 12, height = 6)
ggsave("./PRS/img/kin_boxplot_adj.png", plot = kin_boxplot_adj, width = 12, height = 6)

kin_scatter_adj <- ggplot(related12, aes(x=Kinship, y=adj_diff_PRS)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  labs(x="kinship coefficient", y = "adjusted mtDNA-CN PRS")
p1 <- kin_scatter_adj + stat_cor(label.x = 0.2, r.digits = 2, p.digits = 1)

kin_scatter_raw <- ggplot(related12, aes(x=Kinship, y=raw_diff_PRS)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  labs(x="kinship coefficient", y = "raw mtDNA-CN PRS")
p2 <- kin_scatter_raw + stat_cor(label.x = 0.2, r.digits = 2, p.digits = 1)

p_both <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
ggsave("./PRS/img/kin_scatter.png", plot = p_both, width = 12, height = 6)

load(file = "./mt.all.RData")
mt_PRS <- merge(mt, PRS, by.x = "ID", by.y = "IID")
#### PRS vs mitoCN estimates: whole blood samples ####
c2_min <- median(c(min(mt_PRS$mt[mt_PRS$study.cluster == "PPMI.2"]),min(mt_PRS$mt[mt_PRS$study.cluster == "PDBP.2"]), 
                   min(mt_PRS$mt[mt_PRS$study.cluster == "STEADY-PD3.2"]), min(mt_PRS$mt[mt_PRS$study.cluster == "SURE-PD.2"])))
c2_max <- median(c(max(mt_PRS$mt[mt_PRS$study.cluster == "PPMI.2"]),max(mt_PRS$mt[mt_PRS$study.cluster == "PDBP.2"]), 
                   max(mt_PRS$mt[mt_PRS$study.cluster == "STEADY-PD3.2"]), max(mt_PRS$mt[mt_PRS$study.cluster == "SURE-PD.2"])))

data <- mt_PRS[mt_PRS$cluster==2 & mt_PRS$mt > c2_min & mt_PRS$mt < c2_max,]
data$int_mt <- RankNorm(data$mt)

raw_c2 <- ggplot(data, aes(x=int_mt, y=raw_PRS) ) +
  geom_point() + 
  geom_smooth(method=lm, se=T, fullrange=F) + 
  facet_wrap(~study.cluster, scales = "free_y") + 
  labs(x="INT(mtDNA-CN)", y = "raw mtDNA-CN PRS")

adj_c2 <- ggplot(data, aes(x=int_mt, y=adj_PRS) ) +
  geom_point() + 
  geom_smooth(method=lm, se=T, fullrange=F) +
  facet_wrap(~study.cluster, scales = "free_y") + 
    labs(x="INT(mtDNA-CN)", y = "adjusted mtDNA-CN PRS")

c2_both <- ggarrange(raw_c2, adj_c2, nrow = 2, labels = c("A", "B"))

ggsave(filename = "./PRS/img/PRS_mt.c2.png", plot = c2_both, width = 12, height = 12)

sample_size <- data %>% group_by(study.cluster) %>% tally()
  
raw_mt_cor <- data %>% group_by(study.cluster) %>% summarize(cor=cor(raw_PRS, int_mt))
raw_mt_cor_p <- data %>% group_by(study.cluster) %>% summarize(cor.test(int_mt, raw_PRS)[["p.value"]])
raw_mt <- merge(raw_mt_cor, raw_mt_cor_p, by = "study.cluster")
colnames(raw_mt) <- c("study", "R_raw", "p_raw")

adj_mt_cor <- data %>% group_by(study.cluster) %>% summarize(cor=cor(adj_PRS, int_mt))
adj_mt_cor_p <- data %>% group_by(study.cluster) %>% summarize(cor.test(int_mt, adj_PRS)[["p.value"]])
adj_mt <- merge(adj_mt_cor, adj_mt_cor_p, by = "study.cluster")
colnames(adj_mt) <- c("study", "R_adj", "p_adj")

mt_cor_c2 <- cbind(merge(raw_mt, adj_mt, by = "study"), sample_size[,2])

#### PRS vs mitoCN estimates: cell line samples ####
c1_min <- median(c(min(mt_PRS$mt[mt_PRS$study.cluster == "BioFIND.1"]),min(mt_PRS$mt[mt_PRS$study.cluster == "PDBP.1"]),
                   min(mt_PRS$mt[mt_PRS$study.cluster == "PPMI.1"])))
c1_max <- median(c(max(mt_PRS$mt[mt_PRS$study.cluster == "BioFIND.1"]),max(mt_PRS$mt[mt_PRS$study.cluster == "PDBP.1"]),
                   max(mt_PRS$mt[mt_PRS$study.cluster == "PPMI.1"])))

data <- mt_PRS[mt_PRS$cluster==1 & mt_PRS$mt>c1_min & mt_PRS$mt<c1_max,]
data$int_mt <- RankNorm(data$mt)

sample_size <- data %>% group_by(study.cluster) %>% tally()

raw_mt_cor <- data %>% group_by(study.cluster) %>% summarize(cor=cor(raw_PRS, int_mt))
raw_mt_cor_p <- data %>% group_by(study.cluster) %>% summarize(cor.test(int_mt, raw_PRS)[["p.value"]])
raw_mt <- merge(raw_mt_cor, raw_mt_cor_p, by = "study.cluster")
colnames(raw_mt) <- c("study", "R_raw", "p_raw")

adj_mt_cor <- data %>% group_by(study.cluster) %>% summarize(cor=cor(adj_PRS, int_mt))
adj_mt_cor_p <- data %>% group_by(study.cluster) %>% summarize(cor.test(int_mt, adj_PRS)[["p.value"]])
adj_mt <- merge(adj_mt_cor, adj_mt_cor_p, by = "study.cluster")
colnames(adj_mt) <- c("study", "R_adj", "p_adj")

mt_cor_c1 <- cbind(merge(raw_mt, adj_mt, by = "study"), sample_size[,2])

#### PRS vs mitoCN estimates: brain samples ####
data <- mt_PRS[mt_PRS$cluster == 3 | mt_PRS$cluster == 4,]
data$int_mt <- RankNorm(data$mt)

sample_size <- nrow(data)

raw_core <- cor.test(data$raw_PRS, data$int_mt)
raw_mt <- c("LBD",raw_core$estimate, raw_core$p.value)
names(raw_mt) <- c("study", "R_raw", "p_raw")

adj_core <- cor.test(data$adj_PRS, data$int_mt)
adj_mt <- c("LBD",adj_core$estimate, adj_core$p.value)
names(adj_mt) <- c("study", "R_adj", "p_adj")

mt_cor_c3 <- c(raw_mt, adj_mt[-1], sample_size)

mt_cor <- rbind(mt_cor_c1,mt_cor_c2)
write.table(mt_cor, "./PRS/results/mtCN/PRS_mt_cor.tsv", row.names = F, quote = F, sep = "\t")

