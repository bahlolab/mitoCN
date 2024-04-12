rm(list = ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/mtSwirl/mtcn/")
library(ggplot2)
library(plyr)
library(ggpubr) # for multiple graphs

BioFIND <- read.table(paste0("BioFIND","/mtcn.txt"), header = T, sep = "\t")
Sure <- read.table(paste0("Sure","/mtcn.txt"), header = T, sep = "\t")
Steady <- read.table(paste0("Steady","/mtcn.txt"), header = T, sep = "\t")
PPMI <- read.table(paste0("PPMI","/mtcn.txt"), header = T, sep = "\t")
PDBP <- read.table(paste0("PDBP","/mtcn.txt"), header = T, sep = "\t")
HBS <- read.table(paste0("HBS","/mtcn.txt"), header = T, sep = "\t")
mtSwirl <- rbind(BioFIND,Sure,Steady,PPMI,PDBP,HBS)
write.table(mtSwirl, "mtSwirl_mtcn.txt", row.names = F, quote = F, sep = "\t")
load("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/meta_analysis/mt.all.RData")
mitoCN <- mt[,c("Np","mt","ID","cluster")]
data <- merge(mitoCN, mtSwirl, by.x = "ID", by.y = "ID")

### concordance
cor_test <- cor.test(data$mt, data$mtcn)
p <- cor_test$p.value
R2 <- round(cor_test$estimate,3)

### scatter plot
cor_all <- ggplot(data, aes(x=mt, y=mtcn)) + 
  geom_point() + 
  geom_smooth(method=lm, fullrange=TRUE) +
  xlab("mitoCN") + ylab("mtSwirl") + 
  ggtitle(paste0("AMP-PD (R2 = ",R2,")"))

### percentage mtcn change with mtSwirl
data$diff <- data$mtcn - data$mt
data$pct_diff <- data$diff/data$mt

## histogram by cluster
mu_cluster <- ddply(data, "cluster", summarise, grp.mean=mean(pct_diff))
diff_cluster <- ggplot(data, aes(x=pct_diff, fill=cluster, color=cluster)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(data=mu_cluster, aes(xintercept=grp.mean, color=cluster), linetype="dashed") + 
  xlab("percentage mtcn change with mtSwirl") + 
  theme(legend.position="top")

t.test(data$pct_diff[data$cluster==1])
t.test(data$pct_diff[data$cluster==2])

### compare by ancestry
ancestry <- read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/ancestry/version2/PC_table_with_estimated_superpop_PC123.txt")[,c("sample","ancestry")]
cmp_ancestry <- merge(data, ancestry, by.x = "ID", by.y = "sample")
cmp_ancestry <- cmp_ancestry[cmp_ancestry$ancestry != "ND",]
table(cmp_ancestry$ancestry)
mu_ancestry <- ddply(cmp_ancestry, c("ancestry","cluster"), summarise, grp.mean=mean(pct_diff))

diff_ancestry <-ggplot(cmp_ancestry, aes(x=pct_diff, fill=cluster, color=cluster) ) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(data=mu_ancestry, aes(xintercept=grp.mean, color=cluster), linetype="dashed") + 
  theme_bw() + xlab("percentage mtcn change with mtSwirl") + 
  facet_wrap(~ancestry, scales = "free_y")

mu_ancestry[mu_ancestry$cluster == 2 & mu_ancestry$grp.mean > -0.004,]

p1 <- ggarrange(cor_all, diff_cluster, nrow = 2, labels = c("A", "B"))
fig <- ggarrange(p1, diff_ancestry, ncol = 2, labels = c("","C"), widths = c(0.9,2))
ggsave(filename = "fig.mtswirl.png", plot = fig, width = 13, height = 8)
