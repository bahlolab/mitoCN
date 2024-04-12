rm(list=ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")
library(tidyverse)
library(meta)
library(metaforest)
library(dplyr)
library(ggplot2)

mlabfun <- function(text, x) { list(bquote(paste(.(text),
            " (Q = ", .(fmtx(x$QE, digits=2)), ", ",
            .(fmtp(x$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), ")")))
  }
#################################################################################
#### PD diagnosis
PD.stat <- read.table("./meta_analysis/res/PD.rlm.tsv", header = T, sep = "\t")
PD.re <- rma(yi = Estimate, sei= Std..Error, data = PD.stat)
PD.re

png(filename = "./meta_analysis/img/meta_DX.png", width = 650, height = 400)
forest(PD.re, xlim=c(-1, 0.4), slab = PD.stat$study,  digits = 2,
       cex=1, ylim=c(-1, 17), order=PD.stat$cluster, rows=c(11:13,4:6),
       psize=1, header="Study", mlab="", xlab = "", 
       ilab = cbind(PD.stat$Case, PD.stat$Control), ilab.xpos = c(-0.7,-0.4))

par(cex=1, font=2)
text(c(-0.7,-0.4), 16, c("Case","Control"))
text(-0.55, 17, "mtDNA-CN: mean(SD)")

par(font=4)
text(-1, c(14,7), pos=4, c("Cluster 1", "Cluster 2"))

### fit random-effects model in the two subgroups
PD1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=PD.stat)
PD2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=PD.stat)

par(font=1)
text(-1, pos = 4, c(9.5,2.5), 
     c(paste0("RE Model: p-value = ",round(PD1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(PD2.re$pval,3)))
)

### add summary polygons for the two subgroups
addpoly(PD1.re, row=8.5, mlab=mlabfun("Heterogeneity ", PD1.re))
addpoly(PD2.re, row=1.5, mlab=mlabfun("Heterogeneity ", PD2.re))

dev.off()

rma(yi = Estimate, sei= Std..Error, data=PD.stat[PD.stat$cluster == 2 & PD.stat$study != "PDBP", ])

#################################################################################
#### MDS_UPDRS_III
MDS_UPDRS.stat <- read.table("./meta_analysis/res/MDS_UPDRS_III.rlm.tsv", header = T, sep = "\t")
MDS_UPDRS.re <- rma(yi = Estimate, sei= Std..Error, data = MDS_UPDRS.stat)
MDS_UPDRS.re

png(filename = "./meta_analysis/img/meta_MDS_UPDRS_III.png", width = 650, height = 400)
forest(MDS_UPDRS.re, xlim=c(-0.02, 0.012), slab = MDS_UPDRS.stat$study,  digits = 3,
       cex=1, ylim=c(-1, 18), order=MDS_UPDRS.stat$cluster, rows=c(12:14,4:7),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(MDS_UPDRS.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.011, 17, c("P-value"))
par(font=4)
text(-0.02, c(15,8), pos=4, c("Cluster 1", "Cluster 2"))

### fit random-effects model in the two subgroups
MDS_UPDRS1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=MDS_UPDRS.stat)
MDS_UPDRS2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=MDS_UPDRS.stat)

par(font=1)
text(-0.02, pos = 4, c(10.5,2.5), 
     c(paste0("RE Model: p-value = ",round(MDS_UPDRS1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(MDS_UPDRS2.re$pval,3)))
)

### add summary polygons for the two subgroups
addpoly(MDS_UPDRS1.re, row=9.5, mlab=mlabfun("Heterogeneity ", MDS_UPDRS1.re))
addpoly(MDS_UPDRS2.re, row=1.5, mlab=mlabfun("Heterogeneity ", MDS_UPDRS2.re))

dev.off()

rma(yi = Estimate, sei= Std..Error, data=MDS_UPDRS.stat[MDS_UPDRS.stat$cluster == 2 & MDS_UPDRS.stat$study != "BioFIND", ])

#################################################################################
#### ADL
ADL.stat <- read.table("./meta_analysis/res/ADL.rlm.tsv", header = T, sep = "\t")
ADL.re <- rma(yi = Estimate, sei= Std..Error, data = ADL.stat)
ADL.re

png(filename = "./meta_analysis/img/meta_ADL.png", width = 650, height = 400)
forest(ADL.re, xlim=c(-0.02, 0.02), slab = ADL.stat$study,  digits = 3,
       cex=1, ylim=c(-1, 18), order=ADL.stat$cluster, rows=c(12:14,4:7),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(ADL.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.009, 17, c("P-value"))
par(font=4)
text(-0.02, c(15,8), pos=4, c("Cluster 1", "Cluster 2"))

### fit random-effects model in the two subgroups
ADL1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=ADL.stat)
ADL2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=ADL.stat)

par(font=1)
text(-0.02, pos = 4, c(10.5,2.5), 
     c(paste0("RE Model: p-value = ",round(ADL1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(ADL2.re$pval,3)))
)

### add summary polygons for the two subgroups
addpoly(ADL1.re, row=9.5, mlab=mlabfun("Heterogeneity ", ADL1.re))
addpoly(ADL2.re, row=1.5, mlab=mlabfun("Heterogeneity ", ADL2.re))

dev.off()

#################################################################################
#### MOCA
MOCA.stat <- read.table("./meta_analysis/res/MOCA.rlm.tsv", header = T, sep = "\t")
MOCA.re <- rma(yi = Estimate, sei= Std..Error, data = MOCA.stat)
MOCA.re

png(filename = "./meta_analysis/img/meta_MOCA.png", width = 650, height = 400)
forest(MOCA.re, xlim=c(-0.05, 0.05), slab = MOCA.stat$study,  digits = 3,
       cex=1, ylim=c(-1, 18), order=MOCA.stat$cluster, rows=c(12:14,4:7),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(MOCA.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.026, 17, c("P-value"))
par(font=4)
text(-0.05, c(15,8), pos=4, c("Cluster 1", "Cluster 2"))

### fit random-effects model in the two subgroups
MOCA1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=MOCA.stat)
MOCA2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=MOCA.stat)

par(font=1)
text(-0.05, pos = 4, c(10.5,2.5), 
     c(paste0("RE Model: p-value = ",round(MOCA1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(MOCA2.re$pval,3)))
)

### add summary polygons for the two subgroups
addpoly(MOCA1.re, row=9.5, mlab=mlabfun("Heterogeneity ", MOCA1.re))
addpoly(MOCA2.re, row=1.5, mlab=mlabfun("Heterogeneity ", MOCA2.re))

dev.off()

#################################################################################
#### UPSIT
UPSIT.stat <- read.table("./meta_analysis/res/UPSIT.rlm.tsv", header = T, sep = "\t")
UPSIT.re <- rma(yi = Estimate, sei= Std..Error, data = UPSIT.stat)
UPSIT.re

png(filename = "./meta_analysis/img/meta_UPSIT.png", width = 650, height = 400)
forest(UPSIT.re, xlim=c(-0.02, 0.02), slab = UPSIT.stat$study,  digits = 3,
       cex=1, ylim=c(-1, 15), order=UPSIT.stat$cluster, rows=c(10:11,4:5),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(UPSIT.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.009, 14, c("P-value"))
par(font=4)
text(-0.02, c(12,6), pos=4, c("Cluster 1", "Cluster 2"))

### fit random-effects model in the two subgroups
UPSIT1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=UPSIT.stat)
UPSIT2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=UPSIT.stat)

par(font=1)
text(-0.02, pos = 4, c(8.5,2.5), 
     c(paste0("RE Model: p-value = ",round(UPSIT1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(UPSIT2.re$pval,3)))
)

### add summary polygons for the two subgroups
addpoly(UPSIT1.re, row=7.5, mlab=mlabfun("Heterogeneity ", UPSIT1.re))
addpoly(UPSIT2.re, row=1.5, mlab=mlabfun("Heterogeneity ", UPSIT2.re))

dev.off()

