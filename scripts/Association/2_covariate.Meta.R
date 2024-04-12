rm(list=ls())
setwd("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/")
library(tidyverse)
library(meta)
library(metaforest)
library(dplyr)
library(ggplot2)

mlabfun <- function(text, x) {
  list(bquote(paste(.(text),
                    " (Q = ", .(fmtx(x$QE, digits=2)), ", ",
                    .(fmtp(x$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), ")")))}


#### age
age.stat <- read.table("./meta_analysis/res/age.rlm.tsv", header = T, sep = "\t")
age.re <- rma(yi = Estimate, sei= Std..Error, data = age.stat)
age.re

png(filename = "./meta_analysis/img/meta_age.png", width = 800, height = 500)
forest(age.re, xlim=c(-0.02, 0.015), slab = age.stat$study, digits = 3,
       cex=1, ylim=c(-1, 19), order=age.stat$cluster, rows=c(13:15,4:8),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(age.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.01, 18, c("P-value"))

par(font=4)
text(-0.02, c(16,9), pos=4, c("Cluster 1", "Cluster 2"))

age1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=age.stat)
age2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=age.stat)

par(font=1)
text(-0.02, pos = 4, c(11.5,2.5), 
     c(paste0("RE Model: p-value = ",round(age1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(age2.re$pval,3)))
)

addpoly(age1.re, row=10.5, mlab=mlabfun("Heterogeneity ", age1.re))
addpoly(age2.re, row=1.5, mlab=mlabfun("Heterogeneity ", age2.re))

dev.off()

#### sex
sex.stat <- read.table("./meta_analysis/res/sex.rlm.tsv", header = T, sep = "\t")
sex.re <- rma(yi = Estimate, sei= Std..Error, data = sex.stat)
sex.re

png(filename = "./meta_analysis/img/meta_sex.png", width = 800, height = 500)
forest(sex.re, xlim=c(-0.3, 0.15), slab = sex.stat$study,  digits = 2,
       cex=1, ylim=c(-1, 19), order=sex.stat$cluster, rows=c(13:15,4:8),
       psize=1, header="Study", mlab="", xlab = "", ilab = round(sex.stat$Pr...t.., 3))

par(cex=1, font=2)
text(-0.18, 18, c("P-value"))
par(font=4)
text(-0.3, c(16,9), pos=4, c("Cluster 1", "Cluster 2"))

sex1.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==1), data=sex.stat)
sex2.re <- rma(yi = Estimate, sei= Std..Error, subset=(cluster==2), data=sex.stat)

par(font=1)
text(-0.3, pos = 4, c(11.5,2.5), 
     c(paste0("RE Model: p-value = ",round(sex1.re$pval,3)), 
       paste0("RE Model: p-value = ",round(sex2.re$pval,3)))
)

addpoly(sex1.re, row=10.5, mlab=mlabfun("Heterogeneity ", sex1.re))
addpoly(sex2.re, row=1.5, mlab=mlabfun("Heterogeneity ", sex2.re))

dev.off()
