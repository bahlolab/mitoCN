rm(list = ls())
##### generate main Figure #####
library(forestplot)
library(dplyr)
library(readxl)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(ggpubr) #for multiple plots
library(ggplotify)
library(patchwork)

### a. forest plot
MR <- read_excel("/Users/wang.lo/Documents/MR.xlsx", sheet = "Causal effect estimation", range = "E3:K10", col_names = F)
IVW <- MR[,1:4]; colnames(IVW) <- c("SNP", "beta", "se", "p")
Egger <- MR[,c(1,5:7)]; colnames(Egger) <- c("SNP", "beta", "se", "p")
est <- rbind(IVW, Egger)
est$Method <- rep(c("IVW","MR Egger"), each = 8)
mtcn_names <- c("mtDNA-CN (WGS,raw)", "mtDNA-CN (WGS,adjusted)", "mtDNA-CN (genotyping)", "mtDNA-CN (WES+genotyping)")
est$exposure <- rep(c(mtcn_names, rep("PD", 4)), 2)
est$outcome <- rep(c(rep("PD", 4), mtcn_names), 2)
est$mean <- est$beta
est$lower <- est$beta - 1.96*est$se
est$upper <- est$beta + 1.96*est$se
est$p <- round(est$p, 2)
est$beta <- round(est$beta, 2)
data <- est[est$Method == "IVW",]

fig.a <- data |> 
  filter(outcome == "PD") |> 
  forestplot(labeltext = c(exposure, outcome, SNP, beta, p)) |>
  fp_add_lines() |> 
  fp_set_style(box = "royalblue", line = "darkblue") |> 
  fp_add_header(exposure = "exposure",
                outcome = "outcome",
                SNP = "#SNPs",
                beta = "beta",
                p = "p value") |>
  fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                    graph.pos = 5) |> 
  fp_set_zebra_style("#EFEFEF")

fig.b <- data |> 
  filter(exposure == "PD") |> 
  forestplot(labeltext = c(exposure, outcome, SNP, beta, p)) |>
  fp_add_lines() |> 
  fp_set_style(box = "royalblue", line = "darkblue") |> 
  fp_add_header(exposure = "exposure",
                outcome = "outcome",
                SNP = "#SNPs",
                beta = "beta",
                p = "p value") |>
  fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                    graph.pos = 5) |> 
  fp_set_zebra_style("#EFEFEF")

p1 <- grid2grob(print(fig.a))
p2 <- grid2grob(print(fig.b))
p_both <- ggarrange(p1, p2, nrow = 2, labels = c("A", "B"))
ggsave(filename = "/Users/wang.lo/Documents/MR.main.fig.png", plot = p_both, dpi = 500)
