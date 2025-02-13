library(tidyverse)
library(lme4)
library(ggplot2)
library(ggeffects)
library(dplyr)
library(patchwork)
library(grid)
library(ggtext)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

young <- "orange" 
old <- "#924900" 

Fo <- "#EC5E27"
FvFm <- "#0FA379"
PhiPSII <- "#A39F1C"
PhiNPQ <- "#4B6637"

    
# load df
pam_df <- read.csv(file = "../data/Cuscuta_stem_ages_PAM.csv")


#### Hierarchical mixed linear model (random effect = sample) ####

## Fv.Fm
HLM <-lmer(Fv.Fm ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM)

## PSII
HLM2 <-lmer(PSII ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM2)

## phiNPQ
HLM3 <-lmer(phiNPQ ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM3)

## Fo
HLM4 <-lmer(Fo ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM3)

#### plot_Fv/Fm, PSII, Fo, and phiNPQ with lmer predictions ####

# select plotting data and convert to long
pam <- ggplot() +
  ylim (0,0.7) +
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fv.Fm), shape = 19, size = 2.5, stroke = 0, color = FvFm, alpha = 0.5) +
  geom_line(data = pam_df, color = FvFm, aes(x=Distance.from.AM, y=predict(HLM), group=Sample), linewidth = 0.4) + 
  
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=PSII), shape = 19, size = 2.5, stroke = 0, color = PhiPSII, alpha = 0.5) +
  geom_line(data = pam_df, color = PhiPSII, aes(x=Distance.from.AM, y=predict(HLM2), group=Sample), linewidth = 0.4) +

  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fo), shape = 19, size = 2.5, stroke = 0, color = Fo, alpha = 0.5) +
  # geom_line(data = pam_df, color = Fo, aes(x=Distance.from.AM, y=predict(HLM4), group=Sample), linewidth = 0.2) +
  
  theme_minimal() +
  theme(text = element_text(size = 14),
        strip.text.x = element_text(angle = 0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        legend.position = "none") +
  labs(x = "Distance from apical meristem (cm)", 
       y = "Fluorescence") 
pam 

pdf("../output/boxplots/pam_along_stem.pdf", width=6,height=3.5) 
pam
dev.off()


#### PAM v2 ####
pamv2 <- ggplot() +
  ylim (0,0.7) +
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fv.Fm), shape = 19, size = 2.5, stroke = 0, color = FvFm, alpha = 0.5) +
  geom_line(data = pam_df, color = FvFm, aes(x=Distance.from.AM, y=predict(HLM), group=Sample), linewidth = 0.4) + 
  
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=PSII), shape = 19, size = 2.5, stroke = 0, color = PhiPSII, alpha = 0.5) +
  geom_line(data = pam_df, color = PhiPSII, aes(x=Distance.from.AM, y=predict(HLM2), group=Sample), linewidth = 0.4) +
  
  theme_minimal() +
  theme(text = element_text(size = 14),
        strip.text.x = element_text(angle = 0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        legend.position = "none") +
  labs(x = "Distance from apical meristem (cm)", 
       y = "Fluorescence") 
pamv2 


pdf("../output/boxplots/pam_along_stem_v2.pdf", width=6,height=3.5) 
pamv2
dev.off()


#### npq ####


phiNPQ <- ggplot() +
  ylim (0,0.5) +
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=phiNPQ), shape = 19, size = 2.5, color = PhiNPQ, alpha = 0.5, stroke = 0) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        strip.text.x = element_text(angle = 0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        legend.position = "none") +
  labs(x = "Distance from apical meristem (cm)", 
       y = "PhiNPQ") 
phiNPQ

pdf("../output/boxplots/npq_along_stem.pdf", width=6,height=3.5)  
phiNPQ
dev.off()
