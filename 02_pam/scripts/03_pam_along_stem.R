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
  # # add background rectangles
  # # 1 - 5 cm used for fluorescence analyses, 1 - 10 for pigments
  # annotate("rect", xmin = 1, xmax = 10, ymin = 0, ymax = 0.7, color = young, fill = young, alpha = 0) +
  # annotate("rect", xmin = 19, xmax = 71, ymin = 0, ymax = 0.7, color = old, fill = old, alpha = 0) +
  # 
  # annotate("text", x = 5.5, y = 0.735, size = 3.5, label = "Young stem", color = young) +
  # annotate("text", x = 45, y = 0.735, size = 3.5, label = "Old stem", color = old) +
  # 
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


# PSII <- ggplot() +
#   ylim (0,0.7) +
#   geom_point(data = pam_df, aes(x=Distance.from.AM, y=PSII), size=0.75, color = "gray") +
#   geom_line(data = pam_df, aes(x=Distance.from.AM, y=predict(HLM2), group=Sample), size = 0.3) + 
#   theme_light() +
#   labs(x = "Distance from apical meristem (cm)",
#        y = "phiPSII") 
# PSII
# 
# pdf("../output/boxplots/phipsii_along_stem.pdf", width=4,height=4) 
# PSII
# dev.off()


phiNPQ <- ggplot() +
  ylim (0,0.5) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0, ymax = 0.5,
           alpha = .2,fill = young) +
  annotate("rect", xmin = 5, xmax = 10, ymin = 0, ymax = 0.5,
           alpha = .075,fill = young) +
  annotate("text", x = 5.5, y = 0.08, label = "Young stem", color = young) +
  
  annotate("rect", xmin = 19, xmax = 71, ymin = 0, ymax = 0.5,
           alpha = .2,fill = old) +
  annotate("text", x = 45, y = 0.08, label = "Old stem", color = old) +
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=phiNPQ), size=1.5) +
  # geom_line(aes(y=predict(HLM3), group=Sample), size = 1, color="#ffbd00") + 
  theme_light() +
  theme(axis.title.y = element_blank())+
  labs(x = "Distance from apical meristem (cm)")
phiNPQ

pdf("../output/boxplots/npq_along_stem.pdf", width=5,height=4)  
phiNPQ
dev.off()
