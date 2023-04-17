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

young<-"orange" 
old<-"#924900" #Brown
    
# load df
pam_df <- read.csv(file = "../data/Cuscuta_stem_ages_PAM.csv")


#### Hierarchical mixed linear model (random effect = sample) ####
## Fv.Fm
HLM <-lmer(Fv.Fm ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM)

## now repeat for PSII
HLM2 <-lmer(PSII ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM2)

## now repeat for phiNPQ
HLM3 <-lmer(phiNPQ ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM3)

#### plot_Fv/Fm, PSII, and phiNPQ with lmer predictions ####
# select plotting data and convert to long


pam <- ggplot() +
  ylim (0,0.7) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0, ymax = 0.7,
           alpha = .2,fill = young) +
  annotate("rect", xmin = 5, xmax = 10, ymin = 0, ymax = 0.7,
           alpha = .075,fill = young) +
  annotate("text", x = 5.5, y = 0.08, label = "Young stem", color = young) +
  
  annotate("rect", xmin = 19, xmax = 71, ymin = 0, ymax = 0.7,
           alpha = .2,fill = old) +
  annotate("text", x = 45, y = 0.08, label = "Old stem", color = old) +
  
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fv.Fm), size=1.5, color = "darkgray") +
  geom_line(data = pam_df, aes(x=Distance.from.AM, y=predict(HLM), group=Sample), linewidth = 0.5) + 
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=PSII), shape = 17, size=1.6, color = "darkgray") +
  geom_line(data = pam_df, aes(x=Distance.from.AM, y=predict(HLM2), group=Sample), linetype = 'dashed', linewidth = 0.5) + 
  theme_light() +
  theme(axis.title.y = element_blank()) +
  labs(x = "Distance from apical meristem (cm)") 
pam 

pdf("../output/boxplots/pam_along_stem.pdf", width=6,height=4) 
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
