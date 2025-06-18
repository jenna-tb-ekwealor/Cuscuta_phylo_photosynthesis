library(tidyverse)
library(lme4)
library(ggeffects)
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
# Linear mixed model fit by REML ['lmerMod']
# Formula: Fv.Fm ~ log(Distance.from.AM) + (1 | Biorep/Sample)
# Data: pam_df
# 
# REML criterion at convergence: -138
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.98695 -0.46294 -0.05664  0.47688  2.07186 
# 
# Random effects:
#   Groups        Name        Variance  Std.Dev.
# Sample:Biorep (Intercept) 0.0004156 0.02039 
# Biorep        (Intercept) 0.0002910 0.01706 
# Residual                  0.0010113 0.03180 
# Number of obs: 40, groups:  Sample:Biorep, 7; Biorep, 2
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            0.630617   0.017546   35.94
# log(Distance.from.AM) -0.032074   0.003616   -8.87
# 
# Correlation of Fixed Effects:
#   (Intr)
# lg(Dst..AM) -0.496


## PSII
HLM2 <-lmer(PSII ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM2)
# Linear mixed model fit by REML ['lmerMod']
# Formula: PSII ~ log(Distance.from.AM) + (1 | Biorep/Sample)
# Data: pam_df
# 
# REML criterion at convergence: -123.8
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -3.16000 -0.71542  0.09017  0.70156  2.07891 
# 
# Random effects:
#   Groups        Name        Variance  Std.Dev.
# Sample:Biorep (Intercept) 9.189e-06 0.003031
# Biorep        (Intercept) 1.237e-03 0.035173
# Residual                  1.692e-03 0.041134
# Number of obs: 40, groups:  Sample:Biorep, 7; Biorep, 2
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            0.318426   0.028098  11.333
# log(Distance.from.AM) -0.022608   0.004666  -4.845
# 
# Correlation of Fixed Effects:
#   (Intr)
# lg(Dst..AM) -0.399


## phiNPQ
HLM3 <-lmer(phiNPQ ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM3)
# # poor model fit
# Linear mixed model fit by REML ['lmerMod']
# Formula: phiNPQ ~ log(Distance.from.AM) + (1 | Biorep/Sample)
# Data: pam_df
# 
# REML criterion at convergence: -92.3
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.9584 -0.5530 -0.1368  0.8331  2.0378 
# 
# Random effects:
#   Groups        Name        Variance Std.Dev.
# Sample:Biorep (Intercept) 0.000000 0.00000 
# Biorep        (Intercept) 0.001223 0.03497 
# Residual                  0.003968 0.06300 
# Number of obs: 40, groups:  Sample:Biorep, 7; Biorep, 2
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            0.311721   0.031764   9.814
# log(Distance.from.AM) -0.012427   0.007138  -1.741
# 
# Correlation of Fixed Effects:
#   (Intr)
# lg(Dst..AM) -0.541
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')

## Fo
HLM4 <-lmer(Fo ~ log(Distance.from.AM) + (1|Biorep/Sample), data=pam_df)
summary(HLM4)
# # poor model fit
# Linear mixed model fit by REML ['lmerMod']
# Formula: Fo ~ log(Distance.from.AM) + (1 | Biorep/Sample)
# Data: pam_df
# 
# REML criterion at convergence: -253.6
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.5540 -0.7208 -0.1064  0.7208  1.8562 
# 
# Random effects:
#   Groups        Name        Variance  Std.Dev. 
# Sample:Biorep (Intercept) 6.142e-05 7.837e-03
# Biorep        (Intercept) 2.233e-13 4.725e-07
# Residual                  4.221e-05 6.497e-03
# Number of obs: 40, groups:  Sample:Biorep, 7; Biorep, 2
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            0.0587587  0.0036158   16.25
# log(Distance.from.AM) -0.0088097  0.0007394  -11.91
# 
# Correlation of Fixed Effects:
#   (Intr)
# lg(Dst..AM) -0.496
# optimizer (nloptwrap) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')



#### plot_Fv/Fm, PSII, Fo, and phiNPQ with appropriate lmer predictions ####
# select plotting data and convert to long
pam <- ggplot() +
  ylim (0,0.7) +
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fv.Fm), shape = 19, size = 2.5, stroke = 0, color = FvFm, alpha = 0.5) +
  geom_line(data = pam_df, color = FvFm, aes(x=Distance.from.AM, y=predict(HLM), group=Sample), linewidth = 0.4) + 
  
  geom_point(data = pam_df, aes(x=Distance.from.AM, y=PSII), shape = 19, size = 2.5, stroke = 0, color = PhiPSII, alpha = 0.5) +
  geom_line(data = pam_df, color = PhiPSII, aes(x=Distance.from.AM, y=predict(HLM2), group=Sample), linewidth = 0.4) +

  geom_point(data = pam_df, aes(x=Distance.from.AM, y=Fo), shape = 19, size = 2.5, stroke = 0, color = Fo, alpha = 0.5) +
  # geom_line(data = pam_df, color = Fo, aes(x=Distance.from.AM, y=predict(HLM4), group=Sample), linewidth = 0.2) +
  # Not plotting HLM4 predictions due to poor model fit
  
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
