# library(devtools)
# devtools::install_github('cmt2/RevGadgets')
# devtools::install_github("YuLab-SMU/ggtree") 
# devtools::install_github("GuangchuangYu/treeio")
# devtools::install_github("lfabreti/convenience")
# install.packages("MCMCtrace")
library(RevGadgets)
library(coda)
library(ggplot2)
library(convenience) 
library(dplyr)
library(ggpubr)
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

lifestages <- c("seedling", "young", "old", "haustorium", "flower", "seed")

young<-"orange" 
old<-"#924900" #Brown
haustorium<- "red"
seedling<-"#009292" #PersianGreen
flower<-"#ff6db6" #HotPink
seed<- "#6db6ff" #Malibu

colors <- c(seedling, young, old, haustorium, flower, seed)

samples_combined <- list()
tissue_name <- list()
BF_correlated <- list()
BF_positive <- list()
BF_negative <- list()


for(i in 1:length(lifestages)) {
  
  filename <- paste0("../output/neo_phipsii_",lifestages[i],"/multivariate_BM.log")
  print(filename)
  
  samples <- RevGadgets::readTrace(filename)
  
  samples_combined[i] <- samples
  
  # plot densities 
  # print(RevGadgets::plotTrace(samples, vars="correlations[1]") ) 
  
  
  #### Hypothesis Tests for Correlation Parameters ####
  eta   <- 1
  c     <- 2 # number of characters
  alpha <- eta + (c - 2) / 2
  
  # choose the correlation to test
  corr  <- 1 # this is the correlation parameter (the only one in this analysis)
  
  # read the samples from the posterior
  samples <- RevGadgets::readTrace(filename)
  correlation_samples <- samples[[1]][,paste0("correlations[",corr,"]")]
  
  # fit a density to the samples
  posterior_density <- density(correlation_samples)
  
  # compute the approximate posterior probability
  # of the point hypothesis that rho_ij = 0
  post <- approxfun(posterior_density)(0)
  
  # compute the prior probability of the uncorrelated hypothesis
  # we use x = (0 + 1) / 2 = 0.5 because rho = 0 corresponds to the middle
  # of the beta distribution
  prior <- dbeta(0.5, alpha, alpha)
  
  # bayes factor for the uncorrelated hypothesis
  print("twolnBF_uncorrelated")
  twolnBF_uncorrelated <- 2*(log(post / prior))
  print(twolnBF_uncorrelated)
  
  # bayes factor for the correlated hypothesis
  print("twolnBF_correlated")
  twolnBF_correlated <- -1*twolnBF_uncorrelated
  print(twolnBF_correlated)
  
  # bayes factor for positive correlation
  # BF (positive) = (p > 0) / (1  - p > 0) (ratio of the posterior probabilities of two competing models)  all over 0.5/0.5 (the prior odds ratio: the prior probability of a particular correlation divided by the prior of the alternative, )                                                    (the prior odds ratio: the prior probability of a particular correlation divided by the prior of the alternative, he prior you are using is symmetrical, so 0.5/0.5)
  print("twolnBF_positively correlated")
  prop_pos <- sum(correlation_samples > 0)/length(correlation_samples)
  twolnBF_positive <- 2*(log((prop_pos / (1-prop_pos))/(0.5/0.5)))
  print(twolnBF_positive)
  
  # test if correlation is negative
  print("twolnBF_negatively correlated")
  prop_neg <- sum(correlation_samples < 0)/length(correlation_samples)
  twolnBF_negative <- 2*(log((prop_neg / (1-prop_neg))/(0.5/0.5)))
  print(twolnBF_negative)
  
  # save specific BF outputs to a list
  tissue_name[i] <- filename
  BF_correlated[i] <- twolnBF_correlated
  BF_positive[i] <- twolnBF_positive
  BF_negative[i] <- twolnBF_negative
}

# make BFs into data frames
file_name <- as.vector(unlist(tissue_name))
BF_correlated <- as.vector(unlist(BF_correlated))
BF_positive <- as.vector(unlist(BF_positive))
BF_negative <- as.vector(unlist(BF_negative))

stats <- data.frame(file_name, BF_correlated, BF_positive, BF_negative)

stats %>%
  dplyr::mutate("Life stage" = case_when(
    file_name == "../output/neo_phipsii_seedling/multivariate_BM.log" ~ "1", 
    file_name == "../output/neo_phipsii_young/multivariate_BM.log" ~ "2", 
    file_name == "../output/neo_phipsii_old/multivariate_BM.log" ~ "3", 
    file_name == "../output/neo_phipsii_haustorium/multivariate_BM.log" ~ "4", 
    file_name == "../output/neo_phipsii_flower/multivariate_BM.log" ~ "5", 
    file_name == "../output/neo_phipsii_seed/multivariate_BM.log" ~ "6"), .before = file_name) %>% dplyr::select("Life stage", BF_correlated, BF_positive, BF_negative) -> stats

# add significance asterisks to stats table, a la  kass and raftery 1995
# 0 to 2 ns
# 2 to 6 *
# 6 to 10 **
# > 10 ***

stats$sig_pos[stats$BF_positive >0 & stats$BF_positive <2] <- "n.s."
stats$sig_pos[stats$BF_positive >2 & stats$BF_positive <6] <- "*"
stats$sig_pos[stats$BF_positive >6 & stats$BF_positive <10] <- "**"
stats$sig_pos[stats$BF_positive >10] <- "***"

stats$sig_neg[stats$BF_negative >0 & stats$BF_negative <2] <- "n.s"
stats$sig_neg[stats$BF_negative >2 & stats$BF_negative <6] <- "*"
stats$sig_neg[stats$BF_negative >6 & stats$BF_negative <10] <- "**"
stats$sig_neg[stats$BF_negative >10] <- "***"


# combine data into one df for a density plot of all life stages 
seedling_corr <- samples_combined[[1]][["correlations[1]"]]
young_corr <- samples_combined[[2]][["correlations[1]"]]
old_corr <- samples_combined[[3]][["correlations[1]"]]
haustorium_corr <- samples_combined[[4]][["correlations[1]"]]
flower_corr <- samples_combined[[5]][["correlations[1]"]]
seed_corr <- samples_combined[[6]][["correlations[1]"]]

corrs_all <- data.frame(seedling_corr,young_corr,old_corr,haustorium_corr,flower_corr,seed_corr)
colnames(corrs_all) <- c("Seedling", "Young stem", "Old stem", "Haustorium", "Flower", "Seed")

# density_hist <- RevGadgets::plotTrace(list(corrs_all), vars=c(colnames(corrs_all)), color = colors)[[1]]
# 
# density_hist +
#   xlab("Correlation") +
#   scale_fill_manual(values=alpha(colors, .5)) +
#   scale_color_manual(values=colors) +
#   geom_vline(xintercept=0) +
#   guides(color=guide_legend(title="Life stage")) +
#   ggtitle(NULL) +
#   theme(legend.position = c(0.9,0.75)) -> density_hist
# 
# pdf("correlations_density_hist.pdf", width=8,height=5)
# density_hist
# dev.off()


# plot instead as violin plot using the same color palette RevGadgets::colFun() 
corrs_all_long <- tidyr::gather(corrs_all, Life.stage, Correlation, Seedling:Seed, factor_key=TRUE)
# corrs_all_long$Life.stage <- factor(corrs_all_long$Life.stage, levels = c("Seed", "Flower", "Haustorium", "Old stem", "Young stem", "Seedling"))

corrs_all_long %>% group_by(Life.stage) %>% ggplot(., aes(x=Correlation, y = Life.stage, fill = Life.stage, color = Life.stage)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_violin(draw_quantiles = c(0.025, 0.975)) + 
  ggthemes::theme_few() +
  theme(legend.position = "none") +
  scale_fill_manual(values=alpha(colors, .5)) + 
  scale_color_manual(values=colors) + 
  ylab("Life stage") -> density_viol


# prepare data for plotting violins

df <- corrs_all_long
df %>%
  dplyr::mutate("grp" = case_when(
    Life.stage == "Seedling" ~ "1", 
    Life.stage == "Young stem" ~ "2", 
    Life.stage == "Old stem" ~ "3", 
    Life.stage == "Haustorium" ~ "4", 
    Life.stage == "Flower" ~ "5", 
    Life.stage == "Seed" ~"6")) %>% dplyr::select(grp, Correlation) -> df
colnames(df) <- c("grp", "val")

# write.csv(df, file = "violins_data.csv", row.names = F) # sample data set for violins.R


# now that we have df, add y.position to stats table for asterisks
stats$y.position_pos <- (df %>% group_by(grp) %>% summarize(max = max(val)))$max + 0.055
stats$y.position_neg <- (df %>% group_by(grp) %>% summarize(min = min(val)))$min - 0.04

# prepare for stat_pvalue_manual
names(stats)[names(stats) == 'Life stage'] <- 'group1'
stats$group2 <- rep("null model")





#### PLOT ####
g <- ggplot(df) +
  geom_violin(aes(x = as.factor(grp), y = val, group = grp, color = grp), fill = NA) +
  geom_hline(yintercept = 0.0, linetype = "dashed")  # build the base violins, make them invisble fill!

coords <- ggplot_build(g)$data        # use ggbuild to get the outline co-ords
d <- coords[[1]]                      # this gets the df in a usable form
groups <- unique(d$group)             # get the unique "violin" ids

# function to create geom_ploygon calls
fill_viol<-function(v,gr){
  quants<-mutate(v, x.l = x - violinwidth*0.45, x.r = x + violinwidth*0.45, 
                 cuts = cut(y, quantile(df[df$grp == gr,"val"],
                                        probs = c(0, 0.025, 0.975, 1)))) # add 1/2 (changed from 0.5 to 0.45) width each way to each x value, change to only 95% CI
  plotquants<-data.frame(x=c(quants$x.l,rev(quants$x.r)),   # left x bottom to top, then right x top to bottom
                         y=c(quants$y,rev(quants$y)),       # double up the y values to match
                         id=c(quants$cuts,rev(quants$cuts)),# cut by quantile to create polygon id
                         grp = as.factor(rep(gr)))
  plotquants$interaction <- interaction(as.numeric(plotquants$id), plotquants$grp)
  plotquants <- dplyr::filter(plotquants, !grepl("1.",interaction)) # drop bottom 2.5% 
  plotquants <- dplyr::filter(plotquants, !grepl("3.",interaction)) # drop top 2.5% 
  geom_polygon(aes(x,y,fill=interaction(as.numeric(plotquants$id), plotquants$grp)),data=plotquants) # return the geom_ploygon object
}


p <- g +                                                      # plot g (empty violin plots)
  lapply(groups,function(x)fill_viol(d[d$group==x,],x)) +     # plus 95% polygon object for each violin
  scale_fill_manual(values=alpha(colors, .5)) # plus fill


p +                                 # add style specifications to p
  coord_flip() +
  ggthemes::theme_few() +
  theme(legend.position = "none") +
  scale_color_manual(values=colors) + 
  ylab("Correlation") +
  xlab("Life stage") + 
  scale_x_discrete(labels=c("1" = "Seedling", 
                            "2" = "Young stem",
                            "3" = "Old stem",
                            "4" = "Haustorium",
                            "5" = "Flower",
                            "6" = "Seed")) -> p

p + 
  ggpubr::stat_pvalue_manual(data = stats, label = "sig_pos", x = "group1", y.position = "y.position_pos") +
  ggpubr::stat_pvalue_manual(data = stats, label = "sig_neg", x = "group1", y.position = "y.position_neg") -> density_viol
density_viol





pdf("../output/plots/correlations_neo_phipsii_density_violin.pdf", width=6,height=5) 
density_viol
dev.off()

