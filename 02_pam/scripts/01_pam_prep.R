# http://www.sthda.com/english/wiki/r-xlsx-package-a-quick-start-guide-to-manipulate-excel-files-in-r
# https://www.oracle.com/java/technologies/javase-jdk16-downloads.html
# library(rJava)
library(xlsx)
library(tidyverse)
library(plyr)
library(rstatix)
library(gtools)
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )


# load data
samples <- read.csv("../data/pam_data_tab.csv", header = T, stringsAsFactors = FALSE) 

# remove columns of only na
samples <- samples[,colSums(is.na(samples))<nrow(samples)]

# replace "*" with "NA" # need to check if this should instead be 0?
samples[ samples == "*" ] <- NA


#### data cleaning ####
# exclude all rows that have "exclude" in the exclude column (grayed out in original excel file)
samples %>% dplyr::filter(exclude != "exclude") -> samples

# replace 'CUMO ' with 'CUMO'
samples$Accession. <- gsub('CUMO ', 'CUMO', samples$Accession.)

# replace 'CEUP ' with 'CEUP'
samples$Accession. <- gsub('CUEP ', 'CUEP', samples$Accession.)

# replace accession ss17105 17105 with ss17-105
samples$Accession. <- gsub('ss.17.105', 'ss17-105', samples$Accession.)
samples$Accession. <- gsub('ss17-105 ', 'ss17-105', samples$Accession.)
samples$Accession. <- gsub('17105', 'ss17-105', samples$Accession.)

# replace 1677 with ss16-17
samples$Accession. <- gsub('1677', 'ss16-17', samples$Accession.)
samples$Accession. <- gsub('ss.16.77', 'ss16-17', samples$Accession.)

# replace ss.13.46 with ss13-46
samples$Accession. <- gsub('ss.13.46', 'ss13-46', samples$Accession.)

# replace ss.17.15 with 17-15
samples$Accession. <- gsub('ss.17.15', '17-15', samples$Accession.)
samples$Accession. <- gsub('ss.17.15 ', '17-15', samples$Accession.)
samples$Accession. <- gsub('17-15 ', '17-15', samples$Accession.)

# replace 1346 with ss13-46
samples$Accession. <- gsub('1346', 'ss13-46', samples$Accession.)

# create a new column that simplifies the existing Tissue column, can modify if need to be more precise
samples %>%
  mutate(Tissue.edit = case_when(
    Tissue == "seedling"	~ "sdlg", # no change
    Tissue == "leaf"	~ "l", # no change
    Tissue == "young"	~ "y", # no change
    Tissue == "old"	~ "o", # no change
    Tissue == "haustorium"	~ "h", # no change
    Tissue == "flower"	~ "f", # no change
    Tissue == "seed"	~ "s", # no change
    Tissue == "fruit"	~ "fruit", # no change
    Tissue == "flower (old)"	~ "flower (old)", 
    Tissue == "seed (old)"	~ "s",
    Tissue == "Leaf"	~ "l",
    Tissue == "Corolla"	~ "Corolla", # to count as flower or omit? omitting for now
    Tissue == "seedling "	~ "sdlg")) -> samples

# check if we are really to omit corolla!
samples %>% dplyr::filter(Tissue.edit != "Corolla") -> samples

# check if we are really to omit fruit!
samples %>% dplyr::filter(Tissue.edit != "fruit") -> samples

samples %>% dplyr::filter(Tissue.edit != "flower (old)") -> samples
samples %>% dplyr::filter(Species != "C_gronovii") -> samples


# add taxonomic subgenus column, except C. purpurata which is in Grammica but kept at species level 
samples %>%
  mutate(Subgenus = case_when(
    Species == "Ipomoea_nil"	~ "Ipomoea_nil", 
    Species == "C_lupuliformis"	~ "Monogynella",
    Species == "C_monogyna"	~ "Monogynella",
    Species == "C_epithymum"	~ "Cuscuta",
    Species == "C_purpurata"	~ "C_purpurata",
    Species == "C_gracillima"	~ "Grammica",
    Species == "C_indecora"	~ "Grammica",
    Species == "C_californica"	~ "Grammica",
    Species == "C_australis"	~ "Grammica",
    Species == "C_cephalanthii"	~ "Grammica", 
    Species == "C_compacta"	~ "Grammica",
    Species == "C_denticulata"	~ "Grammica",
    Species == "C_polygonorum"	~ "Grammica", 
    Species == "C_sandwichiana"	~ "Grammica",
    Species == "C_tasmanica"	~ "Grammica",
    Species == "C_costaricensis"	~ "Grammica"), .before = Species) -> samples


# set levels to treatments
samples$Tissue.edit <- factor(samples$Tissue.edit, levels = c("l", "sdlg", "y", "o", "h", "f", "s"))

#### remove erroneous rows ####
samples <- samples[!is.na(samples$Species), ] 


#### check for outliers ####

# list of the metrics
metrics.used <- c("Fm", "Ft", "Fo", "Fo.", "Fv.Fm", "φPSII", "ΦNPQ", "Φf.D", "NPQ", "X1.qP")
# An observation with Cook’s distance larger than 4 times the mean Cook’s distance might be an outlier
# http://r-statistics.co/Outlier-Treatment-With-R.html

# set up loop
outliers <- c()
data.loop <- data.frame()
met <- NULL
Subgenus.list <- unique(samples$Subgenus)
samples$ID <- seq.int(nrow(samples)) # generate a sample ID for removing outlier

# # drop C_lupuliformis for now; not enough data to check for outliers
# Subgenus.list <- Subgenus.list[!Subgenus.list %in% "C_lupuliformis"]
# # drop C_gronovii for now; not enough data to check for outliers
# Subgenus.list <- Subgenus.list[!Subgenus.list %in% "C_gronovii"]


# loop through to find outliers per F metric per Subgenus
# first create a data frame for each Subgenus, then look for outliers in nested loop
for (sub in Subgenus.list) {
  samples %>% group_by(Subgenus) %>%
    dplyr::filter(Subgenus == sub) -> data.Subgenus
  data.Subgenus %>% group_by(Individual)

  for (met in metrics.used) {
    data.loop <- as.data.frame(cbind(ID=data.Subgenus$ID, Subgenus=data.Subgenus$Subgenus, Tissue=data.Subgenus$Tissue.edit, data.Subgenus[met]))
    data.loop <- na.omit(data.loop) # omit NA readings for data in the loop
    # proceed if >0 rows in data.loop
    if (nrow(data.loop)>1) {
    formula <- as.formula(paste(met," ~ Tissue",sep = ""))
    mod <- lm(formula, data=data.loop)
    cooksd <- cooks.distance(mod)
    # # plot
    # plot(cooksd, pch="*", cex=2, main=paste("Influential Obs by Cooks distance: ",met,sep = ""))  # plot cook's distance
    # abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
    # text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
    # n <- length(data.loop$Individual)
    influential <- as.numeric(names(cooksd)[(cooksd > (4)*mean(cooksd, na.rm=T))])  # influential row numbers
    y <- data.loop[influential,1] # identify the sample "name"ID" for the outlier reading
    outliers[[sub]][[met]] <- y
  }}}

# how many outliers for Fv.Fm, PhiPSII, Φf.D, and NPQ?
total_outliers_count <- length(outliers[["Cuscuta"]][["ΦNPQ"]]) + length(outliers[["Cuscuta"]][["Fv.Fm"]]) + length(outliers[["Cuscuta"]][["φPSII"]]) + length(outliers[["Cuscuta"]][["Φf.D"]]) +
  length(outliers[["Grammica"]][["ΦNPQ"]]) + length(outliers[["Grammica"]][["Fv.Fm"]]) + length(outliers[["Grammica"]][["φPSII"]]) + length(outliers[["Grammica"]][["Φf.D"]]) +
  length(outliers[["C_purpurata"]][["ΦNPQ"]]) + length(outliers[["C_purpurata"]][["Fv.Fm"]]) + length(outliers[["C_purpurata"]][["φPSII"]]) + length(outliers[["C_purpurata"]][["Φf.D"]]) +
  length(outliers[["Monogynella"]][["ΦNPQ"]]) + length(outliers[["Monogynella"]][["Fv.Fm"]]) + length(outliers[["Monogynella"]][["φPSII"]]) + length(outliers[["Monogynella"]][["Φf.D"]]) +
  length(outliers[["Ipomoea_nil"]][["ΦNPQ"]]) + length(outliers[["Ipomoea_nil"]][["Fv.Fm"]]) + length(outliers[["Ipomoea_nil"]][["φPSII"]]) + length(outliers[["Ipomoea_nil"]][["Φf.D"]]) 
total_outliers_count
# [1] 280 NOW 344 with Φf.D!

# will remove outliers after data is long

# convert to long format
data_long <- gather(samples, Metric, Value, Fo:X1.qP, factor_key=TRUE)

# omit NA metrics (e.g. some samples that had truly 0 signal had measurement for 1-qP)
data_long <- na.omit(data_long)
data_long$Value <- as.numeric(data_long$Value)
  
# how many total measurements? 
total_measurement_count <- length((data_long %>% dplyr::filter(Metric == "Fv.Fm"|Metric == "φPSII"| Metric == "ΦNPQ" | Metric == "Φf.D"))$Value)
total_measurement_count
# [1] 3027 NOW 4036 with Φf.D!


# loop through outliers results to remove from data_long to create data_no_outliers
data_no_outliers_list <- c()
removed.loop <- data.frame()
met <- NULL

for (sub in Subgenus.list) {
  for (met in metrics.used) {
    data_long %>% subset(Subgenus == sub & Metric == met) %>% dplyr::filter(., !(ID %in% outliers[[sub]][[met]])) -> removed.loop
    data_no_outliers_list[[sub]][[met]] <- removed.loop
  }
}


do.call(rbind, (do.call(rbind, data_no_outliers_list))) -> data_no_outliers


# some mets were causing an error because "x all 'x' values are identical"
# i identified them manually with:
data_no_outliers %>% group_by(Metric, Subgenus) %>% dplyr::summarize(sd = sd(Value)) -> test_tissues
test_tissues %>% dplyr::filter(sd == 0) -> test_tissues_0


# create table of just ipomoea lutein epoxide for correlation analysis
ipo_fvfm <- data_no_outliers %>% dplyr::filter(Species == "Ipomoea_nil") %>% dplyr::filter(Metric == "Fv.Fm") 
# save a copy of this dataset
write.csv(ipo_fvfm, file = "../output/stat_results/ipomoea_fvfm.csv", row.names = F)



#### statistical tests ####
# shapiro test for normality per Metric per subgenus
data_no_outliers %>% ungroup() %>% 
  dplyr::filter(!(Metric == "Fv.Fm" & Subgenus == "C_purpurata")) %>% 
  dplyr::filter(!(Metric == "φPSII" & Subgenus == "C_purpurata")) %>% 
  dplyr::filter(!(Metric == "NPQ" & Subgenus == "C_purpurata")) %>% 
  dplyr::filter(!(Metric == "Φf.D" & Subgenus == "C_purpurata")) %>% 
  dplyr::filter(!(Metric == "ΦNPQ" & Subgenus == "C_purpurata")) %>% 
  group_by(Metric, Subgenus) %>% 
  rstatix::shapiro_test(Value, data = .) -> shapiro_met_sub
shapiro_met_sub %>%   add_significance("p") -> shapiro_met_sub
shapiro_met_sub
# most are not normal, use kruskal-wallace

  
# kruskal wallace per Metric per subgenus
data_no_outliers %>% ungroup() %>% 
  group_by(Metric, Subgenus) %>% 
  rstatix::kruskal_test(Value ~ Tissue.edit, data = .) -> kruskal_met_sub
kruskal_met_sub

# ajust p-values
kruskal_met_sub$p.adj <- p.adjust(kruskal_met_sub$p, method = "BH", n = length(kruskal_met_sub$p))
kruskal_met_sub$p.adj.signif <- stars.pval(kruskal_met_sub$p.adj)

# save a copy of this anova
write.csv(kruskal_met_sub, file = "../output/stat_results/kruskal_met_sub.csv", row.names = F)


# filter for signficant Metrics x main effects for posthoc comparisons
kruskal_met_sub_sig <- as.data.frame(kruskal_met_sub) %>% dplyr::filter(p.adj < 0.05)

# save a copy of this anova
write.csv(kruskal_met_sub_sig, file = "../output/stat_results/kruskal_met_sub_sig.csv", row.names = F)


# perform posthoc pairwise comparisons on the sig anovas
# first create a vector on which to filter data_no_outliers by
kruskal_met_sub_sig_vector <- paste0(kruskal_met_sub_sig$Subgenus, "__", kruskal_met_sub_sig$Metric)
# create a combined column in data_no_outliers too
data_no_outliers$Subgenus.Metric <- paste0(data_no_outliers$Subgenus, "__", data_no_outliers$Metric)
# filter data long calcs by the Subgenus.Metric column, based on vector
filter(data_no_outliers, Subgenus.Metric %in% kruskal_met_sub_sig_vector) -> only_sig_met

# perform posthoc on these
# filter for one Metric x subgenus then do a pairwise wilcox on those
# after getting it working for one, set up loop to do for all combos (for those that had a positive kruskal-wallace at least)
# add these to individual plots that then get arranged with patchwork package

only_sig_met$Metric <- droplevels(only_sig_met$Metric)
subgenusMetric_list <- unique(only_sig_met$Subgenus.Metric)

wilcox_list <- list()
for(submet in subgenusMetric_list) {
  data_loop_wilcox <- only_sig_met %>% dplyr::filter(Subgenus.Metric == submet) 
  data_loop_wilcox$Tissue.edit <- droplevels(data_loop_wilcox$Tissue.edit)
  data_loop_wilcox$Metric <- droplevels(data_loop_wilcox$Metric)
  pairwise.wilcox.test(data_loop_wilcox$Value, data_loop_wilcox$Tissue.edit, p.adjust.method = "BH") -> wilcox_list[[submet]]
}

# save this list for later use 
saveRDS(wilcox_list, file="../output/stat_results/wilcox_list.RData")


#### all tissues x subgenera comparison, per metric ####
# within each metric, test for differences among genera and among tissues. achieved using the metric__subgenus
# add a column with full tissue names spelled out
data_no_outliers$Tissue.edit_names <- data_no_outliers$Tissue.edit
plyr::revalue(data_no_outliers$Tissue.edit_names, c("sdlg" = "Seedling", "l" = "Leaf", "y" = "Young", "o" = "Old", "h" = "Haustorium", "f" = "Flower", "s" = "Seed")) -> data_no_outliers$Tissue.edit_names

# kruskal wallace per metric by Subgenus.Tissue
data_no_outliers$Subgenus.Tissue <- paste0(data_no_outliers$Subgenus, "__", data_no_outliers$Tissue.edit_names)

# kruskal wallace in loop, per metric among all subgenera-tissues
kruskal_big_list <- list()
metric_list <- unique(data_no_outliers$Metric)
for(met in metric_list) {
  data_loop_kruskal <- data_no_outliers %>% dplyr::filter(Metric == met) 
  data_loop_kruskal$Tissue.edit_names <- droplevels(data_loop_kruskal$Tissue.edit_names)
  rstatix::kruskal_test(Value ~ Subgenus.Tissue, data = data_loop_kruskal) -> kruskal_big_list[[met]]
}

# all highly sig, move on to post hoc 
wilcox_big_list <- list()
metric_list <- unique(data_no_outliers$Metric)
for(met in metric_list) {
  data_loop_wilcox <- data_no_outliers %>% dplyr::filter(Metric == met) 
  data_loop_wilcox$Tissue.edit_names <- droplevels(data_loop_wilcox$Tissue.edit_names)
  pairwise.wilcox.test(data_loop_wilcox$Value, data_loop_wilcox$Subgenus.Tissue, p.adjust.method = "BH") -> wilcox_big_list[[met]]
}

# export only plotted metrics
plotted_mets <- c("Fv.Fm", "φPSII", "ΦNPQ", "Φf.D")
# format p-values, NAs, and column headers
excel_list <- list()
for(met in plotted_mets) {
  sheet_loop <- wilcox_big_list[[met]][["p.value"]] 
  rownames(sheet_loop) <- gsub(x = rownames(sheet_loop), pattern = "__", replacement = ".Tissue.")  
  colnames(sheet_loop) <- gsub(x = colnames(sheet_loop), pattern = "__", replacement = ".Tissue.")
  rownames(sheet_loop) <- gsub(x = rownames(sheet_loop), pattern = "_", replacement = ". ")  
  colnames(sheet_loop) <- gsub(x = colnames(sheet_loop), pattern = "_", replacement = ". ") 
  sheet_loop[is.na(sheet_loop)] <- "-"
  rownames(sheet_loop) <- gsub(x = rownames(sheet_loop), pattern = ".Tissue.", replacement = "__")  
  colnames(sheet_loop) <- gsub(x = colnames(sheet_loop), pattern = ".Tissue.", replacement = "__")  
  rbind(colnames(sheet_loop), sheet_loop) -> sheet_loop
  cbind(rownames(sheet_loop), sheet_loop) -> sheet_loop
  setNames(rbind(names(sheet_loop), sheet_loop), names(sheet_loop)) -> sheet_loop
  sheet_loop <- separate(as.data.frame(sheet_loop), 1, into = c("Subgenus", "Tissue"), sep = "__")
  sheet_loop <- as.data.frame(data.table::transpose(sheet_loop))
  sheet_loop <- separate(sheet_loop, 1, into = c("Subgenus", "Tissue"), sep = "__")
  sheet_loop <- data.frame(lapply(sheet_loop, function(x) {gsub("Ipomoea.", "Ipomoea", x)})) 
  sheet_loop <- data.table::transpose(sheet_loop)
  excel_list[[met]] <- sheet_loop
}

# save pairwise wilcox test with a different sheet for each metric
write.xlsx(excel_list[["Fv.Fm"]], file="../output/stat_results/pairwise_allmetsxtissues.xlsx", sheetName="Fv.Fm", row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["φPSII"]], file="../output/stat_results/pairwise_allmetsxtissues.xlsx", sheetName="φPSII", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["ΦNPQ"]], file="../output/stat_results/pairwise_allmetsxtissues.xlsx", sheetName="ΦNPQ", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Φf.D"]], file="../output/stat_results/pairwise_allmetsxtissues.xlsx", sheetName="Φf.D", append=TRUE, row.names=FALSE, col.names=FALSE)

#### summary ####


# summary stats per Individual (summarizes those that have replicates): mean, n, and std dev  
summary_individual <- data_no_outliers %>%
  group_by(Subgenus, Species, Accession., Individual, Tissue.edit, Metric)  %>%
  dplyr::summarize(mean = mean(as.numeric(Value)), n = n(), sd = sd(as.numeric(Value)))

# write to csv
write_csv(summary_individual, file = "../output/stat_results/fluorescence_individual_summary.csv")


# summary stats per species: mean, n, and std dev  
summary_species <- data_no_outliers %>%
  group_by(Subgenus, Species, Accession., Tissue.edit, Metric)  %>%
  dplyr::summarize(mean = mean(as.numeric(Value)), n = n(), sd = sd(as.numeric(Value)))

# write to csv
write_csv(summary_species, file = "../output/stat_results/fluorescence_species_summary.csv")



# summary stats per subgenus x tissue x Metric: mean, n, and std dev  
summary_subgenus <- data_no_outliers %>%
  group_by(Subgenus, Tissue.edit, Metric)  %>%
  dplyr::summarize(mean = mean(as.numeric(Value)), n = n(), sd = sd(as.numeric(Value)))

# write to csv
write_csv(summary_subgenus, file = "../output/stat_results/fluorescence_subgenus_summary.csv")




# will need to add kruskal wallace results to each box of the subgenus boxplot from kruskal_met_sub (including not significant global result)
# edit scientific notation 
kruskal_met_sub$p.adj <- kruskal_met_sub$p.adj %>% p_round(digits = 2)

# # will need to add wilcox test asterisks to each tissue in each box of the subgenus boxplot from posthoc_met_sub (including ns?)
# # edit scientific notation
# wilcox_df$p.adj <- wilcox_df$p.adj %>% p_round(digits = 2)


# #### add significance groups to wilcox results, will want to peform per metmnet x subgenus (in a loop or pipe) ####
# # for now, pick one
# wilcox_df %>% dplyr::filter(Metric == "Chl.a" & Subgenus == "Monogynella") -> wilcox_df_mono_Chl.a
# ?multcompLetters



# create a dataframe using the kruskal wallace results
# omit p of NaN
kruskal_met_sub <- kruskal_met_sub %>% dplyr::filter(p.adj != "NaN")
dat_text_plot_kruskal <- kruskal_met_sub[, c("Subgenus", "Metric", "method", "p.adj")]

# add equal sign to p.adj
dat_text_plot_kruskal$p.adj_eq <- paste0(" = ", dat_text_plot_kruskal$p.adj)
dat_text_plot_kruskal$label <- paste0('paste(italic("P"),"', dat_text_plot_kruskal$p.adj_eq, "\")") 

dat_text_plot_kruskal <- dat_text_plot_kruskal %>% select(Subgenus, Metric, label)

dat_text_plot_kruskal$Subgenus <- factor(dat_text_plot_kruskal$Subgenus, 
                                         levels = c("Ipomoea_nil",
                                                    "Monogynella",
                                                    "Cuscuta",
                                                    "C_purpurata",
                                                    "Grammica"),
                                         labels = c("Ipomoea nil",
                                                    "Monogynella",
                                                    "Cuscuta",
                                                    "C. purpurata",
                                                    "Grammica"))


# save a copy of this plot text df
write.csv(dat_text_plot_kruskal, file = "../output/stat_results/dat_text_plot_kruskal.csv", row.names = F)


# write this to csv for plotting
write_csv(data_no_outliers, file = "../output/stat_results/data_no_outliers_for_plots.csv")





#### GRAMMICA ONLY ####

#### stats to be used for Grammica only pigment plots ####

# remove species with n = 1 in Grammica 
Grammica <- c("C_australis",
              "C_polygonorum",
              "C_sandwichiana",
              "C_californica",
              "C_compacta",
              "C_cephalanthii",
              "C_denticulata",
              "C_tasmanica",
              "C_costaricensis",
              "C_gracilima",
              "C_indecora")

data_no_outliers_Grammica_plot <- dplyr::filter(data_no_outliers, Species %in% Grammica )

# identify species with n = 1
summary_species_Grammica <- dplyr::filter(summary_species, Species %in% Grammica )
average_ns <- summary_species_Grammica %>% group_by(Subgenus, Species, Accession., Tissue.edit) %>% dplyr::summarize(average_n = mean(n))

summary_accession_n1 <- dplyr::filter(average_ns, average_n == 1 )
summary_accession_n1 %>% select(Subgenus, Species, Accession., Tissue.edit) -> summary_accession_n1

# these are the specimens being removed for n = 1
Grammica_n1 <- unique(summary_accession_n1)
# add variable
Grammica_n1$drop <- paste0(Grammica_n1$Accession., Grammica_n1$Tissue.edit)

# remove species with n = 1 from plot data
data_no_outliers_Grammica_plot$drop <- paste0(data_no_outliers_Grammica_plot$Accession., data_no_outliers_Grammica_plot$Tissue.edit)

data_no_outliers_Grammica_plot %>% dplyr::filter(!drop %in% Grammica_n1$drop) -> data_no_outliers_Grammica_plot


data_no_outliers_Grammica_plot$Species <- factor(data_no_outliers_Grammica_plot$Species,
                                                 levels = c("C_australis",
                                                            "C_polygonorum",
                                                            "C_sandwichiana",
                                                            "C_californica",
                                                            "C_compacta",
                                                            "C_cephalanthii",
                                                            "C_denticulata",
                                                            "C_tasmanica",
                                                            "C_costaricensis",
                                                            "C_gracilima",
                                                            "C_indecora"))

# shapiro test for normality per fluorescence per species in Grammica only
data_no_outliers_Grammica_plot %>% ungroup() %>% 
  dplyr::filter(., Subgenus == "Grammica") %>%
  group_by(Metric, Species) %>% 
  rstatix::shapiro_test(Value, data = .) -> shapiro_met_spe_Grammica

shapiro_met_spe_Grammica %>%   add_significance("p") -> shapiro_met_spe_Grammica
shapiro_met_spe_Grammica
# many are not normal, use kruskal-wallace

# kruskal wallace per fluorescence per species in Grammica only
data_no_outliers_Grammica_plot %>% ungroup() %>%
  dplyr::filter(., Subgenus == "Grammica") %>%
  group_by(Metric, Species) %>% 
  rstatix::kruskal_test(Value ~ Tissue.edit, data = .) -> kruskal_Tissue_met_spe_Grammica
kruskal_Tissue_met_spe_Grammica

# ajust p-values
kruskal_Tissue_met_spe_Grammica$p.adj <- p.adjust(kruskal_Tissue_met_spe_Grammica$p, method = "BH", n = length(kruskal_Tissue_met_spe_Grammica$p))
kruskal_Tissue_met_spe_Grammica$p.adj.signif <- stars.pval(kruskal_Tissue_met_spe_Grammica$p.adj)

# save a copy of this anova
write.csv(kruskal_Tissue_met_spe_Grammica, file = "../output/stat_results/kruskal_Tissue_met_spe_Grammica.csv", row.names = F)


# filter for signficant fluorescence x main effects for posthoc comparisons, but round p-values first
kruskal_Tissue_met_spe_Grammica$p.adj <- kruskal_Tissue_met_spe_Grammica$p.adj %>% p_round(digits = 3)
kruskal_Tissue_met_spe_sig_Grammica <- as.data.frame(kruskal_Tissue_met_spe_Grammica) %>% dplyr::filter(p.adj <= 0.05)



# perform posthoc pairwise comparisons on the sig anovas
# first create a vector on which to filter data_no_outliers_Grammica_plot by
kruskal_Tissue_met_spe_sig_vector_Grammica <- paste0(kruskal_Tissue_met_spe_sig_Grammica$Species, "__", kruskal_Tissue_met_spe_sig_Grammica$Metric)
# create a combined column in data_no_outliers_Grammica_plot too
data_no_outliers_Grammica_plot$Species.Metric <- paste0(data_no_outliers_Grammica_plot$Species, "__", data_no_outliers_Grammica_plot$Metric)
# filter data long calcs by the Species.Metric column, based on vector
filter(data_no_outliers_Grammica_plot, Species.Metric %in% kruskal_Tissue_met_spe_sig_vector_Grammica) -> only_sig_met_Grammica

####### #######
# perform posthoc on these
# filter for one fluorescence x subgenus then do a pairwise wilcox on those
# after getting it working for one, set up loop to do for all combos (for those that had a positive kruskal-wallace at least)
# add these to individual plots that then get arranged with patchwork package

only_sig_met_Grammica$Metric <- droplevels(only_sig_met_Grammica$Metric)
specisfluorescence_list <- unique(only_sig_met_Grammica$Species.Metric)

wilcox_list_Grammica <- list()
for(spepig in specisfluorescence_list) {
  data_loop_wilcox <- only_sig_met_Grammica %>% dplyr::filter(Species.Metric == spepig) 
  data_loop_wilcox$Tissue.edit <- droplevels(data_loop_wilcox$Tissue.edit)
  data_loop_wilcox$Metric <- droplevels(data_loop_wilcox$Metric)
  pairwise.wilcox.test(data_loop_wilcox$Value, data_loop_wilcox$Tissue.edit, p.adjust.method = "BH") -> wilcox_list_Grammica[[spepig]]
}

# save this list for later use 
saveRDS(wilcox_list_Grammica, file="../output/stat_results/wilcox_list_Grammica.RData")


# will need to add kruskal wallace results to each box of the subgenus boxplot from kruskal_Tissue_met_spe (including not significant global result)
# edit scientific notation 
kruskal_Tissue_met_spe_Grammica$p.adj <- kruskal_Tissue_met_spe_Grammica$p.adj %>% p_round(digits = 3)

# # will need to add wilcox test asterisks to each tissue in each box of the subgenus boxplot from posthoc_met_spe (including ns?)


# create a dataframe using the kruskal wallace results
dat_text_plot_kruskal_Grammica <- kruskal_Tissue_met_spe_Grammica[, c("Species", "Metric", "method", "p.adj")]
# add equal sign to p.adj
dat_text_plot_kruskal_Grammica$p.adj_eq <- paste0(" = ", dat_text_plot_kruskal_Grammica$p.adj)
dat_text_plot_kruskal_Grammica$label <- paste0('paste(italic("P"),"', dat_text_plot_kruskal_Grammica$p.adj_eq, "\")") 

dat_text_plot_kruskal_Grammica <- dat_text_plot_kruskal_Grammica %>% select(Species, Metric, label)

dat_text_plot_kruskal_Grammica$Species <- factor(dat_text_plot_kruskal_Grammica$Species, 
                                                 levels = c("C_australis",
                                                            "C_polygonorum",
                                                            "C_sandwichiana",
                                                            "C_californica",
                                                            "C_compacta",
                                                            "C_cephalanthii",
                                                            "C_denticulata",
                                                            "C_tasmanica",
                                                            "C_costaricensis",
                                                            "C_gracilima",
                                                            "C_indecora"))

# save a copy of this plot text df
write.csv(dat_text_plot_kruskal_Grammica, file = "../output/stat_results/dat_text_plot_kruskal_Grammica.csv", row.names = F)

# write this to csv for plotting
write_csv(data_no_outliers_Grammica_plot, file = "../output/stat_results/data_no_outliers_Grammica_plot.csv")

