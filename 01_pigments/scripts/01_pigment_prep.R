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
# skip the first row which has units info (check raw csv for units)
samples <- read.csv("../data/pigment_organized_data_tab_AS_edit.csv", skip = 1, header = T, stringsAsFactors = FALSE, na.strings = c("NA","NaN","","#DIV/0!")) 


#### data cleaning ####
# replace 'seedling' with 'sdlg'
samples$Tissue.code <- gsub('seedling', 'sdlg', samples$Tissue.code)

# replace 'CUMO ' with 'CUMO'
samples$Accession.No <- gsub('CUMO ', 'CUMO', samples$Accession.No)

# replace 'CEUP ' with 'CEUP'
samples$Accession.No <- gsub('CUEP ', 'CUEP', samples$Accession.No)

# replace accession ss17105 17105 with ss17-105
samples$Accession.No <- gsub('ss17105', 'ss17-105', samples$Accession.No)
samples$Accession.No <- gsub('17105', 'ss17-105', samples$Accession.No)

# replace 1677 with ss16-17
samples$Accession.No <- gsub('1677', 'ss16-17', samples$Accession.No)

# replace 1715 with 17-15
samples$Accession.No <- gsub('1715', '17-15', samples$Accession.No)

# replace 1346 with ss13-46
samples$Accession.No <- gsub('1346', 'ss13-46', samples$Accession.No)



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
samples$Tissue.code <- factor(samples$Tissue.code, levels = c("l", "sdlg", "y", "o", "h", "f", "s"))


#### select ng pigments, remove peak.11.9 column, omit d-Tocopherol because no standard and it's in area not ng ####
samples %>% select( -peak.11.9) %>% select( -d.Tocopherol) %>% select(Date.Extracted:Chl.a, Car.Chl31) -> samples

#### remove erroneous rows ####
samples_naomit <- samples[!is.na(samples$Species), ]


#### normalize to fresh weight ####
# convert to long format
data_long <- gather(samples_naomit, Pigment, ng, a.Tocopherol:Car.Chl31, factor_key=TRUE)

# normalize to fresh weight, ng/mg FW
data_long %>% mutate(., FW.norm = ng*Conversion) -> data_long


#### create relevant ratios and pigment sums ####
# drop ng column
data_long %>% select(-ng) -> data_long

# convert to wide again for easier transformations of columns 
data_wide <- spread(data_long, Pigment, FW.norm)

# add total chlorophyll column
data_wide %>% mutate(Tot.Chl = Chl.a + Chl.b) -> data_wide

# add chlorophyll ratio
data_wide %>% mutate(Chl.a.b = Chl.a/Chl.b) -> data_wide

# add total carotenoid
data_wide %>% mutate(Tot.Car = Violaxanthin + Neoxanthin + Antheraxanthin + Lutein + Zeaxanthin + a.Carotene + b.Carotene) -> data_wide

# add VAZ
data_wide %>% mutate(VAZ = Violaxanthin + Antheraxanthin + Zeaxanthin) -> data_wide

# add VAZ/Chl
data_wide %>% mutate(VAZ.Chl = VAZ/Tot.Chl) -> data_wide

# add VAZ/Car
data_wide %>% mutate(VAZ.Car = VAZ/Tot.Car) -> data_wide

# add NVZ
data_wide %>% mutate(NVZ = Neoxanthin + Violaxanthin + Zeaxanthin) -> data_wide

# add NVZ/Car
data_wide %>% mutate(NVZ.Car = NVZ/Tot.Car) -> data_wide

# add Lutein/Chl
data_wide %>% mutate(Lut.Chl = Lutein/Tot.Chl) -> data_wide

# add Lutein/Car
data_wide %>% mutate(Lut.Car = Lutein/Tot.Car) -> data_wide

# add Lutein epoxide/Car
data_wide %>% mutate(Lut.epo.Car = Lutein.epoxide/Tot.Car) -> data_wide

# add a.Carotene/Car
data_wide %>% mutate(a.Car.Car = a.Carotene/Tot.Car) -> data_wide

# add b.Carotene/Car
data_wide %>% mutate(b.Car.Car = b.Carotene/Tot.Car) -> data_wide

# add Zeaxanthin/Car
data_wide %>% mutate(Zea.Car = Zeaxanthin/Tot.Car) -> data_wide

# add Violaxanthin/Car
data_wide %>% mutate(Vio.Car = Violaxanthin/Tot.Car) -> data_wide

# add Neoxanthin/Car
data_wide %>% mutate(Neo.Car = Neoxanthin/Tot.Car) -> data_wide

# add Antheraxanthin/Car
data_wide %>% mutate(Ant.Car = Antheraxanthin/Tot.Car) -> data_wide

# add DEPS
data_wide %>% mutate(DEPS = (0.5*Antheraxanthin + Zeaxanthin)/VAZ) -> data_wide
# convert DEPS from fraction to %
data_wide$DEPS <- data_wide$DEPS*100

# convert data back to long
data_long_calcs <- gather(data_wide, Pigment, FW.norm, a.Tocopherol:DEPS, factor_key=TRUE)


# convert to long format
data_long <- gather(samples_naomit, Pigment, ng, a.Tocopherol:Car.Chl31, factor_key=TRUE)

# normalize to fresh weight, ng/mg FW
data_long %>% mutate(., FW.norm = ng*Conversion) -> data_long

# replace Inf with NA (for 0/0 ratios)
is.na(data_long_calcs$FW.norm) <- sapply(data_long_calcs$FW.norm, is.infinite)

# replace NaN with NA (for 0/0 ratios)
is.na(data_long_calcs$FW.norm) <- sapply(data_long_calcs$FW.norm, is.nan)

# omit rows that are NA in FW.norm
data_long_calcs %>% drop_na(FW.norm) -> data_long_calcs 


# create table of just ipomoea lutein epoxide for correlation analysis
ipo_le <- data_long_calcs %>% filter(Species == "Ipomoea_nil") %>% filter(Pigment == "Lutein.epoxide")
# save a copy of this dataset
write.csv(ipo_le, file = "../output/stat_results/ipomoea_le.csv", row.names = F)

#### overall patterns stats ####
# is pigment composition different across subgenera (and tissues)?
# shapiro test for normality per pigment 
data_long_calcs[sample(nrow(data_long_calcs), 5000), ] %>% ungroup() %>% 
  rstatix::shapiro_test(FW.norm, data = .) -> shapiro_sub
shapiro_sub %>%  add_significance("p") -> shapiro_sub
shapiro_sub
# not normal, use kruskal-wallace

# kruskal wallace per pigment 
data_long_calcs %>% ungroup() %>% 
  rstatix::kruskal_test(FW.norm ~ Subgenus, data = .) -> kruskal_sub
kruskal_sub

# ajust p-values
kruskal_sub$p.adj <- p.adjust(kruskal_sub$p, method = "BH", n = length(kruskal_sub$p))
kruskal_sub$p.adj.signif <- stars.pval(kruskal_sub$p.adj)

# save a copy of this anova
write.csv(kruskal_sub, file = "../output/stat_results/kruskal_sub.csv", row.names = F)

# is pigment composition different across tissues (and subgenera)?
# shapiro test for normality per pigment 
data_long_calcs[sample(nrow(data_long_calcs), 5000), ] %>% ungroup() %>% 
  rstatix::shapiro_test(FW.norm, data = .) -> shapiro_tiss
shapiro_tiss %>%  add_significance("p") -> shapiro_tiss
shapiro_tiss
# not normal, use kruskal-wallace

# kruskal wallace tissues
data_long_calcs %>% ungroup() %>% 
  rstatix::kruskal_test(FW.norm ~ Tissue.code, data = .) -> kruskal_tiss
kruskal_tiss

# ajust p-values
kruskal_tiss$p.adj <- p.adjust(kruskal_tiss$p, method = "BH", n = length(kruskal_tiss$p))
kruskal_tiss$p.adj.signif <- stars.pval(kruskal_tiss$p.adj)

# save a copy of this anova
write.csv(kruskal_tiss, file = "../output/stat_results/kruskal_tiss.csv", row.names = F)


# is each pigment different across subgenera?
# shapiro test for normality per pigment 
data_long_calcs %>% ungroup() %>% 
  group_by(Pigment) %>% 
  rstatix::shapiro_test(FW.norm, data = .) -> shapiro_perpigment_sub
shapiro_perpigment_sub %>%  add_significance("p") -> shapiro_perpigment_sub
shapiro_perpigment_sub
# most are not normal, use kruskal-wallace

# kruskal wallace per pigment 
data_long_calcs %>% ungroup() %>% 
  group_by(Pigment) %>% 
  rstatix::kruskal_test(FW.norm ~ Subgenus, data = .) -> kruskal_Subgenus_perpigment
kruskal_Subgenus_perpigment

# ajust p-values
kruskal_Subgenus_perpigment$p.adj <- p.adjust(kruskal_Subgenus_perpigment$p, method = "BH", n = length(kruskal_Subgenus_perpigment$p))
kruskal_Subgenus_perpigment$p.adj.signif <- stars.pval(kruskal_Subgenus_perpigment$p.adj)

# save a copy of this anova
write.csv(kruskal_Subgenus_perpigment, file = "../output/stat_results/kruskal_Subgenus_perpigment.csv", row.names = F)


#### all tissues x subgenera comparison, per pigment ####
# within each metric, test for differences among genera and among tissues. achieved using the metric__subgenus
# add a column with full tissue names spelled out
data_long_calcs$Tissue.code_names <- data_long_calcs$Tissue.code
plyr::revalue(data_long_calcs$Tissue.code_names, c("sdlg" = "Seedling", "l" = "Leaf", "y" = "Young", "o" = "Old", "h" = "Haustorium", "f" = "Flower", "s" = "Seed")) -> data_long_calcs$Tissue.code_names

# kruskal wallace per pigment by Subgenus.Tissue
data_long_calcs$Subgenus.Tissue <- paste0(data_long_calcs$Subgenus, "__", data_long_calcs$Tissue.code_names)

# kruskal wallace in loop, per pigment among all subgenera-tissues
kruskal_big_list <- list()
pigment_list <- unique(data_long_calcs$Pigment)
for(pig in pigment_list) {
  data_loop_kruskal <- data_long_calcs %>% filter(Pigment == pig) 
  data_loop_kruskal$Tissue.code_names <- droplevels(data_loop_kruskal$Tissue.code_names)
  rstatix::kruskal_test(FW.norm ~ Subgenus.Tissue, data = data_loop_kruskal) -> kruskal_big_list[[pig]]
}

# all highly sig, move on to post hoc 
wilcox_big_list <- list()
pigment_list <- unique(data_long_calcs$Pigment)
for(pig in pigment_list) {
  data_loop_wilcox <- data_long_calcs %>% filter(Pigment == pig) 
  data_loop_wilcox$Tissue.code_names <- droplevels(data_loop_wilcox$Tissue.code_names)
  pairwise.wilcox.test(data_loop_wilcox$FW.norm, data_loop_wilcox$Subgenus.Tissue, p.adjust.method = "BH") -> wilcox_big_list[[pig]]
}

# export only plotted pigments
plotted_pigs <- c("Chl.a", "Chl.b", "Tot.Chl", "Chl.a.b", "VAZ", "Neoxanthin", "Lutein.epoxide", "Lutein", "a.Carotene", "b.Carotene", "Tot.Car", "NVZ.Car")
# format p-values, NAs, and column headers
excel_list <- list()
for(pig in plotted_pigs) {
  sheet_loop <- wilcox_big_list[[pig]][["p.value"]] 
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
  excel_list[[pig]] <- sheet_loop
}

# save pairwise wilcox test with a different sheet for each pigment
write.xlsx(excel_list[["Chl.a"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Chl.a", row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Chl.b"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Chl.b", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Tot.Chl"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Tot.Chl", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Chl.a.b"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Chl.a.b", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["VAZ"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="VAZ", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Neoxanthin"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Neoxanthin", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Lutein.epoxide"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Lutein.epoxide", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Lutein"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Lutein", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["a.Carotene"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="a.Carotene", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["b.Carotene"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="b.Carotene", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["Tot.Car"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="Tot.Car", append=TRUE, row.names=FALSE, col.names=FALSE)
write.xlsx(excel_list[["NVZ.Car"]], file="../output/stat_results/pairwise_allpigsxtissues.xlsx", sheetName="NVZ.Car", append=TRUE, row.names=FALSE, col.names=FALSE)



#### stats to be used for all pigment plots ####

# some pigs were causing an error because it is 0 across all tissues in this species
# i identified them manually with:
data_long_calcs %>% group_by(Pigment, Subgenus) %>% dplyr::summarize(mean = mean(FW.norm)) %>% filter(mean == 0) -> test_means
# # A tibble: 5 × 3
# # Groups:   Pigment [3]
# Pigment      Subgenus     mean
# <fct>        <chr>       <dbl>
#   1 a.Tocopherol Monogynella     0
# 2 Neoxanthin   C_purpurata     0
# 3 Neoxanthin   Cuscuta         0
# 4 Neo.Car      C_purpurata     0
# 5 Neo.Car      Cuscuta         0

# some pigs were causing an error because in that pigment x subgenus combo there was only 1 tissue
# i identified them manually with:
data_long_calcs %>% na.omit() %>% group_by(Pigment, Subgenus, Tissue.code) %>% summarize(mean = mean(FW.norm)) %>% dplyr::summarize(n = n()) %>% filter(n == 1) -> test_tissues
# Chl.a.b	C_purpurata	1

# shapiro test for normality per pigment per subgenus
# removing those combos with all 0s
data_long_calcs %>% ungroup() %>% 
  filter(!(Pigment == "a.Tocopherol" & Subgenus == "Monogynella")) %>% 
  filter(!(Pigment == "Neoxanthin" & Subgenus == "C_purpurata")) %>% 
  filter(!(Pigment == "Neoxanthin" & Subgenus == "Cuscuta")) %>% 
  filter(!(Pigment == "Neo.Car" & Subgenus == "C_purpurata")) %>% 
  filter(!(Pigment == "Neo.Car" & Subgenus == "Cuscuta")) %>% 
  filter(!(Pigment == "Chl.a.b" & Subgenus == "C_purpurata")) %>% 
  group_by(Pigment, Subgenus) %>% 
  rstatix::shapiro_test(FW.norm, data = .) -> shapiro_pig_sub
shapiro_pig_sub %>%   add_significance("p") -> shapiro_pig_sub
shapiro_pig_sub
# most are not normal, use kruskal-wallace

# kruskal wallace per pigment per subgenus
data_long_calcs %>% ungroup() %>% 
  filter(!(Pigment == "a.Tocopherol" & Subgenus == "Monogynella")) %>% 
  filter(!(Pigment == "Neoxanthin" & Subgenus == "C_purpurata")) %>% 
  filter(!(Pigment == "Neoxanthin" & Subgenus == "Cuscuta")) %>% 
  filter(!(Pigment == "Neo.Car" & Subgenus == "C_purpurata")) %>% 
  filter(!(Pigment == "Neo.Car" & Subgenus == "Cuscuta")) %>% 
  filter(!(Pigment == "Chl.a.b" & Subgenus == "C_purpurata")) %>% 
  group_by(Pigment, Subgenus) %>% 
  rstatix::kruskal_test(FW.norm ~ Tissue.code, data = .) -> kruskal_Tissue_pig_sub
kruskal_Tissue_pig_sub
 
# ajust p-values
kruskal_Tissue_pig_sub$p.adj <- p.adjust(kruskal_Tissue_pig_sub$p, method = "BH", n = length(kruskal_Tissue_pig_sub$p))
kruskal_Tissue_pig_sub$p.adj.signif <- stars.pval(kruskal_Tissue_pig_sub$p.adj)

# save a copy of this anova
write.csv(kruskal_Tissue_pig_sub, file = "../output/stat_results/kruskal_Tissue_pig_sub.csv", row.names = F)


# filter for signficant pigments x main effects for posthoc comparisons, but round p-values first
kruskal_Tissue_pig_sub$p.adj <- kruskal_Tissue_pig_sub$p.adj %>% p_round(digits = 3)
kruskal_Tissue_pig_sub_sig <- as.data.frame(kruskal_Tissue_pig_sub) %>% dplyr::filter(p.adj <= 0.05)



# perform posthoc pairwise comparisons on the sig anovas
# first create a vector on which to filter data_long_calcs by
kruskal_Tissue_pig_sub_sig_vector <- paste0(kruskal_Tissue_pig_sub_sig$Subgenus, "__", kruskal_Tissue_pig_sub_sig$Pigment)
# create a combined column in data_long_calcs too
data_long_calcs$Subgenus.Pigment <- paste0(data_long_calcs$Subgenus, "__", data_long_calcs$Pigment)
# filter data long calcs by the Subgenus.Pigment column, based on vector
filter(data_long_calcs, Subgenus.Pigment %in% kruskal_Tissue_pig_sub_sig_vector) -> only_sig_pig

####### #######
# perform posthoc on these
# filter for one pigment x subgenus then do a pairwise wilcox on those
# after getting it working for one, set up loop to do for all combos (for those that had a positive kruskal-wallace at least)
# add these to individual plots that then get arranged with patchwork package

only_sig_pig$Pigment <- droplevels(only_sig_pig$Pigment)
subgenuspigment_list <- unique(only_sig_pig$Subgenus.Pigment)
only_sig_pig$log.FW.norm <- log(only_sig_pig$FW.norm)

wilcox_list <- list()
for(subpig in subgenuspigment_list) {
data_loop_wilcox <- only_sig_pig %>% filter(Subgenus.Pigment == subpig) 
data_loop_wilcox$Tissue.code <- droplevels(data_loop_wilcox$Tissue.code)
data_loop_wilcox$Pigment <- droplevels(data_loop_wilcox$Pigment)
pairwise.wilcox.test(data_loop_wilcox$FW.norm, data_loop_wilcox$Tissue.code, p.adjust.method = "BH") -> wilcox_list[[subpig]]
}

# save this list for later use 
saveRDS(wilcox_list, file="../output/stat_results/wilcox_list.RData")




# create a dataframe using the kruskal wallace results
dat_text_plot_kruskal <- kruskal_Tissue_pig_sub[, c("Subgenus", "Pigment", "method", "p.adj")]
# add equal sign to p.adj
dat_text_plot_kruskal$p.adj_eq <- paste0(" = ", dat_text_plot_kruskal$p.adj)
dat_text_plot_kruskal$label <- paste0('paste(italic("P"),"', dat_text_plot_kruskal$p.adj_eq, "\")") 

dat_text_plot_kruskal <- dat_text_plot_kruskal %>% select(Subgenus, Pigment, label)

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



#### summary  ####
# summary stats per Accession.No (summarizes those that have replicates), which is also per species: mean, n, and std dev  

summary_accession <- data_long_calcs %>%
  group_by(Subgenus, Species, Accession.No, Tissue.code, Pigment)  %>%
  dplyr::summarize(Mean = mean(FW.norm, na.rm = TRUE), n = sum(!is.na(FW.norm)), sd = sd(FW.norm, na.rm = TRUE))
# replace NaN with NA (for 0/0 ratios)
is.na(summary_accession$Mean) <- sapply(summary_accession$Mean, is.nan)


# write to csv
write_csv(summary_accession, file = "../output/stat_results/pigments_species_summary.csv")


# summary stats per subgenus x tissue x pigment: mean, n, and std dev  
summary_subgenus <- data_long_calcs %>%
  group_by(Subgenus, Tissue.code, Pigment)  %>%
  dplyr::summarize(Mean = mean(FW.norm, na.rm = TRUE), n = sum(!is.na(FW.norm)), sd = sd(FW.norm, na.rm = TRUE))

# replace NaN with NA (for 0/0 ratios)
is.na(summary_subgenus$Mean) <- sapply(summary_subgenus$Mean, is.nan)

# write to csv
write_csv(summary_subgenus, file = "../output/stat_results/pigments_subgenus_summary.csv")


# convert FW.norm to log(FW.norm) for plotting (when necesary), when > 0
data_long_calcs$logFW.norm <- log(data_long_calcs$FW.norm + 0.0000001)
data_long_calcs$logFW.norm <- ifelse(data_long_calcs$logFW.norm < 0, 0, data_long_calcs$logFW.norm)

# write this to csv for plotting
write_csv(data_long_calcs, file = "../output/stat_results/data_long_calcs_for_plots.csv")





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

data_long_calcs_Grammica_plot <- dplyr::filter(data_long_calcs, Species %in% Grammica )

# identify species with n = 1
summary_accession_Grammica <- dplyr::filter(summary_accession, Species %in% Grammica )
average_ns <- summary_accession_Grammica %>% group_by(Subgenus, Species, Accession.No, Tissue.code) %>% dplyr::summarize(average_n = mean(n))

summary_accession_n1 <- dplyr::filter(average_ns, average_n == 1 )
summary_accession_n1 %>% select(Subgenus, Species, Accession.No, Tissue.code) -> summary_accession_n1

# these are the specimens being removed for n = 1
Grammica_n1 <- unique(summary_accession_n1)
# add variable
Grammica_n1$drop <- paste0(Grammica_n1$Accession.No, Grammica_n1$Tissue.code)

# remove species with n = 1 from plot data
data_long_calcs_Grammica_plot$drop <- paste0(data_long_calcs_Grammica_plot$Accession.No, data_long_calcs_Grammica_plot$Tissue.code)

data_long_calcs_Grammica_plot %>% dplyr::filter(!drop %in% Grammica_n1$drop) -> data_long_calcs_Grammica_plot


data_long_calcs_Grammica_plot$Species <- factor(data_long_calcs_Grammica_plot$Species,
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

# some pigs were causing an error because it is 0 across all tissues in this species
# i identified them manually with:
data_long_calcs_Grammica_plot %>% group_by(Pigment, Species) %>% dplyr::summarize(mean = mean(FW.norm)) %>% dplyr::filter(mean == 0) -> test_means_Grammica
# # A tibble: 7 x 3
# # Groups:   Pigment [1]
# Pigment    Species          mean
# <fct>      <fct>           <dbl>
# 1 Neoxanthin C_sandwichiana      0
# 2 Neoxanthin C_californica       0
# 3 Neoxanthin C_compacta          0
# 4 Neoxanthin C_denticulata       0
# 5 Neoxanthin C_tasmanica         0
# 6 Neoxanthin C_costaricensis     0
# 7 Neoxanthin C_indecora          0
# 8 Neo.Car    C_sandwichiana      0
# 9 Neo.Car    C_californica       0
# 10 Neo.Car    C_compacta          0
# 11 Neo.Car    C_denticulata       0
# 12 Neo.Car    C_tasmanica         0
# 13 Neo.Car    C_costaricensis     0
# 14 Neo.Car    C_indecora          0

# shapiro test for normality per pigment per species in Grammica only
data_long_calcs_Grammica_plot %>% ungroup() %>% 
  dplyr::filter(., Subgenus == "Grammica") %>%
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_sandwichiana")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_californica")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_compacta")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_denticulata")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_tasmanica")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_costaricensis")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_indecora")) %>%
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_sandwichiana")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_californica")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_compacta")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_denticulata")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_tasmanica")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_costaricensis")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_indecora")) %>%
  group_by(Pigment, Species) %>% 
  rstatix::shapiro_test(FW.norm, data = .) -> shapiro_pig_spe_Grammica

shapiro_pig_spe_Grammica %>%   add_significance("p") -> shapiro_pig_spe_Grammica
shapiro_pig_spe_Grammica
# many are not normal, use kruskal-wallace

# kruskal wallace per pigment per species in Grammica only
data_long_calcs_Grammica_plot %>% ungroup() %>%
  dplyr::filter(., Subgenus == "Grammica") %>%
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_sandwichiana")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_californica")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_compacta")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_denticulata")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_tasmanica")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_costaricensis")) %>% 
  dplyr::filter(!(Pigment == "Neoxanthin" & Species == "C_indecora")) %>%
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_sandwichiana")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_californica")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_compacta")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_denticulata")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_tasmanica")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_costaricensis")) %>% 
  dplyr::filter(!(Pigment == "Neo.Car" & Species == "C_indecora")) %>%
  group_by(Pigment, Species) %>% 
  rstatix::kruskal_test(FW.norm ~ Tissue.code, data = .) -> kruskal_Tissue_pig_spe_Grammica
kruskal_Tissue_pig_spe_Grammica

# ajust p-values
kruskal_Tissue_pig_spe_Grammica$p.adj <- p.adjust(kruskal_Tissue_pig_spe_Grammica$p, method = "BH", n = length(kruskal_Tissue_pig_spe_Grammica$p))
kruskal_Tissue_pig_spe_Grammica$p.adj.signif <- stars.pval(kruskal_Tissue_pig_spe_Grammica$p.adj)

# save a copy of this anova
write.csv(kruskal_Tissue_pig_spe_Grammica, file = "../output/stat_results/kruskal_Tissue_pig_spe_Grammica.csv", row.names = F)


# dplyr::filter for signficant pigments x main effects for posthoc comparisons, but round p-values first
kruskal_Tissue_pig_spe_Grammica$p.adj <- kruskal_Tissue_pig_spe_Grammica$p.adj %>% p_round(digits = 3)
kruskal_Tissue_pig_spe_sig_Grammica <- as.data.frame(kruskal_Tissue_pig_spe_Grammica) %>% dplyr::filter(p.adj <= 0.05)



# perform posthoc pairwise comparisons on the sig anovas
# first create a vector on which to dplyr::filter data_long_calcs_Grammica_plot by
kruskal_Tissue_pig_spe_sig_vector_Grammica <- paste0(kruskal_Tissue_pig_spe_sig_Grammica$Species, "__", kruskal_Tissue_pig_spe_sig_Grammica$Pigment)
# create a combined column in data_long_calcs_Grammica_plot too
data_long_calcs_Grammica_plot$Species.Pigment <- paste0(data_long_calcs_Grammica_plot$Species, "__", data_long_calcs_Grammica_plot$Pigment)
# dplyr::filter data long calcs by the Species.Pigment column, based on vector
filter(data_long_calcs_Grammica_plot, Species.Pigment %in% kruskal_Tissue_pig_spe_sig_vector_Grammica) -> only_sig_pig_Grammica

# write this to csv for plotting
write_csv(data_long_calcs_Grammica_plot, file = "../output/stat_results/data_long_calcs_for_Grammica_plots.csv")



####### #######
# perform posthoc on these
# dplyr::filter for one pigment x subgenus then do a pairwise wilcox on those
# after getting it working for one, set up loop to do for all combos (for those that had a positive kruskal-wallace at least)
# add these to individual plots that then get arranged with patchwork package

only_sig_pig_Grammica$Pigment <- droplevels(only_sig_pig_Grammica$Pigment)
specispigment_list <- unique(only_sig_pig_Grammica$Species.Pigment)
only_sig_pig_Grammica$log.FW.norm <- log(only_sig_pig_Grammica$FW.norm)

wilcox_list_Grammica <- list()
for(spepig in specispigment_list) {
  data_loop_wilcox <- only_sig_pig_Grammica %>% dplyr::filter(Species.Pigment == spepig) 
  data_loop_wilcox$Tissue.code <- droplevels(data_loop_wilcox$Tissue.code)
  data_loop_wilcox$Pigment <- droplevels(data_loop_wilcox$Pigment)
  pairwise.wilcox.test(data_loop_wilcox$FW.norm, data_loop_wilcox$Tissue.code, p.adjust.method = "BH") -> wilcox_list_Grammica[[spepig]]
}

# save this list for later use 
saveRDS(wilcox_list_Grammica, file="../output/stat_results/wilcox_list_Grammica.RData")




# will need to add kruskal wallace results to each box of the subgenus boxplot from kruskal_Tissue_pig_spe (including not significant global result)
# edit scientific notation 
kruskal_Tissue_pig_spe_Grammica$p.adj <- kruskal_Tissue_pig_spe_Grammica$p.adj %>% p_round(digits = 3)

# # will need to add wilcox test asterisks to each tissue in each box of the subgenus boxplot from posthoc_pig_spe (including ns?)


# create a dataframe using the kruskal wallace results
dat_text_plot_kruskal_Grammica <- kruskal_Tissue_pig_spe_Grammica[, c("Species", "Pigment", "method", "p.adj")]
# add equal sign to p.adj
dat_text_plot_kruskal_Grammica$p.adj_eq <- paste0(" = ", dat_text_plot_kruskal_Grammica$p.adj)
dat_text_plot_kruskal_Grammica$label <- paste0('paste(italic("P"),"', dat_text_plot_kruskal_Grammica$p.adj_eq, "\")") 

dat_text_plot_kruskal_Grammica <- dat_text_plot_kruskal_Grammica %>% select(Species, Pigment, label)

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

