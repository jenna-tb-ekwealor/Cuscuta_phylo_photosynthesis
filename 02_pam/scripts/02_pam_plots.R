library(tidyverse)
library(ggplot2)
library(dplyr)
library(colorBlindness)
library(patchwork)
library(grid)
library(data.table)
# install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(phyloseq)
# install.packages("devtools")
# devtools::install_github("jfq3/QsRutils", build_vignettes = TRUE)
library(QsRutils)
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

# add colors for life stages, most picked from colorBlindness::paletteMartin
# displayAllColors(paletteMartin)

leaf<-"#24ff24" #Harlequin
young<-"orange" 
old<-"#924900" #Brown
haustorium<- "red"
seedling<-"#009292" #PersianGreen
flower<-"#ff6db6" #HotPink
fruit<-"#ffff6d" #LaserLemon
seed<- "#6db6ff" #Malibu
                
# load data
data_no_outliers <- read.csv(file = "../output/stat_results/data_no_outliers_for_plots.csv", stringsAsFactors = T)
dat_text_plot_kruskal <- read.csv(file = "../output/stat_results/dat_text_plot_kruskal.csv")
data_no_outliers_Grammica_plot <- read.csv(file = "../output/stat_results/data_no_outliers_for_Grammica_plots.csv", stringsAsFactors = T)
dat_text_plot_kruskal_Grammica <- read.csv(file = "../output/stat_results/dat_text_plot_kruskal_Grammica.csv")
summary_accession <- read.csv(file = "../output/stat_results/fluorescence_species_summary.csv", stringsAsFactors = T)
wilcox_list <- readRDS("../output/stat_results/wilcox_list.RData")
wilcox_list_Grammica <- readRDS("../output/stat_results/wilcox_list_Grammica.RData")
              

#### FLUORESCENCE PLOTS ####

#### Fv.Fm Ipomoea_nil ####
data_ipomoea <-  dplyr::filter(data_no_outliers, Subgenus == "Ipomoea_nil")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_ipomoea$Tissue.edit <- factor(data_ipomoea$Tissue.edit, levels = c("l", "y", "o", "f", "s"))

data_ipomoea_Fv.Fm <- dplyr::filter(data_ipomoea, Metric == "Fv.Fm")
# calculate minimum Fv.Fm for apporpriate bottom limit of Fv.Fm y axes
min_Fv.Fm <-as.numeric( data_no_outliers %>% dplyr::filter(., Metric == "Fv.Fm") %>% dplyr::summarize(., min(Value)))
min_Fv.Fm

Fv.Fm_ipomoea_boxplot <- ggplot(data_ipomoea_Fv.Fm, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c( "l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("l" = "L", "y" = "Y", "o" = "O", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5),        
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(min_Fv.Fm, 1)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "Fv.Fm" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

Fv.Fm_ipomoea_boxplot 

# use base R boxplot to get the coordinates of the boxes
box.rslt_Fv.Fm_ipomoea <- with(data_ipomoea_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_ipomoea)
boxplot_positions_Fv.Fm_ipomoea <- as.data.frame(box.rslt_Fv.Fm_ipomoea$stats)

# what are these column tissue codes?
tissues_Fv.Fm_ipomoea <- levels(data_ipomoea_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_ipomoea) <- tissues_Fv.Fm_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_ipomoea <- boxplot_positions_Fv.Fm_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_ipomoea, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_ipomoea %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_ipomoea
# now join these positions to cbd
left_join(cbd_Fv.Fm_ipomoea, top_positions_Fv.Fm_ipomoea, by = "Tissue.edit") -> cbd_Fv.Fm_ipomoea

# calculate how much to nudge
data_ipomoea_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_ipomoea
cbd_Fv.Fm_ipomoea$nudged <- max_Fv.Fm_ipomoea$max * 1.05


# add CLDs to plot
Fv.Fm_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Fv.Fm_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_ipomoea_boxplot

Fv.Fm_ipomoea_boxplot





#### φPSII Ipomoea_nil ####
data_ipomoea_φPSII <- dplyr::filter(data_ipomoea, Metric == "φPSII")
# calculate minimum φPSII for apporpriate bottom limit of φPSII y axes
min_φPSII <-as.numeric( data_no_outliers %>% dplyr::filter(., Metric == "φPSII") %>% dplyr::summarize(., min(Value)))
min_φPSII

φPSII_ipomoea_boxplot <- ggplot(data_ipomoea_φPSII, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c( "l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("l" = "L", "y" = "Y", "o" = "O", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(min_φPSII, 1)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "φPSII" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_φPSII_ipomoea <- with(data_ipomoea_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_ipomoea)
boxplot_positions_φPSII_ipomoea <- as.data.frame(box.rslt_φPSII_ipomoea$stats)

# what are these column tissue codes?
tissues_φPSII_ipomoea <- levels(data_ipomoea_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_ipomoea) <- tissues_φPSII_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_ipomoea <- boxplot_positions_φPSII_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_ipomoea, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_ipomoea %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_ipomoea
# now join these positions to cbd
left_join(cbd_φPSII_ipomoea, top_positions_φPSII_ipomoea, by = "Tissue.edit") -> cbd_φPSII_ipomoea

# calculate how much to nudge
data_ipomoea_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_ipomoea
cbd_φPSII_ipomoea$nudged <- max_φPSII_ipomoea$max* 1.05


# add CLDs to plot
φPSII_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_φPSII_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_ipomoea_boxplot


φPSII_ipomoea_boxplot


#### NPQ Ipomoea_nil plot with x axis ####
data_ipomoea_NPQ <- dplyr::filter(data_ipomoea, Metric == "NPQ")
# calculate minimum NPQ for apporpriate bottom limit of NPQ y axes
max_NPQ <-as.numeric( data_no_outliers %>% dplyr::filter(., Metric == "NPQ") %>% dplyr::summarize(., max(Value)))
max_NPQ


NPQ_ipomoea_boxplot <- ggplot(data_ipomoea_NPQ, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c( "l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("l" = "L", "y" = "Y", "o" = "O", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 4))+ 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "NPQ" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_NPQ_ipomoea <- with(data_ipomoea_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_ipomoea)
boxplot_positions_NPQ_ipomoea <- as.data.frame(box.rslt_NPQ_ipomoea$stats)

# what are these column tissue codes?
tissues_NPQ_ipomoea <- levels(data_ipomoea_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_ipomoea) <- tissues_NPQ_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_ipomoea <- boxplot_positions_NPQ_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_ipomoea, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_ipomoea %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_ipomoea
# now join these positions to cbd
left_join(cbd_NPQ_ipomoea, top_positions_NPQ_ipomoea, by = "Tissue.edit") -> cbd_NPQ_ipomoea

# calculate how much to nudge
data_ipomoea_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_ipomoea
cbd_NPQ_ipomoea$nudged <- max_NPQ_ipomoea$max * 1.05


# add CLDs to plot
NPQ_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NPQ_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_ipomoea_boxplot


NPQ_ipomoea_boxplot


#### Fv.Fm loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# filter for Fv.Fm
data_Fv.Fm_plots_cuscutasub <- dplyr::filter(data_no_outliers, Metric == "Fv.Fm")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Fv.Fm_plots_cuscutasub$Tissue.edit <- factor(data_Fv.Fm_plots_cuscutasub$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Fv.Fm = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Fv.Fm_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "left", limits = c(min_Fv.Fm, 1))
  plot_list_Cuscuta_Fv.Fm[[sub]] = p
  
}


#### Fv.Fm Monogynella ####
plot_list_Cuscuta_Fv.Fm[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "Fv.Fm" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Fv.Fm_Monogynella <- with(data_Monogynella_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_Monogynella)
boxplot_positions_Fv.Fm_Monogynella <- as.data.frame(box.rslt_Fv.Fm_Monogynella$stats)

# what are these column tissue codes?
tissues_Fv.Fm_Monogynella <- levels(data_Monogynella_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_Monogynella) <- tissues_Fv.Fm_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_Monogynella <- boxplot_positions_Fv.Fm_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_Monogynella, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_Monogynella %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_Monogynella
# now join these positions to cbd
left_join(cbd_Fv.Fm_Monogynella, top_positions_Fv.Fm_Monogynella, by = "Tissue.edit") -> cbd_Fv.Fm_Monogynella

# calculate how much to nudge
data_Monogynella_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_Monogynella
cbd_Fv.Fm_Monogynella$nudged <- max_Fv.Fm_Monogynella$max * 1.05


# add CLDs to plot
Fv.Fm_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Fv.Fm_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_Monogynella_boxplot


Fv.Fm_Monogynella_boxplot


#### Fv.Fm Cuscuta ####
plot_list_Cuscuta_Fv.Fm[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "Fv.Fm" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_Cuscuta_boxplot


# KW not significant so no post hoc or CLDs
Fv.Fm_Cuscuta_boxplot






#### Fv.Fm Grammica ####
plot_list_Cuscuta_Fv.Fm[["Grammica"]]  + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "Fv.Fm" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Fv.Fm_Grammica <- with(data_Grammica_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_Grammica)
boxplot_positions_Fv.Fm_Grammica <- as.data.frame(box.rslt_Fv.Fm_Grammica$stats)

# what are these column tissue codes?
tissues_Fv.Fm_Grammica <- levels(data_Grammica_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_Grammica) <- tissues_Fv.Fm_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_Grammica <- boxplot_positions_Fv.Fm_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_Grammica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_Grammica %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_Grammica
# now join these positions to cbd
left_join(cbd_Fv.Fm_Grammica, top_positions_Fv.Fm_Grammica, by = "Tissue.edit") -> cbd_Fv.Fm_Grammica

# calculate how much to nudge
data_Grammica_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_Grammica
cbd_Fv.Fm_Grammica$nudged <- max_Fv.Fm_Grammica$max * 1.05

# add CLDs to plot
Fv.Fm_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Fv.Fm_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_Grammica_boxplot


Fv.Fm_Grammica_boxplot


#### Fv.Fm C_purpurata alone ####

data_C_purpurata_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_cuscutasub, Subgenus == "C_purpurata")

Fv.Fm_C_purpurata_boxplot <- ggplot(data_C_purpurata_Fv.Fm, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(min_Fv.Fm, 1)) 

Fv.Fm_C_purpurata_boxplot

## not sig


#### φPSII loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# filter for φPSII
data_φPSII_plots_cuscutasub <- dplyr::filter(data_no_outliers, Metric == "φPSII")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_φPSII_plots_cuscutasub$Tissue.edit <- factor(data_φPSII_plots_cuscutasub$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_φPSII = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_φPSII_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(min_φPSII, 1))
  plot_list_Cuscuta_φPSII[[sub]] = p
  
}

#### φPSII. Monogynella ####
plot_list_Cuscuta_φPSII[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "φPSII" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_φPSII <- dplyr::filter(data_φPSII_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_φPSII_Monogynella <- with(data_Monogynella_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_Monogynella)
boxplot_positions_φPSII_Monogynella <- as.data.frame(box.rslt_φPSII_Monogynella$stats)

# what are these column tissue codes?
tissues_φPSII_Monogynella <- levels(data_Monogynella_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_Monogynella) <- tissues_φPSII_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_Monogynella <- boxplot_positions_φPSII_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_Monogynella, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_Monogynella %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_Monogynella
# now join these positions to cbd
left_join(cbd_φPSII_Monogynella, top_positions_φPSII_Monogynella, by = "Tissue.edit") -> cbd_φPSII_Monogynella

# calculate how much to nudge
data_Monogynella_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_Monogynella
cbd_φPSII_Monogynella$nudged <- max_φPSII_Monogynella$max * 1.05


# add CLDs to plot
φPSII_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_φPSII_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_Monogynella_boxplot


φPSII_Monogynella_boxplot


#### φPSII Cuscuta ####
plot_list_Cuscuta_φPSII[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "φPSII" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_φPSII <- dplyr::filter(data_φPSII_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_φPSII_Cuscuta <- with(data_Cuscuta_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_Cuscuta)
boxplot_positions_φPSII_Cuscuta <- as.data.frame(box.rslt_φPSII_Cuscuta$stats)

# what are these column tissue codes?
tissues_φPSII_Cuscuta <- levels(data_Cuscuta_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_Cuscuta) <- tissues_φPSII_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_Cuscuta <- boxplot_positions_φPSII_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_Cuscuta, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_Cuscuta %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_Cuscuta
# now join these positions to cbd
left_join(cbd_φPSII_Cuscuta, top_positions_φPSII_Cuscuta, by = "Tissue.edit") -> cbd_φPSII_Cuscuta

# calculate how much to nudge
data_Cuscuta_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_Cuscuta
cbd_φPSII_Cuscuta$nudged <- max_φPSII_Cuscuta$max * 1.05


# add CLDs to plot
φPSII_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_φPSII_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_Cuscuta_boxplot


φPSII_Cuscuta_boxplot



#### φPSII Grammica ####
plot_list_Cuscuta_φPSII[["Grammica"]]  + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "φPSII" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_Grammica_boxplot

φPSII_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_φPSII <- dplyr::filter(data_φPSII_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_φPSII_Grammica <- with(data_Grammica_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_Grammica)
boxplot_positions_φPSII_Grammica <- as.data.frame(box.rslt_φPSII_Grammica$stats)

# what are these column tissue codes?
tissues_φPSII_Grammica <- levels(data_Grammica_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_Grammica) <- tissues_φPSII_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_Grammica <- boxplot_positions_φPSII_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_Grammica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_Grammica %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_Grammica
# now join these positions to cbd
left_join(cbd_φPSII_Grammica, top_positions_φPSII_Grammica, by = "Tissue.edit") -> cbd_φPSII_Grammica

# calculate how much to nudge
data_Grammica_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_Grammica
cbd_φPSII_Grammica$nudged <- max_φPSII_Grammica$max * 1.05

# add CLDs to plot
φPSII_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_φPSII_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_Grammica_boxplot


φPSII_Grammica_boxplot




#### φPSII C. purpurata alone ####

data_C_purpurata_φPSII <- dplyr::filter(data_φPSII_plots_cuscutasub, Subgenus == "C_purpurata")

φPSII_C_purpurata_boxplot <- ggplot(data_C_purpurata_φPSII, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(min_φPSII, 1)) 




#### NPQ loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# filter for NPQ
data_NPQ_plots_cuscutasub <- dplyr::filter(data_no_outliers, Metric == "NPQ")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_NPQ_plots_cuscutasub$Tissue.edit <- factor(data_NPQ_plots_cuscutasub$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_NPQ = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_NPQ_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 4))
  plot_list_Cuscuta_NPQ[[sub]] = p
  
}

#### NPQ Monogynella ####
plot_list_Cuscuta_NPQ[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "NPQ" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_NPQ <- dplyr::filter(data_NPQ_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_NPQ_Monogynella <- with(data_Monogynella_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_Monogynella)
boxplot_positions_NPQ_Monogynella <- as.data.frame(box.rslt_NPQ_Monogynella$stats)

# what are these column tissue codes?
tissues_NPQ_Monogynella <- levels(data_Monogynella_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_Monogynella) <- tissues_NPQ_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_Monogynella <- boxplot_positions_NPQ_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_Monogynella, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_Monogynella %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_Monogynella
# now join these positions to cbd
left_join(cbd_NPQ_Monogynella, top_positions_NPQ_Monogynella, by = "Tissue.edit") -> cbd_NPQ_Monogynella

# calculate how much to nudge
data_Monogynella_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_Monogynella
cbd_NPQ_Monogynella$nudged <- (max_NPQ_Monogynella$max + 0.00) * 1.05

# add CLDs to plot
NPQ_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NPQ_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_Monogynella_boxplot


NPQ_Monogynella_boxplot



#### NPQ Cuscuta ####
plot_list_Cuscuta_NPQ[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "NPQ" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_NPQ <- dplyr::filter(data_NPQ_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_NPQ_Cuscuta <- with(data_Cuscuta_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_Cuscuta)
boxplot_positions_NPQ_Cuscuta <- as.data.frame(box.rslt_NPQ_Cuscuta$stats)

# what are these column tissue codes?
tissues_NPQ_Cuscuta <- levels(data_Cuscuta_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_Cuscuta) <- tissues_NPQ_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_Cuscuta <- boxplot_positions_NPQ_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_Cuscuta, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_Cuscuta %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_Cuscuta
# now join these positions to cbd
left_join(cbd_NPQ_Cuscuta, top_positions_NPQ_Cuscuta, by = "Tissue.edit") -> cbd_NPQ_Cuscuta

# calculate how much to nudge
data_Cuscuta_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_Cuscuta
cbd_NPQ_Cuscuta$nudged <- (max_NPQ_Cuscuta$max + 0.00) * 1.05

# add CLDs to plot
NPQ_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NPQ_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_Cuscuta_boxplot


NPQ_Cuscuta_boxplot



#### NPQ Grammica ####
plot_list_Cuscuta_NPQ[["Grammica"]]  + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Metric == "NPQ" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_Grammica_boxplot 

# use base R boxplot to get the coordinates of the boxes
data_Grammica_NPQ <- dplyr::filter(data_NPQ_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_NPQ_Grammica <- with(data_Grammica_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_Grammica)
boxplot_positions_NPQ_Grammica <- as.data.frame(box.rslt_NPQ_Grammica$stats)

# what are these column tissue codes?
tissues_NPQ_Grammica <- levels(data_Grammica_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_Grammica) <- tissues_NPQ_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_Grammica <- boxplot_positions_NPQ_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_Grammica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_Grammica %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_Grammica
# now join these positions to cbd
left_join(cbd_NPQ_Grammica, top_positions_NPQ_Grammica, by = "Tissue.edit") -> cbd_NPQ_Grammica

# calculate how much to nudge
data_Grammica_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_Grammica
cbd_NPQ_Grammica$nudged <- (max_NPQ_Grammica$max + 0.00) * 1.05

# add CLDs to plot
NPQ_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NPQ_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_Grammica_boxplot


NPQ_Grammica_boxplot



#### NPQ C_purpurata alone ####
data_C_purpurata_NPQ <- dplyr::filter(data_NPQ_plots_cuscutasub, Subgenus == "C_purpurata")

NPQ_C_purpurata_boxplot <- ggplot(data_C_purpurata_NPQ, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 4))


#### add "absent" to those with mean Value = 0 in loop ####
Subgenus_list_Subgenus <- unique(data_no_outliers$Subgenus)
Metric_list_Subgenus <- unique(data_no_outliers$Metric)
label <- "absent"
absent <- data.frame(label)

Subgenus_fluorescence_absent_list <- list()

for (subgenus in Subgenus_list_Subgenus) {
  for (metric in Metric_list_Subgenus) {
    data_loop <- data_no_outliers %>% dplyr::filter(Subgenus == subgenus & Metric == metric)
    if (nrow(data_loop) > 0) { if (mean(data_loop$Value) == 0) {
      print("absent") -> Subgenus_fluorescence_absent_list[[metric]][[subgenus]]
    }}}}

# Subgenus_fluorescence_absent_list
# $Fv.Fm
# $Fv.Fm$C_purpurata
# [1] "absent"
# 
# 
# $φPSII
# $φPSII$C_purpurata
# [1] "absent"
# 
# 
# $ΦNPQ
# $ΦNPQ$C_purpurata
# [1] "absent"
# 
# 
# $NPQ
# $NPQ$C_purpurata
# [1] "absent"


# add the following geom_text to the above plots
# + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) 


#### COMBINED (faceted) fluorescence plot: subgenus ####
wrap_elements(gridtext::richtext_grob('*F*<sub>v</sub>/*F*<sub>m</sub>', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Fv.Fm_ipomoea_boxplot + ggtitle('Ipomoea nil') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Fv.Fm_Monogynella_boxplot + ggtitle('Monogynella') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Fv.Fm_Cuscuta_boxplot + ggtitle('Cuscuta') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Fv.Fm_Grammica_boxplot + ggtitle('Grammica') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) +  
  Fv.Fm_C_purpurata_boxplot + ggtitle('C. purpurata') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  wrap_elements(gridtext::richtext_grob('φPSII', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  φPSII_ipomoea_boxplot +
  φPSII_Monogynella_boxplot +  
  φPSII_Cuscuta_boxplot + 
  φPSII_Grammica_boxplot + 
  φPSII_C_purpurata_boxplot + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) + 
  wrap_elements(gridtext::richtext_grob('NPQ', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  NPQ_ipomoea_boxplot + 
  NPQ_Monogynella_boxplot +  
  NPQ_Cuscuta_boxplot + 
  NPQ_Grammica_boxplot + 
  NPQ_C_purpurata_boxplot + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) + plot_layout(nrow = 3, byrow = T) -> fluorescence_boxplot_new

fluorescence_boxplot_new 

pdf("../output/stat_results/boxplots/fluorescence_boxplot_new.pdf", width=7,height=7) 
fluorescence_boxplot_new
dev.off()








#### Grammica only FLUORESCENCE PLOTS  ####

#### Fv.Fm C_australis ####
data_C_australis <-  dplyr::filter(data_no_outliers_Grammica_plot, Species == "C_australis")
data_C_australis_Fv.Fm <- dplyr::filter(data_C_australis, Metric == "Fv.Fm")

# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_C_australis$Tissue.edit <- factor(data_C_australis$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))

Fv.Fm_C_australis_boxplot <- ggplot(data_C_australis_Fv.Fm, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5),        
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 1)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

Fv.Fm_C_australis_boxplot 

# use base R boxplot to get the coordinates of the boxes
box.rslt_Fv.Fm_C_australis <- with(data_C_australis_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_australis)
boxplot_positions_Fv.Fm_C_australis <- as.data.frame(box.rslt_Fv.Fm_C_australis$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_australis <- levels(data_C_australis_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_australis) <- tissues_Fv.Fm_C_australis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_australis <- boxplot_positions_Fv.Fm_C_australis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_australis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_australis__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_australis)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_australis, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_australis %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_australis
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_australis, top_positions_Fv.Fm_C_australis, by = "Tissue.edit") -> cbd_Fv.Fm_C_australis

# calculate how much to nudge
data_C_australis_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_australis
cbd_Fv.Fm_C_australis$nudged <- max_Fv.Fm_C_australis$max * 1.05


# add CLDs to plot
Fv.Fm_C_australis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_australis,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_australis_boxplot


Fv.Fm_C_australis_boxplot


#### φPSII C_australis ####
# calculate minimum φPSII for apporpriate bottom limit of φPSII y axes
min_φPSII <-as.numeric( data_no_outliers %>% dplyr::filter(., Metric == "φPSII") %>% dplyr::summarize(., min(Value)))
min_φPSII

data_C_australis_φPSII <- dplyr::filter(data_C_australis, Metric == "φPSII")


φPSII_C_australis_boxplot <- ggplot(data_C_australis_φPSII, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(min_φPSII, 1))+
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
box.rslt_φPSII_C_australis <- with(data_C_australis_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_australis)
boxplot_positions_φPSII_C_australis <- as.data.frame(box.rslt_φPSII_C_australis$stats)

# what are these column tissue codes?
tissues_φPSII_C_australis <- levels(data_C_australis_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_australis) <- tissues_φPSII_C_australis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_australis <- boxplot_positions_φPSII_C_australis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_australis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_australis__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_australis)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_australis, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_australis %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_australis
# now join these positions to cbd
left_join(cbd_φPSII_C_australis, top_positions_φPSII_C_australis, by = "Tissue.edit") -> cbd_φPSII_C_australis

# calculate how much to nudge
data_C_australis_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_australis
cbd_φPSII_C_australis$nudged <- max_φPSII_C_australis$max * 1.05


# add CLDs to plot
φPSII_C_australis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_australis,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_australis_boxplot


φPSII_C_australis_boxplot


#### NPQ C_australis plot with x axis ####
data_C_australis_NPQ <- dplyr::filter(data_C_australis, Metric == "NPQ")


NPQ_C_australis_boxplot <- ggplot(data_C_australis_NPQ, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 4)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

NPQ_C_australis_boxplot


#### Fv.Fm loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# filter for Fv.Fm
data_Fv.Fm_plots_grammicaspe <- dplyr::filter(data_no_outliers_Grammica_plot, Metric == "Fv.Fm")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Fv.Fm_plots_grammicaspe$Tissue.edit <- factor(data_Fv.Fm_plots_grammicaspe$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Fv.Fm = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "left", limits = c(0, 1))
  plot_list_Grammica_Fv.Fm[[sub]] = p
  
}


#### Fv.Fm C_polygonorum ####
plot_list_Grammica_Fv.Fm[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_polygonorum_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Fv.Fm_C_polygonorum <- with(data_C_polygonorum_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_polygonorum)
boxplot_positions_Fv.Fm_C_polygonorum <- as.data.frame(box.rslt_Fv.Fm_C_polygonorum$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_polygonorum <- levels(data_C_polygonorum_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_polygonorum) <- tissues_Fv.Fm_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_polygonorum <- boxplot_positions_Fv.Fm_C_polygonorum[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_polygonorum, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_polygonorum %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_polygonorum
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_polygonorum, top_positions_Fv.Fm_C_polygonorum, by = "Tissue.edit") -> cbd_Fv.Fm_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_polygonorum
cbd_Fv.Fm_C_polygonorum$nudged <- max_Fv.Fm_C_polygonorum$max * 1.05


# add CLDs to plot
Fv.Fm_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_polygonorum_boxplot


Fv.Fm_C_polygonorum_boxplot



#### Fv.Fm C_sandwichiana ####
plot_list_Grammica_Fv.Fm[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_sandwichiana_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_Fv.Fm_C_sandwichiana <- with(data_C_sandwichiana_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_sandwichiana)
boxplot_positions_Fv.Fm_C_sandwichiana <- as.data.frame(box.rslt_Fv.Fm_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_sandwichiana <- levels(data_C_sandwichiana_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_sandwichiana) <- tissues_Fv.Fm_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_sandwichiana <- boxplot_positions_Fv.Fm_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_sandwichiana, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_sandwichiana %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_sandwichiana
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_sandwichiana, top_positions_Fv.Fm_C_sandwichiana, by = "Tissue.edit") -> cbd_Fv.Fm_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_sandwichiana
cbd_Fv.Fm_C_sandwichiana$nudged <- max_Fv.Fm_C_sandwichiana$max * 1.05


# add CLDs to plot
Fv.Fm_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_sandwichiana_boxplot


Fv.Fm_C_sandwichiana_boxplot



#### Fv.Fm C_californica ####
plot_list_Grammica_Fv.Fm[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_californica")

box.rslt_Fv.Fm_C_californica <- with(data_C_californica_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_californica)
boxplot_positions_Fv.Fm_C_californica <- as.data.frame(box.rslt_Fv.Fm_C_californica$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_californica <- levels(data_C_californica_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_californica) <- tissues_Fv.Fm_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_californica <- boxplot_positions_Fv.Fm_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_californica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_californica %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_californica
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_californica, top_positions_Fv.Fm_C_californica, by = "Tissue.edit") -> cbd_Fv.Fm_C_californica

# calculate how much to nudge
data_C_californica_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_californica
cbd_Fv.Fm_C_californica$nudged <- max_Fv.Fm_C_californica$max * 1.05


# add CLDs to plot
Fv.Fm_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_californica_boxplot


Fv.Fm_C_californica_boxplot


#### Fv.Fm C_compacta ####
plot_list_Grammica_Fv.Fm[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_compacta")

box.rslt_Fv.Fm_C_compacta <- with(data_C_compacta_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_compacta)
boxplot_positions_Fv.Fm_C_compacta <- as.data.frame(box.rslt_Fv.Fm_C_compacta$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_compacta <- levels(data_C_compacta_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_compacta) <- tissues_Fv.Fm_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_compacta <- boxplot_positions_Fv.Fm_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_compacta, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_compacta %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_compacta
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_compacta, top_positions_Fv.Fm_C_compacta, by = "Tissue.edit") -> cbd_Fv.Fm_C_compacta

# calculate how much to nudge
data_C_compacta_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_compacta
cbd_Fv.Fm_C_compacta$nudged <- max_Fv.Fm_C_compacta$max * 1.05


# add CLDs to plot
Fv.Fm_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_compacta_boxplot


Fv.Fm_C_compacta_boxplot



#### Fv.Fm C_cephalanthii ####
plot_list_Grammica_Fv.Fm[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Fv.Fm_C_cephalanthii <- with(data_C_cephalanthii_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_cephalanthii)
boxplot_positions_Fv.Fm_C_cephalanthii <- as.data.frame(box.rslt_Fv.Fm_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_cephalanthii <- levels(data_C_cephalanthii_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_cephalanthii) <- tissues_Fv.Fm_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_cephalanthii <- boxplot_positions_Fv.Fm_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_cephalanthii, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_cephalanthii %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_cephalanthii, top_positions_Fv.Fm_C_cephalanthii, by = "Tissue.edit") -> cbd_Fv.Fm_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_cephalanthii
cbd_Fv.Fm_C_cephalanthii$nudged <- max_Fv.Fm_C_cephalanthii$max * 1.05


# add CLDs to plot
Fv.Fm_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_cephalanthii_boxplot


Fv.Fm_C_cephalanthii_boxplot


#### Fv.Fm C_denticulata ####
plot_list_Grammica_Fv.Fm[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_denticulata_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_denticulata")

box.rslt_Fv.Fm_C_denticulata <- with(data_C_denticulata_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_denticulata)
boxplot_positions_Fv.Fm_C_denticulata <- as.data.frame(box.rslt_Fv.Fm_C_denticulata$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_denticulata <- levels(data_C_denticulata_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_denticulata) <- tissues_Fv.Fm_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_denticulata <- boxplot_positions_Fv.Fm_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_denticulata, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_denticulata %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_denticulata
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_denticulata, top_positions_Fv.Fm_C_denticulata, by = "Tissue.edit") -> cbd_Fv.Fm_C_denticulata

# calculate how much to nudge
data_C_denticulata_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_denticulata
cbd_Fv.Fm_C_denticulata$nudged <- max_Fv.Fm_C_denticulata$max * 1.05


# add CLDs to plot
Fv.Fm_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_denticulata_boxplot


Fv.Fm_C_denticulata_boxplot



#### Fv.Fm C_tasmanica ####
plot_list_Grammica_Fv.Fm[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Fv.Fm_C_tasmanica <- with(data_C_tasmanica_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_tasmanica)
boxplot_positions_Fv.Fm_C_tasmanica <- as.data.frame(box.rslt_Fv.Fm_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_tasmanica <- levels(data_C_tasmanica_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_tasmanica) <- tissues_Fv.Fm_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_tasmanica <- boxplot_positions_Fv.Fm_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_tasmanica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_tasmanica %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_tasmanica
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_tasmanica, top_positions_Fv.Fm_C_tasmanica, by = "Tissue.edit") -> cbd_Fv.Fm_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_tasmanica
cbd_Fv.Fm_C_tasmanica$nudged <- max_Fv.Fm_C_tasmanica$max * 1.05


# add CLDs to plot
Fv.Fm_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_tasmanica_boxplot


Fv.Fm_C_tasmanica_boxplot


#### Fv.Fm C_costaricensis ####
plot_list_Grammica_Fv.Fm[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Fv.Fm_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Fv.Fm_C_costaricensis <- with(data_C_costaricensis_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_costaricensis)
boxplot_positions_Fv.Fm_C_costaricensis <- as.data.frame(box.rslt_Fv.Fm_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_costaricensis <- levels(data_C_costaricensis_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_costaricensis) <- tissues_Fv.Fm_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_costaricensis <- boxplot_positions_Fv.Fm_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_costaricensis, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_costaricensis %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_costaricensis
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_costaricensis, top_positions_Fv.Fm_C_costaricensis, by = "Tissue.edit") -> cbd_Fv.Fm_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_costaricensis
cbd_Fv.Fm_C_costaricensis$nudged <- max_Fv.Fm_C_costaricensis$max * 1.05


# add CLDs to plot
Fv.Fm_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_costaricensis_boxplot


Fv.Fm_C_costaricensis_boxplot


#### Fv.Fm C_indecora species alone ####

data_Grammica_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_indecora")

Fv.Fm_C_indecora_boxplot <- ggplot(data_Grammica_Fv.Fm, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 1)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "Fv.Fm" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Fv.Fm <- dplyr::filter(data_Fv.Fm_plots_grammicaspe, Species == "C_indecora")

box.rslt_Fv.Fm_C_indecora <- with(data_C_indecora_Fv.Fm, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_Fv.Fm_C_indecora)
boxplot_positions_Fv.Fm_C_indecora <- as.data.frame(box.rslt_Fv.Fm_C_indecora$stats)

# what are these column tissue codes?
tissues_Fv.Fm_C_indecora <- levels(data_C_indecora_Fv.Fm$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_Fv.Fm_C_indecora) <- tissues_Fv.Fm_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Fv.Fm_C_indecora <- boxplot_positions_Fv.Fm_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Fv.Fm_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Fv.Fm"]])[["Letters"]])
colnames(cbd_Fv.Fm_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_Fv.Fm_C_indecora, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_Fv.Fm_C_indecora %>% gather(., Tissue.edit, y.position) -> top_positions_Fv.Fm_C_indecora
# now join these positions to cbd
left_join(cbd_Fv.Fm_C_indecora, top_positions_Fv.Fm_C_indecora, by = "Tissue.edit") -> cbd_Fv.Fm_C_indecora

# calculate how much to nudge
data_C_indecora_Fv.Fm %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_Fv.Fm_C_indecora
cbd_Fv.Fm_C_indecora$nudged <- (max_Fv.Fm_C_indecora$max + 0.00) * 1.05

# add CLDs to plot
Fv.Fm_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Fv.Fm_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> Fv.Fm_C_indecora_boxplot


Fv.Fm_C_indecora_boxplot



#### φPSII loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# filter for φPSII
data_φPSII_plots_grammicaspe <- dplyr::filter(data_no_outliers_Grammica_plot, Metric == "φPSII")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_φPSII_plots_grammicaspe$Tissue.edit <- factor(data_φPSII_plots_grammicaspe$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_φPSII = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on C_sandwichiana plots
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(min_φPSII, 1)) 
  plot_list_Grammica_φPSII[[sub]] = p
  
}


#### φPSII C_polygonorum ####
plot_list_Grammica_φPSII[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_polygonorum_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_polygonorum")

box.rslt_φPSII_C_polygonorum <- with(data_C_polygonorum_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_polygonorum)
boxplot_positions_φPSII_C_polygonorum <- as.data.frame(box.rslt_φPSII_C_polygonorum$stats)

# what are these column tissue codes?
tissues_φPSII_C_polygonorum <- levels(data_C_polygonorum_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_polygonorum) <- tissues_φPSII_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_polygonorum <- boxplot_positions_φPSII_C_polygonorum[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_polygonorum, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_polygonorum %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_polygonorum
# now join these positions to cbd
left_join(cbd_φPSII_C_polygonorum, top_positions_φPSII_C_polygonorum, by = "Tissue.edit") -> cbd_φPSII_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_polygonorum
cbd_φPSII_C_polygonorum$nudged <- max_φPSII_C_polygonorum$max * 1.05


# add CLDs to plot
φPSII_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_polygonorum_boxplot


φPSII_C_polygonorum_boxplot



#### φPSII C_sandwichiana ####
plot_list_Grammica_φPSII[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_sandwichiana_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_φPSII_C_sandwichiana <- with(data_C_sandwichiana_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_sandwichiana)
boxplot_positions_φPSII_C_sandwichiana <- as.data.frame(box.rslt_φPSII_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_φPSII_C_sandwichiana <- levels(data_C_sandwichiana_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_sandwichiana) <- tissues_φPSII_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_sandwichiana <- boxplot_positions_φPSII_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_sandwichiana, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_sandwichiana %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_sandwichiana
# now join these positions to cbd
left_join(cbd_φPSII_C_sandwichiana, top_positions_φPSII_C_sandwichiana, by = "Tissue.edit") -> cbd_φPSII_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_sandwichiana
cbd_φPSII_C_sandwichiana$nudged <- max_φPSII_C_sandwichiana$max * 1.05


# add CLDs to plot
φPSII_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_sandwichiana_boxplot


φPSII_C_sandwichiana_boxplot


#### φPSII C_californica ####
plot_list_Grammica_φPSII[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_californica")

box.rslt_φPSII_C_californica <- with(data_C_californica_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_californica)
boxplot_positions_φPSII_C_californica <- as.data.frame(box.rslt_φPSII_C_californica$stats)

# what are these column tissue codes?
tissues_φPSII_C_californica <- levels(data_C_californica_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_californica) <- tissues_φPSII_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_californica <- boxplot_positions_φPSII_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_californica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_californica %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_californica
# now join these positions to cbd
left_join(cbd_φPSII_C_californica, top_positions_φPSII_C_californica, by = "Tissue.edit") -> cbd_φPSII_C_californica

# calculate how much to nudge
data_C_californica_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_californica
cbd_φPSII_C_californica$nudged <- max_φPSII_C_californica$max * 1.05


# add CLDs to plot
φPSII_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_californica_boxplot


φPSII_C_californica_boxplot


#### φPSII C_compacta ####
plot_list_Grammica_φPSII[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_compacta")

box.rslt_φPSII_C_compacta <- with(data_C_compacta_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_compacta)
boxplot_positions_φPSII_C_compacta <- as.data.frame(box.rslt_φPSII_C_compacta$stats)

# what are these column tissue codes?
tissues_φPSII_C_compacta <- levels(data_C_compacta_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_compacta) <- tissues_φPSII_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_compacta <- boxplot_positions_φPSII_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_compacta, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_compacta %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_compacta
# now join these positions to cbd
left_join(cbd_φPSII_C_compacta, top_positions_φPSII_C_compacta, by = "Tissue.edit") -> cbd_φPSII_C_compacta

# calculate how much to nudge
data_C_compacta_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_compacta
cbd_φPSII_C_compacta$nudged <- max_φPSII_C_compacta$max * 1.05


# add CLDs to plot
φPSII_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_compacta_boxplot


φPSII_C_compacta_boxplot


#### φPSII C_cephalanthii ####
plot_list_Grammica_φPSII[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_φPSII_C_cephalanthii <- with(data_C_cephalanthii_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_cephalanthii)
boxplot_positions_φPSII_C_cephalanthii <- as.data.frame(box.rslt_φPSII_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_φPSII_C_cephalanthii <- levels(data_C_cephalanthii_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_cephalanthii) <- tissues_φPSII_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_cephalanthii <- boxplot_positions_φPSII_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_cephalanthii, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_cephalanthii %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_cephalanthii
# now join these positions to cbd
left_join(cbd_φPSII_C_cephalanthii, top_positions_φPSII_C_cephalanthii, by = "Tissue.edit") -> cbd_φPSII_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_cephalanthii
cbd_φPSII_C_cephalanthii$nudged <- max_φPSII_C_cephalanthii$max * 1.05


# add CLDs to plot
φPSII_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_cephalanthii_boxplot


φPSII_C_cephalanthii_boxplot


#### φPSII C_denticulata ####
plot_list_Grammica_φPSII[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_denticulata_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_denticulata")

box.rslt_φPSII_C_denticulata <- with(data_C_denticulata_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_denticulata)
boxplot_positions_φPSII_C_denticulata <- as.data.frame(box.rslt_φPSII_C_denticulata$stats)

# what are these column tissue codes?
tissues_φPSII_C_denticulata <- levels(data_C_denticulata_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_denticulata) <- tissues_φPSII_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_denticulata <- boxplot_positions_φPSII_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_denticulata, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_denticulata %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_denticulata
# now join these positions to cbd
left_join(cbd_φPSII_C_denticulata, top_positions_φPSII_C_denticulata, by = "Tissue.edit") -> cbd_φPSII_C_denticulata

# calculate how much to nudge
data_C_denticulata_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_denticulata
cbd_φPSII_C_denticulata$nudged <- max_φPSII_C_denticulata$max * 1.05


# add CLDs to plot
φPSII_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_denticulata_boxplot


φPSII_C_denticulata_boxplot



#### φPSII C_tasmanica ####
plot_list_Grammica_φPSII[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_φPSII_C_tasmanica <- with(data_C_tasmanica_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_tasmanica)
boxplot_positions_φPSII_C_tasmanica <- as.data.frame(box.rslt_φPSII_C_tasmanica$stats)

# what are these column tissue codes?
tissues_φPSII_C_tasmanica <- levels(data_C_tasmanica_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_tasmanica) <- tissues_φPSII_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_tasmanica <- boxplot_positions_φPSII_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_tasmanica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_tasmanica %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_tasmanica
# now join these positions to cbd
left_join(cbd_φPSII_C_tasmanica, top_positions_φPSII_C_tasmanica, by = "Tissue.edit") -> cbd_φPSII_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_tasmanica
cbd_φPSII_C_tasmanica$nudged <- max_φPSII_C_tasmanica$max * 1.05


# add CLDs to plot
φPSII_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_tasmanica_boxplot


φPSII_C_tasmanica_boxplot


#### φPSII C_costaricensis ####
plot_list_Grammica_φPSII[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> φPSII_C_costaricensis_boxplot

# KW not sig no post hoc
φPSII_C_costaricensis_boxplot


#### φPSII C_indecora species alone ####

data_Grammica_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_indecora")

φPSII_C_indecora_boxplot <- ggplot(data_Grammica_φPSII, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(min_φPSII, 1)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "φPSII" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_φPSII <- dplyr::filter(data_φPSII_plots_grammicaspe, Species == "C_indecora")

box.rslt_φPSII_C_indecora <- with(data_C_indecora_φPSII, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_φPSII_C_indecora)
boxplot_positions_φPSII_C_indecora <- as.data.frame(box.rslt_φPSII_C_indecora$stats)

# what are these column tissue codes?
tissues_φPSII_C_indecora <- levels(data_C_indecora_φPSII$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_φPSII_C_indecora) <- tissues_φPSII_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_φPSII_C_indecora <- boxplot_positions_φPSII_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_φPSII_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__φPSII"]])[["Letters"]])
colnames(cbd_φPSII_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_φPSII_C_indecora, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_φPSII_C_indecora %>% gather(., Tissue.edit, y.position) -> top_positions_φPSII_C_indecora
# now join these positions to cbd
left_join(cbd_φPSII_C_indecora, top_positions_φPSII_C_indecora, by = "Tissue.edit") -> cbd_φPSII_C_indecora

# calculate how much to nudge
data_C_indecora_φPSII %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_φPSII_C_indecora
cbd_φPSII_C_indecora$nudged <- (max_φPSII_C_indecora$max * 1.05)

# add CLDs to plot
φPSII_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_φPSII_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> φPSII_C_indecora_boxplot


φPSII_C_indecora_boxplot


#### NPQ loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# filter for NPQ
data_NPQ_plots_grammicaspe <- dplyr::filter(data_no_outliers_Grammica_plot, Metric == "NPQ")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_NPQ_plots_grammicaspe$Tissue.edit <- factor(data_NPQ_plots_grammicaspe$Tissue.edit, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_NPQ = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on C_sandwichiana plots
          axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 4))
  plot_list_Grammica_NPQ[[sub]] = p
  
}

#### NPQ C_polygonorum ####
plot_list_Grammica_NPQ[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_polygonorum_boxplot

# KW not sig so no post hoc
NPQ_C_polygonorum_boxplot


#### NPQ C_sandwichiana ####
plot_list_Grammica_NPQ[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_sandwichiana_boxplot


# KW not sig so no post hoc
NPQ_C_sandwichiana_boxplot


#### NPQ C_californica ####
plot_list_Grammica_NPQ[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_californica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_californica_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_californica")

box.rslt_NPQ_C_californica <- with(data_C_californica_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_californica)
boxplot_positions_NPQ_C_californica <- as.data.frame(box.rslt_NPQ_C_californica$stats)

# what are these column tissue codes?
tissues_NPQ_C_californica <- levels(data_C_californica_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_californica) <- tissues_NPQ_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_californica <- boxplot_positions_NPQ_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_californica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_californica %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_californica
# now join these positions to cbd
left_join(cbd_NPQ_C_californica, top_positions_NPQ_C_californica, by = "Tissue.edit") -> cbd_NPQ_C_californica

# calculate how much to nudge
data_C_californica_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_californica
cbd_NPQ_C_californica$nudged <- max_NPQ_C_californica$max * 1.05


# add CLDs to plot
NPQ_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_californica_boxplot


NPQ_C_californica_boxplot



#### NPQ C_compacta ####
plot_list_Grammica_NPQ[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_compacta_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_compacta_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_compacta")

box.rslt_NPQ_C_compacta <- with(data_C_compacta_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_compacta)
boxplot_positions_NPQ_C_compacta <- as.data.frame(box.rslt_NPQ_C_compacta$stats)

# what are these column tissue codes?
tissues_NPQ_C_compacta <- levels(data_C_compacta_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_compacta) <- tissues_NPQ_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_compacta <- boxplot_positions_NPQ_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_compacta, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_compacta %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_compacta
# now join these positions to cbd
left_join(cbd_NPQ_C_compacta, top_positions_NPQ_C_compacta, by = "Tissue.edit") -> cbd_NPQ_C_compacta

# calculate how much to nudge
data_C_compacta_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_compacta
cbd_NPQ_C_compacta$nudged <- max_NPQ_C_compacta$max * 1.05


# add CLDs to plot
NPQ_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_compacta_boxplot


NPQ_C_compacta_boxplot



#### NPQ C_cephalanthii ####
plot_list_Grammica_NPQ[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_NPQ_C_cephalanthii <- with(data_C_cephalanthii_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_cephalanthii)
boxplot_positions_NPQ_C_cephalanthii <- as.data.frame(box.rslt_NPQ_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_NPQ_C_cephalanthii <- levels(data_C_cephalanthii_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_cephalanthii) <- tissues_NPQ_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_cephalanthii <- boxplot_positions_NPQ_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_cephalanthii, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_cephalanthii %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_cephalanthii
# now join these positions to cbd
left_join(cbd_NPQ_C_cephalanthii, top_positions_NPQ_C_cephalanthii, by = "Tissue.edit") -> cbd_NPQ_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_cephalanthii
cbd_NPQ_C_cephalanthii$nudged <- max_NPQ_C_cephalanthii$max * 1.05


# add CLDs to plot
NPQ_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_cephalanthii_boxplot


NPQ_C_cephalanthii_boxplot


#### NPQ C_denticulata ####
plot_list_Grammica_NPQ[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_denticulata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_denticulata")

box.rslt_NPQ_C_denticulata <- with(data_C_denticulata_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_denticulata)
boxplot_positions_NPQ_C_denticulata <- as.data.frame(box.rslt_NPQ_C_denticulata$stats)

# what are these column tissue codes?
tissues_NPQ_C_denticulata <- levels(data_C_denticulata_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_denticulata) <- tissues_NPQ_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_denticulata <- boxplot_positions_NPQ_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_denticulata, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_denticulata %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_denticulata
# now join these positions to cbd
left_join(cbd_NPQ_C_denticulata, top_positions_NPQ_C_denticulata, by = "Tissue.edit") -> cbd_NPQ_C_denticulata

# calculate how much to nudge
data_C_denticulata_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_denticulata
cbd_NPQ_C_denticulata$nudged <- max_NPQ_C_denticulata$max * 1.05


# add CLDs to plot
NPQ_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_denticulata_boxplot


NPQ_C_denticulata_boxplot




#### NPQ C_tasmanica ####
plot_list_Grammica_NPQ[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_NPQ_C_tasmanica <- with(data_C_tasmanica_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_tasmanica)
boxplot_positions_NPQ_C_tasmanica <- as.data.frame(box.rslt_NPQ_C_tasmanica$stats)

# what are these column tissue codes?
tissues_NPQ_C_tasmanica <- levels(data_C_tasmanica_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_tasmanica) <- tissues_NPQ_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_tasmanica <- boxplot_positions_NPQ_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_tasmanica, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_tasmanica %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_tasmanica
# now join these positions to cbd
left_join(cbd_NPQ_C_tasmanica, top_positions_NPQ_C_tasmanica, by = "Tissue.edit") -> cbd_NPQ_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_tasmanica
cbd_NPQ_C_tasmanica$nudged <- max_NPQ_C_tasmanica$max * 1.05


# add CLDs to plot
NPQ_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_tasmanica_boxplot


NPQ_C_tasmanica_boxplot


#### NPQ C_costaricensis ####
plot_list_Grammica_NPQ[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NPQ_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_NPQ_C_costaricensis <- with(data_C_costaricensis_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_costaricensis)
boxplot_positions_NPQ_C_costaricensis <- as.data.frame(box.rslt_NPQ_C_costaricensis$stats)

# what are these column tissue codes?
tissues_NPQ_C_costaricensis <- levels(data_C_costaricensis_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_costaricensis) <- tissues_NPQ_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_costaricensis <- boxplot_positions_NPQ_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_costaricensis, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_costaricensis %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_costaricensis
# now join these positions to cbd
left_join(cbd_NPQ_C_costaricensis, top_positions_NPQ_C_costaricensis, by = "Tissue.edit") -> cbd_NPQ_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_costaricensis
cbd_NPQ_C_costaricensis$nudged <- max_NPQ_C_costaricensis$max * 1.05


# add CLDs to plot
NPQ_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_costaricensis_boxplot


NPQ_C_costaricensis_boxplot


#### NPQ C_indecora species alone ####

data_Grammica_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_indecora")

NPQ_C_indecora_boxplot <- ggplot(data_Grammica_NPQ, aes(x=Tissue.edit, y=Value, color=Tissue.edit)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x = element_text(angle = -30, hjust = 0, size = 4),
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 4)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Metric == "NPQ" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
data_C_indecora_NPQ <- dplyr::filter(data_NPQ_plots_grammicaspe, Species == "C_indecora")

box.rslt_NPQ_C_indecora <- with(data_C_indecora_NPQ, graphics::boxplot(Value ~ Tissue.edit, plot = FALSE))
str(box.rslt_NPQ_C_indecora)
boxplot_positions_NPQ_C_indecora <- as.data.frame(box.rslt_NPQ_C_indecora$stats)

# what are these column tissue codes?
tissues_NPQ_C_indecora <- levels(data_C_indecora_NPQ$Tissue.edit)
# add appropriate tissues to position df
colnames(boxplot_positions_NPQ_C_indecora) <- tissues_NPQ_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NPQ_C_indecora <- boxplot_positions_NPQ_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NPQ_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__NPQ"]])[["Letters"]])
colnames(cbd_NPQ_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.edit
setDT(cbd_NPQ_C_indecora, keep.rownames = "Tissue.edit")


# add a column y.position taken from top_positions based on mtaching up Tissue.edit
# first reshape top_positions so that colnames are a column called Tissue.edit
top_positions_NPQ_C_indecora %>% gather(., Tissue.edit, y.position) -> top_positions_NPQ_C_indecora
# now join these positions to cbd
left_join(cbd_NPQ_C_indecora, top_positions_NPQ_C_indecora, by = "Tissue.edit") -> cbd_NPQ_C_indecora

# calculate how much to nudge
data_C_indecora_NPQ %>% group_by(Tissue.edit) %>% dplyr::summarize(., max = max(Value)) -> max_NPQ_C_indecora
cbd_NPQ_C_indecora$nudged <- (max_NPQ_C_indecora$max * 1.05)

# add CLDs to plot
NPQ_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NPQ_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.edit, y = nudged, label = Letter, vjust = 0)) -> NPQ_C_indecora_boxplot


NPQ_C_indecora_boxplot



#### add "absent" to those with mean Value = 0 in loop ####

Species_list_Grammica <- unique(data_no_outliers_Grammica_plot$Species)
Metric_list_Grammica <- unique(data_no_outliers_Grammica_plot$Metric)
label <- "absent"
absent <- data.frame(label)

Grammica_carotenoid_absent_list <- list()

for (species in Species_list_Grammica) {
  for (metric in Metric_list_Grammica) {
    data_loop <- data_no_outliers_Grammica_plot %>% dplyr::filter(Species == species & Metric == metric)
    if (mean(data_loop$Value) == 0) {
      print("absent") -> Grammica_carotenoid_absent_list[[metric]][[species]]
    }}}

# none

# add the following geom_text to the above plots
# + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) 





#### COMBINED (faceted) fluorescence plot: Grammica ONLY ####
wrap_elements(gridtext::richtext_grob('*F*<sub>v</sub>/*F*<sub>m</sub>', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Fv.Fm_C_australis_boxplot + ggtitle('C. australis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_polygonorum_boxplot + ggtitle('C. polygonorum') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_sandwichiana_boxplot + ggtitle('C. sandwichiana') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_californica_boxplot + ggtitle('C. californica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_compacta_boxplot + ggtitle('C. compacta') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_cephalanthii_boxplot + ggtitle('C. cephalanthii') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_denticulata_boxplot + ggtitle('C. denticulata') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Fv.Fm_C_tasmanica_boxplot + ggtitle('C. tasmanica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Fv.Fm_C_costaricensis_boxplot + ggtitle('C. costaricensis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Fv.Fm_C_indecora_boxplot + ggtitle('C. indecora') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  wrap_elements(gridtext::richtext_grob('φPSII', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  φPSII_C_australis_boxplot +
  φPSII_C_polygonorum_boxplot + 
  φPSII_C_sandwichiana_boxplot +
  φPSII_C_californica_boxplot + 
  φPSII_C_compacta_boxplot +
  φPSII_C_cephalanthii_boxplot +
  φPSII_C_denticulata_boxplot +
  φPSII_C_tasmanica_boxplot +
  φPSII_C_costaricensis_boxplot +
  φPSII_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('NPQ', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  NPQ_C_australis_boxplot +
  NPQ_C_polygonorum_boxplot + 
  NPQ_C_sandwichiana_boxplot +
  NPQ_C_californica_boxplot + 
  NPQ_C_compacta_boxplot +
  NPQ_C_cephalanthii_boxplot +
  NPQ_C_denticulata_boxplot +
  NPQ_C_tasmanica_boxplot +
  NPQ_C_costaricensis_boxplot +
  NPQ_C_indecora_boxplot +
  plot_layout(nrow = 3, byrow = T) -> fluorescence_boxplot_Grammica_new

fluorescence_boxplot_Grammica_new 

pdf("../output/stat_results/boxplots/fluorescence_boxplot_Grammica_new.pdf", width=9,height=5.5) 
fluorescence_boxplot_Grammica_new
dev.off()


