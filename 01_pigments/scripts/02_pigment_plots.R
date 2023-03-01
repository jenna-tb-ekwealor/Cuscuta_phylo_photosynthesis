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
data_long_calcs <- read.csv(file = "../output/stat_results/data_long_calcs_for_plots.csv", stringsAsFactors = T)
dat_text_plot_kruskal <- read.csv(file = "../output/stat_results/dat_text_plot_kruskal.csv")
data_long_calcs_Grammica_plot <- read.csv(file = "../output/stat_results/data_long_calcs_for_Grammica_plots.csv", stringsAsFactors = T)
dat_text_plot_kruskal_Grammica <- read.csv(file = "../output/stat_results/dat_text_plot_kruskal_Grammica.csv")
summary_accession <- read.csv(file = "../output/stat_results/pigments_species_summary.csv", stringsAsFactors = T)
wilcox_list <- readRDS("../output/stat_results/wilcox_list.RData")
wilcox_list_Grammica <- readRDS("../output/stat_results/wilcox_list_Grammica.RData")


#### CHLOROPHYLL PLOTS ####

#### Chl.a Ipomoea_nil ####
data_ipomoea <-  dplyr::filter(data_long_calcs, Subgenus == "Ipomoea_nil")
data_ipomoea_Chl.a <- dplyr::filter(data_ipomoea, Pigment == "Chl.a")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_ipomoea$Tissue.code <- factor(data_ipomoea$Tissue.code, levels = c("l", "y", "o", "f", "s"))

Chl.a_ipomoea_boxplot <- ggplot(data_ipomoea_Chl.a, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 8), breaks = seq(0, 8, by = 2)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

Chl.a_ipomoea_boxplot 

# use base R boxplot to get the coordinates of the boxes
box.rslt_Chl.a_ipomoea <- with(data_ipomoea_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_ipomoea)
boxplot_positions_Chl.a_ipomoea <- as.data.frame(box.rslt_Chl.a_ipomoea$stats)

# what are these column tissue codes?
tissues_Chl.a_ipomoea <- levels(data_ipomoea_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_ipomoea) <- tissues_Chl.a_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_ipomoea <- boxplot_positions_Chl.a_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
data.table::setDT(cbd_Chl.a_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_ipomoea %>% tidyr::gather(., Tissue.code, y.position) -> top_positions_Chl.a_ipomoea
# now join these positions to cbd
left_join(cbd_Chl.a_ipomoea, top_positions_Chl.a_ipomoea, by = "Tissue.code") -> cbd_Chl.a_ipomoea

# calculate how much to nudge
data_ipomoea_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_ipomoea
cbd_Chl.a_ipomoea$nudged <- max_Chl.a_ipomoea$max * 1.05


# add CLDs to plot
Chl.a_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_ipomoea_boxplot

Chl.a_ipomoea_boxplot



#### Chl.b Ipomoea_nil ####
data_ipomoea_Chl.b <- dplyr::filter(data_ipomoea, Pigment == "Chl.b")


Chl.b_ipomoea_boxplot <- ggplot(data_ipomoea_Chl.b, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.b" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
box.rslt_Chl.b_ipomoea <- with(data_ipomoea_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_ipomoea)
boxplot_positions_Chl.b_ipomoea <- as.data.frame(box.rslt_Chl.b_ipomoea$stats)

# what are these column tissue codes?
tissues_Chl.b_ipomoea <- levels(data_ipomoea_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_ipomoea) <- tissues_Chl.b_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_ipomoea <- boxplot_positions_Chl.b_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
data.table::setDT(cbd_Chl.b_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_ipomoea
# now join these positions to cbd
left_join(cbd_Chl.b_ipomoea, top_positions_Chl.b_ipomoea, by = "Tissue.code") -> cbd_Chl.b_ipomoea

# calculate how much to nudge
data_ipomoea_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_ipomoea
cbd_Chl.b_ipomoea$nudged <- max_Chl.b_ipomoea$max * 1.05



# add CLDs to plot
Chl.b_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.b_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_ipomoea_boxplot


Chl.b_ipomoea_boxplot


#### Tot.Chl Ipomoea_nil ####
data_ipomoea_Tot.Chl <- dplyr::filter(data_ipomoea, Pigment == "Tot.Chl")


Tot.Chl_ipomoea_boxplot <- ggplot(data_ipomoea_Tot.Chl, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 10), breaks = seq(0, 10, by = 2)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Chl" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_Tot.Chl_ipomoea <- with(data_ipomoea_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_ipomoea)
boxplot_positions_Tot.Chl_ipomoea <- as.data.frame(box.rslt_Tot.Chl_ipomoea$stats)

# what are these column tissue codes?
tissues_Tot.Chl_ipomoea <- levels(data_ipomoea_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_ipomoea) <- tissues_Tot.Chl_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_ipomoea <- boxplot_positions_Tot.Chl_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_ipomoea
# now join these positions to cbd
left_join(cbd_Tot.Chl_ipomoea, top_positions_Tot.Chl_ipomoea, by = "Tissue.code") -> cbd_Tot.Chl_ipomoea

# calculate how much to nudge
data_ipomoea_Tot.Chl %>% group_by(Tissue.code) %>%  dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_ipomoea
cbd_Tot.Chl_ipomoea$nudged <- max_Tot.Chl_ipomoea$max* 1.05


# add CLDs to plot
Tot.Chl_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Chl_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_ipomoea_boxplot


Tot.Chl_ipomoea_boxplot


#### Chl.a.b Ipomoea_nil plot with x axis ####
data_ipomoea_Chl.a.b <- dplyr::filter(data_ipomoea, Pigment == "Chl.a.b")


Chl.a.b_ipomoea_boxplot <- ggplot(data_ipomoea_Chl.a.b, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Mass ratio") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 10), breaks = seq(0, 10, by = 2))+ 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a.b" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_Chl.a.b_ipomoea <- with(data_ipomoea_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_ipomoea)
boxplot_positions_Chl.a.b_ipomoea <- as.data.frame(box.rslt_Chl.a.b_ipomoea$stats)

# what are these column tissue codes?
tissues_Chl.a.b_ipomoea <- levels(data_ipomoea_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_ipomoea) <- tissues_Chl.a.b_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_ipomoea <- boxplot_positions_Chl.a.b_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_ipomoea
# now join these positions to cbd
left_join(cbd_Chl.a.b_ipomoea, top_positions_Chl.a.b_ipomoea, by = "Tissue.code") -> cbd_Chl.a.b_ipomoea

# calculate how much to nudge
data_ipomoea_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_ipomoea
cbd_Chl.a.b_ipomoea$nudged <- max_Chl.a.b_ipomoea$max * 1.05


# add CLDs to plot
Chl.a.b_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a.b_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_ipomoea_boxplot


Chl.a.b_ipomoea_boxplot


#### Chl.a loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Chl.a
data_Chl.a_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Chl.a")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.a_plots_cuscutasub$Tissue.code <- factor(data_Chl.a_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Chl.a = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
    scale_y_continuous(position = "left", limits = c(0, 8), breaks = seq(0, 8, by = 2))
  plot_list_Cuscuta_Chl.a[[sub]] = p
  
}


#### Chl.a Monogynella ####
plot_list_Cuscuta_Chl.a[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Chl.a <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Chl.a_Monogynella <- with(data_Monogynella_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_Monogynella)
boxplot_positions_Chl.a_Monogynella <- as.data.frame(box.rslt_Chl.a_Monogynella$stats)

# what are these column tissue codes?
tissues_Chl.a_Monogynella <- levels(data_Monogynella_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_Monogynella) <- tissues_Chl.a_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_Monogynella <- boxplot_positions_Chl.a_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_Monogynella
# now join these positions to cbd
left_join(cbd_Chl.a_Monogynella, top_positions_Chl.a_Monogynella, by = "Tissue.code") -> cbd_Chl.a_Monogynella

# calculate how much to nudge
data_Monogynella_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_Monogynella
cbd_Chl.a_Monogynella$nudged <- max_Chl.a_Monogynella$max * 1.05


# add CLDs to plot
Chl.a_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_Monogynella_boxplot


Chl.a_Monogynella_boxplot


#### Chl.a Cuscuta ####
plot_list_Cuscuta_Chl.a[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_Cuscuta_boxplot


# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Chl.a <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Chl.a_Cuscuta <- with(data_Cuscuta_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_Cuscuta)
boxplot_positions_Chl.a_Cuscuta <- as.data.frame(box.rslt_Chl.a_Cuscuta$stats)

# what are these column tissue codes?
tissues_Chl.a_Cuscuta <- levels(data_Cuscuta_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_Cuscuta) <- tissues_Chl.a_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_Cuscuta <- boxplot_positions_Chl.a_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_Cuscuta
# now join these positions to cbd
left_join(cbd_Chl.a_Cuscuta, top_positions_Chl.a_Cuscuta, by = "Tissue.code") -> cbd_Chl.a_Cuscuta

# calculate how much to nudge
data_Cuscuta_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_Cuscuta
cbd_Chl.a_Cuscuta$nudged <- max_Chl.a_Cuscuta$max * 1.05


# add CLDs to plot
Chl.a_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_Cuscuta_boxplot


Chl.a_Cuscuta_boxplot




#### Chl.a Grammica ####
plot_list_Cuscuta_Chl.a[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_Grammica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_Grammica_Chl.a <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Chl.a_Grammica <- with(data_Grammica_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_Grammica)
boxplot_positions_Chl.a_Grammica <- as.data.frame(box.rslt_Chl.a_Grammica$stats)

# what are these column tissue codes?
tissues_Chl.a_Grammica <- levels(data_Grammica_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_Grammica) <- tissues_Chl.a_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_Grammica <- boxplot_positions_Chl.a_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_Grammica
# now join these positions to cbd
left_join(cbd_Chl.a_Grammica, top_positions_Chl.a_Grammica, by = "Tissue.code") -> cbd_Chl.a_Grammica

# calculate how much to nudge
data_Grammica_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_Grammica
cbd_Chl.a_Grammica$nudged <- (max_Chl.a_Grammica$max + 0.05) * 1.05

# add CLDs to plot
Chl.a_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_Grammica_boxplot


Chl.a_Grammica_boxplot

#### Chl.a C_purpurata alone ####

data_grammica_Chl.a <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == "C_purpurata")

Chl.a_C_purpurata_boxplot <- ggplot(data_grammica_Chl.a, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 8), breaks = seq(0, 8, by = 2)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 




# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_Chl.a <- dplyr::filter(data_Chl.a_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_Chl.a_C_purpurata <- with(data_C_purpurata_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_C_purpurata)
boxplot_positions_Chl.a_C_purpurata <- as.data.frame(box.rslt_Chl.a_C_purpurata$stats)

# what are these column tissue codes?
tissues_Chl.a_C_purpurata <- levels(data_C_purpurata_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_C_purpurata) <- tissues_Chl.a_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_C_purpurata <- boxplot_positions_Chl.a_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_C_purpurata
# now join these positions to cbd
left_join(cbd_Chl.a_C_purpurata, top_positions_Chl.a_C_purpurata, by = "Tissue.code") -> cbd_Chl.a_C_purpurata

# calculate how much to nudge
data_C_purpurata_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_C_purpurata
cbd_Chl.a_C_purpurata$nudged <- (max_Chl.a_C_purpurata$max + 0.01) * 1.05

# add CLDs to plot
Chl.a_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_C_purpurata_boxplot


Chl.a_C_purpurata_boxplot



#### Chl.b loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Chl.b
data_Chl.b_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Chl.b")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.b_plots_cuscutasub$Tissue.code <- factor(data_Chl.b_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Chl.b = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Chl.b_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 8), breaks = seq(0, 8, by = 2)) 
  plot_list_Cuscuta_Chl.b[[sub]] = p
  
}


#### Chl.b Monogynella ####
plot_list_Cuscuta_Chl.b[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.b" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_Monogynella_boxplot



# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Chl.b <- dplyr::filter(data_Chl.b_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Chl.b_Monogynella <- with(data_Monogynella_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_Monogynella)
boxplot_positions_Chl.b_Monogynella <- as.data.frame(box.rslt_Chl.b_Monogynella$stats)

# what are these column tissue codes?
tissues_Chl.b_Monogynella <- levels(data_Monogynella_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_Monogynella) <- tissues_Chl.b_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_Monogynella <- boxplot_positions_Chl.b_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_Monogynella
# now join these positions to cbd
left_join(cbd_Chl.b_Monogynella, top_positions_Chl.b_Monogynella, by = "Tissue.code") -> cbd_Chl.b_Monogynella

# calculate how much to nudge
data_Monogynella_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_Monogynella
cbd_Chl.b_Monogynella$nudged <- max_Chl.b_Monogynella$max * 1.05

# add CLDs to plot
Chl.b_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.b_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0))-> Chl.b_Monogynella_boxplot

Chl.b_Monogynella_boxplot


#### Chl.b Cuscuta ####
plot_list_Cuscuta_Chl.b[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.b" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_Cuscuta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Chl.b <- dplyr::filter(data_Chl.b_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Chl.b_Cuscuta <- with(data_Cuscuta_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_Cuscuta)
boxplot_positions_Chl.b_Cuscuta <- as.data.frame(box.rslt_Chl.b_Cuscuta$stats)

# what are these column tissue codes?
tissues_Chl.b_Cuscuta <- levels(data_Cuscuta_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_Cuscuta) <- tissues_Chl.b_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_Cuscuta <- boxplot_positions_Chl.b_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_Cuscuta
# now join these positions to cbd
left_join(cbd_Chl.b_Cuscuta, top_positions_Chl.b_Cuscuta, by = "Tissue.code") -> cbd_Chl.b_Cuscuta

# calculate how much to nudge
data_Cuscuta_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_Cuscuta
cbd_Chl.b_Cuscuta$nudged <- max_Chl.b_Cuscuta$max * 1.05

# add CLDs to plot
Chl.b_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.b_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_Cuscuta_boxplot


Chl.b_Cuscuta_boxplot


#### Chl.b Grammica ####
plot_list_Cuscuta_Chl.b[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.b" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Chl.b <- dplyr::filter(data_Chl.b_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Chl.b_Grammica <- with(data_Grammica_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_Grammica)
boxplot_positions_Chl.b_Grammica <- as.data.frame(box.rslt_Chl.b_Grammica$stats)

# what are these column tissue codes?
tissues_Chl.b_Grammica <- levels(data_Grammica_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_Grammica) <- tissues_Chl.b_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_Grammica <- boxplot_positions_Chl.b_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_Grammica
# now join these positions to cbd
left_join(cbd_Chl.b_Grammica, top_positions_Chl.b_Grammica, by = "Tissue.code") -> cbd_Chl.b_Grammica

# calculate how much to nudge
data_Grammica_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_Grammica
cbd_Chl.b_Grammica$nudged <- max_Chl.b_Grammica$max * 1.05

# add CLDs to plot
Chl.b_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.b_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_Grammica_boxplot


Chl.b_Grammica_boxplot



#### Chl.b C_purpurata alone ####

data_C_purpurata_Chl.b <- dplyr::filter(data_Chl.b_plots_cuscutasub, Subgenus == "C_purpurata")

Chl.b_C_purpurata_boxplot <- ggplot(data_C_purpurata_Chl.b, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 8), breaks = seq(0, 8, by = 2)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.b" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE)

Chl.b_C_purpurata_boxplot





#### Tot.Chl loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Tot.Chl
data_Tot.Chl_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Tot.Chl")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Tot.Chl_plots_cuscutasub$Tissue.code <- factor(data_Tot.Chl_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Tot.Chl = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 10), breaks = seq(0, 10, by = 2))
  plot_list_Cuscuta_Tot.Chl[[sub]] = p
  
}

#### Tot.Chl. Monogynella ####
plot_list_Cuscuta_Tot.Chl[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Chl" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Tot.Chl_Monogynella <- with(data_Monogynella_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_Monogynella)
boxplot_positions_Tot.Chl_Monogynella <- as.data.frame(box.rslt_Tot.Chl_Monogynella$stats)

# what are these column tissue codes?
tissues_Tot.Chl_Monogynella <- levels(data_Monogynella_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_Monogynella) <- tissues_Tot.Chl_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_Monogynella <- boxplot_positions_Tot.Chl_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_Monogynella
# now join these positions to cbd
left_join(cbd_Tot.Chl_Monogynella, top_positions_Tot.Chl_Monogynella, by = "Tissue.code") -> cbd_Tot.Chl_Monogynella

# calculate how much to nudge
data_Monogynella_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_Monogynella
cbd_Tot.Chl_Monogynella$nudged <- max_Tot.Chl_Monogynella$max * 1.05


# add CLDs to plot
Tot.Chl_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Chl_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_Monogynella_boxplot


Tot.Chl_Monogynella_boxplot


#### Tot.Chl Cuscuta ####
plot_list_Cuscuta_Tot.Chl[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Chl" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Tot.Chl_Cuscuta <- with(data_Cuscuta_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_Cuscuta)
boxplot_positions_Tot.Chl_Cuscuta <- as.data.frame(box.rslt_Tot.Chl_Cuscuta$stats)

# what are these column tissue codes?
tissues_Tot.Chl_Cuscuta <- levels(data_Cuscuta_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_Cuscuta) <- tissues_Tot.Chl_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_Cuscuta <- boxplot_positions_Tot.Chl_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_Cuscuta
# now join these positions to cbd
left_join(cbd_Tot.Chl_Cuscuta, top_positions_Tot.Chl_Cuscuta, by = "Tissue.code") -> cbd_Tot.Chl_Cuscuta

# calculate how much to nudge
data_Cuscuta_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_Cuscuta
cbd_Tot.Chl_Cuscuta$nudged <- max_Tot.Chl_Cuscuta$max * 1.05


# add CLDs to plot
Tot.Chl_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Chl_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_Cuscuta_boxplot


Tot.Chl_Cuscuta_boxplot



#### Tot.Chl Grammica ####
plot_list_Cuscuta_Tot.Chl[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Chl" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Tot.Chl_Grammica <- with(data_Grammica_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_Grammica)
boxplot_positions_Tot.Chl_Grammica <- as.data.frame(box.rslt_Tot.Chl_Grammica$stats)

# what are these column tissue codes?
tissues_Tot.Chl_Grammica <- levels(data_Grammica_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_Grammica) <- tissues_Tot.Chl_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_Grammica <- boxplot_positions_Tot.Chl_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_Grammica
# now join these positions to cbd
left_join(cbd_Tot.Chl_Grammica, top_positions_Tot.Chl_Grammica, by = "Tissue.code") -> cbd_Tot.Chl_Grammica

# calculate how much to nudge
data_Grammica_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_Grammica
cbd_Tot.Chl_Grammica$nudged <- (max_Tot.Chl_Grammica$max + .05) * 1.05

# add CLDs to plot
Tot.Chl_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Chl_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_Grammica_boxplot


Tot.Chl_Grammica_boxplot




#### Tot.Chl C_purpurata alone ####

data_C_purpurata_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == "C_purpurata")

Tot.Chl_C_purpurata_boxplot <- ggplot(data_C_purpurata_Tot.Chl, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 10), breaks = seq(0, 10, by = 2)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Chl" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

Tot.Chl_C_purpurata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_Tot.Chl_C_purpurata <- with(data_C_purpurata_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_purpurata)
boxplot_positions_Tot.Chl_C_purpurata <- as.data.frame(box.rslt_Tot.Chl_C_purpurata$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_purpurata <- levels(data_C_purpurata_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_purpurata) <- tissues_Tot.Chl_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_purpurata <- boxplot_positions_Tot.Chl_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_purpurata
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_purpurata, top_positions_Tot.Chl_C_purpurata, by = "Tissue.code") -> cbd_Tot.Chl_C_purpurata

# calculate how much to nudge
data_C_purpurata_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_purpurata
cbd_Tot.Chl_C_purpurata$nudged <- max_Tot.Chl_C_purpurata$max * 1.05

# add CLDs to plot
Tot.Chl_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Chl_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_purpurata_boxplot


Tot.Chl_C_purpurata_boxplot





#### Chl.a.b loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Chl.a.b
data_Chl.a.b_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Chl.a.b")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.a.b_plots_cuscutasub$Tissue.code <- factor(data_Chl.a.b_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Chl.a.b = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Chl.a.b_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 10), breaks = seq(0, 10, by = 2))
  plot_list_Cuscuta_Chl.a.b[[sub]] = p
  
}

#### Chl.a.b Monogynella ####
plot_list_Cuscuta_Chl.a.b[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a.b" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Chl.a.b_Monogynella <- with(data_Monogynella_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_Monogynella)
boxplot_positions_Chl.a.b_Monogynella <- as.data.frame(box.rslt_Chl.a.b_Monogynella$stats)

# what are these column tissue codes?
tissues_Chl.a.b_Monogynella <- levels(data_Monogynella_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_Monogynella) <- tissues_Chl.a.b_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_Monogynella <- boxplot_positions_Chl.a.b_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_Monogynella
# now join these positions to cbd
left_join(cbd_Chl.a.b_Monogynella, top_positions_Chl.a.b_Monogynella, by = "Tissue.code") -> cbd_Chl.a.b_Monogynella

# calculate how much to nudge
data_Monogynella_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_Monogynella
cbd_Chl.a.b_Monogynella$nudged <- (max_Chl.a.b_Monogynella$max + 0.00) * 1.05

# add CLDs to plot
Chl.a.b_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a.b_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_Monogynella_boxplot


Chl.a.b_Monogynella_boxplot



#### Chl.a.b Cuscuta ####
plot_list_Cuscuta_Chl.a.b[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a.b" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Chl.a.b_Cuscuta <- with(data_Cuscuta_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_Cuscuta)
boxplot_positions_Chl.a.b_Cuscuta <- as.data.frame(box.rslt_Chl.a.b_Cuscuta$stats)

# what are these column tissue codes?
tissues_Chl.a.b_Cuscuta <- levels(data_Cuscuta_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_Cuscuta) <- tissues_Chl.a.b_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_Cuscuta <- boxplot_positions_Chl.a.b_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_Cuscuta
# now join these positions to cbd
left_join(cbd_Chl.a.b_Cuscuta, top_positions_Chl.a.b_Cuscuta, by = "Tissue.code") -> cbd_Chl.a.b_Cuscuta

# calculate how much to nudge
data_Cuscuta_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_Cuscuta
cbd_Chl.a.b_Cuscuta$nudged <- (max_Chl.a.b_Cuscuta$max + 0.00) * 1.05

# add CLDs to plot
Chl.a.b_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a.b_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_Cuscuta_boxplot


Chl.a.b_Cuscuta_boxplot



#### Chl.a.b Grammica ####
plot_list_Cuscuta_Chl.a.b[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Chl.a.b" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Chl.a.b_Grammica <- with(data_Grammica_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_Grammica)
boxplot_positions_Chl.a.b_Grammica <- as.data.frame(box.rslt_Chl.a.b_Grammica$stats)

# what are these column tissue codes?
tissues_Chl.a.b_Grammica <- levels(data_Grammica_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_Grammica) <- tissues_Chl.a.b_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_Grammica <- boxplot_positions_Chl.a.b_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_Grammica
# now join these positions to cbd
left_join(cbd_Chl.a.b_Grammica, top_positions_Chl.a.b_Grammica, by = "Tissue.code") -> cbd_Chl.a.b_Grammica

# calculate how much to nudge
data_Grammica_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_Grammica
cbd_Chl.a.b_Grammica$nudged <- (max_Chl.a.b_Grammica$max + 0.00) * 1.05

# add CLDs to plot
Chl.a.b_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Chl.a.b_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_Grammica_boxplot


Chl.a.b_Grammica_boxplot




#### Chl.a.b C_purpurata alone ####
data_C_purpurata_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_cuscutasub, Subgenus == "C_purpurata")

Chl.a.b_C_purpurata_boxplot <- ggplot(data_C_purpurata_Chl.a.b, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Mass ratio") +
  scale_y_continuous(position = "right", limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  guides(colour = guide_legend(nrow = 1)) 

Chl.a.b_C_purpurata_boxplot
# no KW, only one sample





#### ####
#### COMBINED (faceted) chlorophyll plot: subgenus ####
wrap_elements(gridtext::richtext_grob('Total chlorophyll', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Tot.Chl_ipomoea_boxplot + ggtitle('Ipomoea nil') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Tot.Chl_Monogynella_boxplot + ggtitle('Monogynella') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Tot.Chl_Cuscuta_boxplot + ggtitle('Cuscuta') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Tot.Chl_Grammica_boxplot + ggtitle('Grammica') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) +
  Tot.Chl_C_purpurata_boxplot + ggtitle('C. purpurata') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  wrap_elements(gridtext::richtext_grob('Chlorophyll *a*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.a_ipomoea_boxplot + 
  Chl.a_Monogynella_boxplot +  
  Chl.a_Cuscuta_boxplot + 
  Chl.a_Grammica_boxplot + 
  Chl.a_C_purpurata_boxplot + 
  wrap_elements(gridtext::richtext_grob('Chlorophyll *b*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.b_ipomoea_boxplot +
  Chl.b_Monogynella_boxplot +  
  Chl.b_Cuscuta_boxplot + 
  Chl.b_Grammica_boxplot + 
  Chl.b_C_purpurata_boxplot + 
  wrap_elements(gridtext::richtext_grob('Chlorophyll *a*:*b*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.a.b_ipomoea_boxplot + 
  Chl.a.b_Monogynella_boxplot +  
  Chl.a.b_Cuscuta_boxplot + 
  Chl.a.b_Grammica_boxplot + 
  Chl.a.b_C_purpurata_boxplot + plot_layout(nrow = 4, byrow = T) -> chlorophyll_boxplot_new

chlorophyll_boxplot_new 

pdf("../output/boxplots/chlorophyll_boxplot_new.pdf", width=7,height=7) 
chlorophyll_boxplot_new
dev.off()


#### CAROTENOID PLOTS ####


#### VAZ.Car Ipomoea_nil ####
data_ipomoea <-  dplyr::filter(data_long_calcs, Subgenus == "Ipomoea_nil")
data_ipomoea_VAZ.Car <- dplyr::filter(data_ipomoea, Pigment == "VAZ.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_ipomoea$Tissue.code <- factor(data_ipomoea$Tissue.code, levels = c("l", "y", "o", "f", "s"))

VAZ.Car_ipomoea_boxplot <- ggplot(data_ipomoea_VAZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) +
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "VAZ.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

VAZ.Car_ipomoea_boxplot 

# use base R boxplot to get the coordinates of the boxes
box.rslt_VAZ.Car_ipomoea <- with(data_ipomoea_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_ipomoea)
boxplot_positions_VAZ.Car_ipomoea <- as.data.frame(box.rslt_VAZ.Car_ipomoea$stats)

# what are these column tissue codes?
tissues_VAZ.Car_ipomoea <- levels(data_ipomoea_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_ipomoea) <- tissues_VAZ.Car_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_ipomoea <- boxplot_positions_VAZ.Car_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_ipomoea
# now join these positions to cbd
left_join(cbd_VAZ.Car_ipomoea, top_positions_VAZ.Car_ipomoea, by = "Tissue.code") -> cbd_VAZ.Car_ipomoea

# calculate how much to nudge
data_ipomoea_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_ipomoea
cbd_VAZ.Car_ipomoea$nudged <- (max_VAZ.Car_ipomoea$max + 0.01) * 1.05


# add CLDs to plot
VAZ.Car_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_VAZ.Car_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_ipomoea_boxplot

VAZ.Car_ipomoea_boxplot


#### Neoxanthin Ipomoea_nil ####
data_ipomoea_Neoxanthin <- dplyr::filter(data_ipomoea, Pigment == "Neo.Car")


Neoxanthin_ipomoea_boxplot <- ggplot(data_ipomoea_Neoxanthin, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100))+
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Neo.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
box.rslt_Neoxanthin_ipomoea <- with(data_ipomoea_Neoxanthin, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Neoxanthin_ipomoea)
boxplot_positions_Neoxanthin_ipomoea <- as.data.frame(box.rslt_Neoxanthin_ipomoea$stats)

# what are these column tissue codes?
tissues_Neoxanthin_ipomoea <- levels(data_ipomoea_Neoxanthin$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Neoxanthin_ipomoea) <- tissues_Neoxanthin_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Neoxanthin_ipomoea <- boxplot_positions_Neoxanthin_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Neoxanthin_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Neo.Car"]])[["Letters"]])
colnames(cbd_Neoxanthin_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Neoxanthin_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Neoxanthin_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Neoxanthin_ipomoea
# now join these positions to cbd
left_join(cbd_Neoxanthin_ipomoea, top_positions_Neoxanthin_ipomoea, by = "Tissue.code") -> cbd_Neoxanthin_ipomoea

# calculate how much to nudge
data_ipomoea_Neoxanthin %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Neoxanthin_ipomoea
cbd_Neoxanthin_ipomoea$nudged <- max_Neoxanthin_ipomoea$max * 1.05



# add CLDs to plot
Neoxanthin_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Neoxanthin_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Neoxanthin_ipomoea_boxplot


Neoxanthin_ipomoea_boxplot


#### Lutein.epoxide Ipomoea_nil ####
data_ipomoea_Lutein.epoxide <- dplyr::filter(data_ipomoea, Pigment == "Lut.epo.Car")


Lutein.epoxide_ipomoea_boxplot <- ggplot(data_ipomoea_Lutein.epoxide, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) +     
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.epo.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# KW not significant so no post hoc

Lutein.epoxide_ipomoea_boxplot



#### Lutein Ipomoea_nil ####
data_ipomoea_Lutein <- dplyr::filter(data_ipomoea, Pigment == "Lut.Car")


Lutein_ipomoea_boxplot <- ggplot(data_ipomoea_Lutein, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_Lutein_ipomoea <- with(data_ipomoea_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_ipomoea)
boxplot_positions_Lutein_ipomoea <- as.data.frame(box.rslt_Lutein_ipomoea$stats)

# what are these column tissue codes?
tissues_Lutein_ipomoea <- levels(data_ipomoea_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_ipomoea) <- tissues_Lutein_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_ipomoea <- boxplot_positions_Lutein_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_ipomoea
# now join these positions to cbd
left_join(cbd_Lutein_ipomoea, top_positions_Lutein_ipomoea, by = "Tissue.code") -> cbd_Lutein_ipomoea

# calculate how much to nudge
data_ipomoea_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_ipomoea
cbd_Lutein_ipomoea$nudged <- max_Lutein_ipomoea$max * 1.05


# add CLDs to plot
Lutein_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_ipomoea_boxplot


Lutein_ipomoea_boxplot




#### a.Carotene Ipomoea_nil ####
data_ipomoea_a.Carotene <- dplyr::filter(data_ipomoea, Pigment == "a.Car.Car")


a.Carotene_ipomoea_boxplot <- ggplot(data_ipomoea_a.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "a.Car.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_a.Carotene_ipomoea <- with(data_ipomoea_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_ipomoea)
boxplot_positions_a.Carotene_ipomoea <- as.data.frame(box.rslt_a.Carotene_ipomoea$stats)

# what are these column tissue codes?
tissues_a.Carotene_ipomoea <- levels(data_ipomoea_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_ipomoea) <- tissues_a.Carotene_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_ipomoea <- boxplot_positions_a.Carotene_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_ipomoea
# now join these positions to cbd
left_join(cbd_a.Carotene_ipomoea, top_positions_a.Carotene_ipomoea, by = "Tissue.code") -> cbd_a.Carotene_ipomoea

# calculate how much to nudge
data_ipomoea_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_ipomoea
cbd_a.Carotene_ipomoea$nudged <-( max_a.Carotene_ipomoea$max + 0.01) * 1.05


# add CLDs to plot
a.Carotene_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_a.Carotene_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_ipomoea_boxplot


a.Carotene_ipomoea_boxplot



#### b.Carotene Ipomoea_nil ####
data_ipomoea_b.Carotene <- dplyr::filter(data_ipomoea, Pigment == "b.Car.Car")


b.Carotene_ipomoea_boxplot <- ggplot(data_ipomoea_b.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "b.Car.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_b.Carotene_ipomoea <- with(data_ipomoea_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_ipomoea)
boxplot_positions_b.Carotene_ipomoea <- as.data.frame(box.rslt_b.Carotene_ipomoea$stats)

# what are these column tissue codes?
tissues_b.Carotene_ipomoea <- levels(data_ipomoea_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_ipomoea) <- tissues_b.Carotene_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_ipomoea <- boxplot_positions_b.Carotene_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_ipomoea
# now join these positions to cbd
left_join(cbd_b.Carotene_ipomoea, top_positions_b.Carotene_ipomoea, by = "Tissue.code") -> cbd_b.Carotene_ipomoea

# calculate how much to nudge
data_ipomoea_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_ipomoea
cbd_b.Carotene_ipomoea$nudged <- max_b.Carotene_ipomoea$max * 1.05


# add CLDs to plot
b.Carotene_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_b.Carotene_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_ipomoea_boxplot


b.Carotene_ipomoea_boxplot


#### Tot.Car Ipomoea_nil ####
data_ipomoea_Tot.Car <- dplyr::filter(data_ipomoea, Pigment == "Tot.Car")
# this one is normalized by FW

Tot.Car_ipomoea_boxplot <- ggplot(data_ipomoea_Tot.Car, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 8)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_Tot.Car_ipomoea <- with(data_ipomoea_Tot.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_ipomoea)
boxplot_positions_Tot.Car_ipomoea <- as.data.frame(box.rslt_Tot.Car_ipomoea$stats)

# what are these column tissue codes?
tissues_Tot.Car_ipomoea <- levels(data_ipomoea_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_ipomoea) <- tissues_Tot.Car_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_ipomoea <- boxplot_positions_Tot.Car_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_ipomoea
# now join these positions to cbd
left_join(cbd_Tot.Car_ipomoea, top_positions_Tot.Car_ipomoea, by = "Tissue.code") -> cbd_Tot.Car_ipomoea

# calculate how much to nudge
data_ipomoea_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_ipomoea
cbd_Tot.Car_ipomoea$nudged <- max_Tot.Car_ipomoea$max * 1.05


# add CLDs to plot
Tot.Car_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Car_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_ipomoea_boxplot


Tot.Car_ipomoea_boxplot


#### NVZ.Car Ipomoea_nil plot####
data_ipomoea_NVZ.Car <- dplyr::filter(data_ipomoea, Pigment == "NVZ.Car")


NVZ.Car_ipomoea_boxplot <- ggplot(data_ipomoea_NVZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  # ylab("% by ng/mg Fresh Weight") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100))+
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "NVZ.Car" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_NVZ.Car_ipomoea <- with(data_ipomoea_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_ipomoea)
boxplot_positions_NVZ.Car_ipomoea <- as.data.frame(box.rslt_NVZ.Car_ipomoea$stats)

# what are these column tissue codes?
tissues_NVZ.Car_ipomoea <- levels(data_ipomoea_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_ipomoea) <- tissues_NVZ.Car_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_ipomoea <- boxplot_positions_NVZ.Car_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_ipomoea
# now join these positions to cbd
left_join(cbd_NVZ.Car_ipomoea, top_positions_NVZ.Car_ipomoea, by = "Tissue.code") -> cbd_NVZ.Car_ipomoea

# calculate how much to nudge
data_ipomoea_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_ipomoea
cbd_NVZ.Car_ipomoea$nudged <- max_NVZ.Car_ipomoea$max * 1.05


# add CLDs to plot
NVZ.Car_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NVZ.Car_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_ipomoea_boxplot


NVZ.Car_ipomoea_boxplot


#### VAZ.Car loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for VAZ.Car
data_VAZ.Car_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "VAZ.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_VAZ.Car_plots_cuscutasub$Tissue.code <- factor(data_VAZ.Car_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_VAZ.Car = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "left", limits = c(0, 100))  
  plot_list_Cuscuta_VAZ.Car[[sub]] = p
  
}


#### VAZ.Car Monogynella ####
plot_list_Cuscuta_VAZ.Car[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "VAZ.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_VAZ.Car_Monogynella <- with(data_Monogynella_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_Monogynella)
boxplot_positions_VAZ.Car_Monogynella <- as.data.frame(box.rslt_VAZ.Car_Monogynella$stats)

# what are these column tissue codes?
tissues_VAZ.Car_Monogynella <- levels(data_Monogynella_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_Monogynella) <- tissues_VAZ.Car_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_Monogynella <- boxplot_positions_VAZ.Car_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_Monogynella
# now join these positions to cbd
left_join(cbd_VAZ.Car_Monogynella, top_positions_VAZ.Car_Monogynella, by = "Tissue.code") -> cbd_VAZ.Car_Monogynella

# calculate how much to nudge
data_Monogynella_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_Monogynella
cbd_VAZ.Car_Monogynella$nudged <- (max_VAZ.Car_Monogynella$max + 0.01) * 1.05


# add CLDs to plot
VAZ.Car_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_VAZ.Car_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_Monogynella_boxplot


VAZ.Car_Monogynella_boxplot


#### VAZ.Car Cuscuta ####
plot_list_Cuscuta_VAZ.Car[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "VAZ.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_Cuscuta_boxplot


# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_VAZ.Car_Cuscuta <- with(data_Cuscuta_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_Cuscuta)
boxplot_positions_VAZ.Car_Cuscuta <- as.data.frame(box.rslt_VAZ.Car_Cuscuta$stats)

# what are these column tissue codes?
tissues_VAZ.Car_Cuscuta <- levels(data_Cuscuta_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_Cuscuta) <- tissues_VAZ.Car_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_Cuscuta <- boxplot_positions_VAZ.Car_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_Cuscuta
# now join these positions to cbd
left_join(cbd_VAZ.Car_Cuscuta, top_positions_VAZ.Car_Cuscuta, by = "Tissue.code") -> cbd_VAZ.Car_Cuscuta

# calculate how much to nudge
data_Cuscuta_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_Cuscuta
cbd_VAZ.Car_Cuscuta$nudged <- max_VAZ.Car_Cuscuta$max  * 1.05


# add CLDs to plot
VAZ.Car_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_VAZ.Car_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_Cuscuta_boxplot


VAZ.Car_Cuscuta_boxplot




#### VAZ.Car Grammica ####
plot_list_Cuscuta_VAZ.Car[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "VAZ.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_Grammica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_Grammica_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_VAZ.Car_Grammica <- with(data_Grammica_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_Grammica)
boxplot_positions_VAZ.Car_Grammica <- as.data.frame(box.rslt_VAZ.Car_Grammica$stats)

# what are these column tissue codes?
tissues_VAZ.Car_Grammica <- levels(data_Grammica_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_Grammica) <- tissues_VAZ.Car_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_Grammica <- boxplot_positions_VAZ.Car_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_Grammica
# now join these positions to cbd
left_join(cbd_VAZ.Car_Grammica, top_positions_VAZ.Car_Grammica, by = "Tissue.code") -> cbd_VAZ.Car_Grammica

# calculate how much to nudge
data_Grammica_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_Grammica
cbd_VAZ.Car_Grammica$nudged <- max_VAZ.Car_Grammica$max  * 1.05

# add CLDs to plot
VAZ.Car_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_VAZ.Car_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_Grammica_boxplot


VAZ.Car_Grammica_boxplot

#### VAZ.Car C_purpurata alone ####

data_grammica_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == "C_purpurata")

VAZ.Car_C_purpurata_boxplot <- ggplot(data_grammica_VAZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "VAZ.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 




# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_VAZ.Car_C_purpurata <- with(data_C_purpurata_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_purpurata)
boxplot_positions_VAZ.Car_C_purpurata <- as.data.frame(box.rslt_VAZ.Car_C_purpurata$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_purpurata <- levels(data_C_purpurata_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_purpurata) <- tissues_VAZ.Car_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_purpurata <- boxplot_positions_VAZ.Car_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_purpurata
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_purpurata, top_positions_VAZ.Car_C_purpurata, by = "Tissue.code") -> cbd_VAZ.Car_C_purpurata

# calculate how much to nudge
data_C_purpurata_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_purpurata
cbd_VAZ.Car_C_purpurata$nudged <- (max_VAZ.Car_C_purpurata$max + 0.01) * 1.05

# add CLDs to plot
VAZ.Car_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_VAZ.Car_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_purpurata_boxplot


VAZ.Car_C_purpurata_boxplot



#### Neoxanthin loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Neoxanthin
data_Neoxanthin_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Neo.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Neoxanthin_plots_cuscutasub$Tissue.code <- factor(data_Neoxanthin_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Neoxanthin = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Neoxanthin_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100)) 
  plot_list_Cuscuta_Neoxanthin[[sub]] = p
  
}


#### Neoxanthin Monogynella ####
plot_list_Cuscuta_Neoxanthin[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Neo.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Neoxanthin_Monogynella_boxplot



# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Neoxanthin <- dplyr::filter(data_Neoxanthin_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Neoxanthin_Monogynella <- with(data_Monogynella_Neoxanthin, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Neoxanthin_Monogynella)
boxplot_positions_Neoxanthin_Monogynella <- as.data.frame(box.rslt_Neoxanthin_Monogynella$stats)

# what are these column tissue codes?
tissues_Neoxanthin_Monogynella <- levels(data_Monogynella_Neoxanthin$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Neoxanthin_Monogynella) <- tissues_Neoxanthin_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Neoxanthin_Monogynella <- boxplot_positions_Neoxanthin_Monogynella[5,]


# # add pairwise significance letter groups (compact letter display; CLD)
# wilcox_list[["Monogynella__Neo.Car"]]
# cbd_Neoxanthin_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Neo.Car"]])[["Letters"]])
# colnames(cbd_Neoxanthin_Monogynella)[1] <- "Letter"
# # turn rownames into first column for Tissue.code
# setDT(cbd_Neoxanthin_Monogynella, keep.rownames = "Tissue.code")


# Monogynella Neoxanthin was causing a problem in this pipeline because of the exact matches.. 
wilcox_list[["Monogynella__Neo.Car"]]
# QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Neo.Car"]])

# create letter groups manually 

# # Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  data_loop_wilcox$FW.norm and data_loop_wilcox$Tissue.code 
# 
#     sdlg   y      o h
#   y 0.0014 -      - -
#   o 0.0014 0.0014 - -
#   h 0.0038 0.0038 - -
#   f 0.0446 0.0446 - -
#   
#   P value adjustment method: BH 

# sdlg diff from everyone
# sdlg a 
# y diff from o, h, f
# y = b
# o, h, f not compared but exact match = c

# make df
cbd_Neoxanthin_Monogynella_Tissue.code <- c("sdlg", "y", "o", "h", "f")
cbd_Neoxanthin_Monogynella_Letters <- c("a", "b", "c", "c", "c")
cbd_Neoxanthin_Monogynella <- data.frame("Tissue.code" = cbd_Neoxanthin_Monogynella_Tissue.code, "Letter" = cbd_Neoxanthin_Monogynella_Letters)

# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Neoxanthin_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Neoxanthin_Monogynella
# now join these positions to cbd
left_join(cbd_Neoxanthin_Monogynella, top_positions_Neoxanthin_Monogynella, by = "Tissue.code") -> cbd_Neoxanthin_Monogynella

# calculate how much to nudge
data_Monogynella_Neoxanthin %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Neoxanthin_Monogynella
cbd_Neoxanthin_Monogynella$nudged <- (max_Neoxanthin_Monogynella$max + 0.01) * 1.05

# add CLDs to plot
Neoxanthin_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Neoxanthin_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0))-> Neoxanthin_Monogynella_boxplot

Neoxanthin_Monogynella_boxplot


#### Neoxanthin Cuscuta ####
plot_list_Cuscuta_Neoxanthin[["Cuscuta"]]  -> Neoxanthin_Cuscuta_boxplot


# no KW and no posthoc

Neoxanthin_Cuscuta_boxplot


#### Neoxanthin Grammica ####
# no KW and no Post hoc
plot_list_Cuscuta_Neoxanthin[["Grammica"]]+ 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Neo.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE)  -> Neoxanthin_Grammica_boxplot





#### Neoxanthin C_purpurata alone ####

data_C_purpurata_Neoxanthin <- dplyr::filter(data_Neoxanthin_plots_cuscutasub, Subgenus == "C_purpurata")

Neoxanthin_C_purpurata_boxplot <- ggplot(data_C_purpurata_Neoxanthin, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) 
Neoxanthin_C_purpurata_boxplot

# not sig





#### Lutein.epoxide loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Lutein.epoxide
data_Lutein.epoxide_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Lut.epo.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Lutein.epoxide_plots_cuscutasub$Tissue.code <- factor(data_Lutein.epoxide_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Lutein.epoxide = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Lutein.epoxide_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Cuscuta_Lutein.epoxide[[sub]] = p
  
}

#### Lutein.epoxide. Monogynella ####
plot_list_Cuscuta_Lutein.epoxide[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.epo.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_Monogynella_boxplot

Lutein.epoxide_Monogynella_boxplot
# not sig

#### Lutein.epoxide Cuscuta ####
plot_list_Cuscuta_Lutein.epoxide[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.epo.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Lutein.epoxide_Cuscuta <- with(data_Cuscuta_Lutein.epoxide, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein.epoxide_Cuscuta)
boxplot_positions_Lutein.epoxide_Cuscuta <- as.data.frame(box.rslt_Lutein.epoxide_Cuscuta$stats)

# what are these column tissue codes?
tissues_Lutein.epoxide_Cuscuta <- levels(data_Cuscuta_Lutein.epoxide$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein.epoxide_Cuscuta) <- tissues_Lutein.epoxide_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein.epoxide_Cuscuta <- boxplot_positions_Lutein.epoxide_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein.epoxide_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Lut.epo.Car"]])[["Letters"]])
colnames(cbd_Lutein.epoxide_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein.epoxide_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein.epoxide_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Lutein.epoxide_Cuscuta
# now join these positions to cbd
left_join(cbd_Lutein.epoxide_Cuscuta, top_positions_Lutein.epoxide_Cuscuta, by = "Tissue.code") -> cbd_Lutein.epoxide_Cuscuta

# calculate how much to nudge
data_Cuscuta_Lutein.epoxide %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein.epoxide_Cuscuta
cbd_Lutein.epoxide_Cuscuta$nudged <- (max_Lutein.epoxide_Cuscuta$max + 0.01) * 1.05


# add CLDs to plot
Lutein.epoxide_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein.epoxide_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein.epoxide_Cuscuta_boxplot


Lutein.epoxide_Cuscuta_boxplot



#### Lutein.epoxide Grammica ####
plot_list_Cuscuta_Lutein.epoxide[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.epo.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Lutein.epoxide_Grammica <- with(data_Grammica_Lutein.epoxide, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein.epoxide_Grammica)
boxplot_positions_Lutein.epoxide_Grammica <- as.data.frame(box.rslt_Lutein.epoxide_Grammica$stats)

# what are these column tissue codes?
tissues_Lutein.epoxide_Grammica <- levels(data_Grammica_Lutein.epoxide$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein.epoxide_Grammica) <- tissues_Lutein.epoxide_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein.epoxide_Grammica <- boxplot_positions_Lutein.epoxide_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein.epoxide_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Lut.epo.Car"]])[["Letters"]])
colnames(cbd_Lutein.epoxide_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein.epoxide_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein.epoxide_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Lutein.epoxide_Grammica
# now join these positions to cbd
left_join(cbd_Lutein.epoxide_Grammica, top_positions_Lutein.epoxide_Grammica, by = "Tissue.code") -> cbd_Lutein.epoxide_Grammica

# calculate how much to nudge
data_Grammica_Lutein.epoxide %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein.epoxide_Grammica
cbd_Lutein.epoxide_Grammica$nudged <- (max_Lutein.epoxide_Grammica$max + 0.01) * 1.05


# add CLDs to plot
Lutein.epoxide_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein.epoxide_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein.epoxide_Grammica_boxplot

Lutein.epoxide_Grammica_boxplot


#### Lutein.epoxide C_purpurata alone ####

data_C_purpurata_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_cuscutasub, Subgenus == "C_purpurata")

Lutein.epoxide_C_purpurata_boxplot <- ggplot(data_C_purpurata_Lutein.epoxide, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.epo.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

Lutein.epoxide_C_purpurata_boxplot

# not sig


#### Lutein loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Lutein
data_Lutein_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Lut.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Lutein_plots_cuscutasub$Tissue.code <- factor(data_Lutein_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Lutein = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Cuscuta_Lutein[[sub]] = p
  
}

#### Lutein. Monogynella ####
plot_list_Cuscuta_Lutein[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Lutein <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Lutein_Monogynella <- with(data_Monogynella_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_Monogynella)
boxplot_positions_Lutein_Monogynella <- as.data.frame(box.rslt_Lutein_Monogynella$stats)

# what are these column tissue codes?
tissues_Lutein_Monogynella <- levels(data_Monogynella_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_Monogynella) <- tissues_Lutein_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_Monogynella <- boxplot_positions_Lutein_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_Monogynella
# now join these positions to cbd
left_join(cbd_Lutein_Monogynella, top_positions_Lutein_Monogynella, by = "Tissue.code") -> cbd_Lutein_Monogynella

# calculate how much to nudge
data_Monogynella_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_Monogynella
cbd_Lutein_Monogynella$nudged <- (max_Lutein_Monogynella$max + 0.01) * 1.05


# add CLDs to plot
Lutein_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_Monogynella_boxplot


Lutein_Monogynella_boxplot


#### Lutein Cuscuta ####
plot_list_Cuscuta_Lutein[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Lutein <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Lutein_Cuscuta <- with(data_Cuscuta_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_Cuscuta)
boxplot_positions_Lutein_Cuscuta <- as.data.frame(box.rslt_Lutein_Cuscuta$stats)

# what are these column tissue codes?
tissues_Lutein_Cuscuta <- levels(data_Cuscuta_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_Cuscuta) <- tissues_Lutein_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_Cuscuta <- boxplot_positions_Lutein_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_Cuscuta
# now join these positions to cbd
left_join(cbd_Lutein_Cuscuta, top_positions_Lutein_Cuscuta, by = "Tissue.code") -> cbd_Lutein_Cuscuta

# calculate how much to nudge
data_Cuscuta_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_Cuscuta
cbd_Lutein_Cuscuta$nudged <- (max_Lutein_Cuscuta$max + 0.01) * 1.05


# add CLDs to plot
Lutein_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_Cuscuta_boxplot


Lutein_Cuscuta_boxplot



#### Lutein Grammica ####
plot_list_Cuscuta_Lutein[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Lutein <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Lutein_Grammica <- with(data_Grammica_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_Grammica)
boxplot_positions_Lutein_Grammica <- as.data.frame(box.rslt_Lutein_Grammica$stats)

# what are these column tissue codes?
tissues_Lutein_Grammica <- levels(data_Grammica_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_Grammica) <- tissues_Lutein_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_Grammica <- boxplot_positions_Lutein_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_Grammica
# now join these positions to cbd
left_join(cbd_Lutein_Grammica, top_positions_Lutein_Grammica, by = "Tissue.code") -> cbd_Lutein_Grammica

# calculate how much to nudge
data_Grammica_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_Grammica
cbd_Lutein_Grammica$nudged <- (max_Lutein_Grammica$max + 0.01) * 1.05

# add CLDs to plot
Lutein_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_Grammica_boxplot


Lutein_Grammica_boxplot




#### Lutein C_purpurata alone ####

data_C_purpurata_Lutein <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == "C_purpurata")

Lutein_C_purpurata_boxplot <- ggplot(data_C_purpurata_Lutein, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Lut.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

Lutein_C_purpurata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_Lutein <- dplyr::filter(data_Lutein_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_Lutein_C_purpurata <- with(data_C_purpurata_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_purpurata)
boxplot_positions_Lutein_C_purpurata <- as.data.frame(box.rslt_Lutein_C_purpurata$stats)

# what are these column tissue codes?
tissues_Lutein_C_purpurata <- levels(data_C_purpurata_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_purpurata) <- tissues_Lutein_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_purpurata <- boxplot_positions_Lutein_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_purpurata
# now join these positions to cbd
left_join(cbd_Lutein_C_purpurata, top_positions_Lutein_C_purpurata, by = "Tissue.code") -> cbd_Lutein_C_purpurata

# calculate how much to nudge
data_C_purpurata_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_purpurata
cbd_Lutein_C_purpurata$nudged <- (max_Lutein_C_purpurata$max + 0.01) * 1.05

# add CLDs to plot
Lutein_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Lutein_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_purpurata_boxplot


Lutein_C_purpurata_boxplot






#### a.Carotene loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for a.Carotene
data_a.Carotene_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "a.Car.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_a.Carotene_plots_cuscutasub$Tissue.code <- factor(data_a.Carotene_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_a.Carotene = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Cuscuta_a.Carotene[[sub]] = p
  
}

#### a.Carotene. Monogynella ####
plot_list_Cuscuta_a.Carotene[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "a.Car.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_a.Carotene <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_a.Carotene_Monogynella <- with(data_Monogynella_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_Monogynella)
boxplot_positions_a.Carotene_Monogynella <- as.data.frame(box.rslt_a.Carotene_Monogynella$stats)

# what are these column tissue codes?
tissues_a.Carotene_Monogynella <- levels(data_Monogynella_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_Monogynella) <- tissues_a.Carotene_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_Monogynella <- boxplot_positions_a.Carotene_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_Monogynella
# now join these positions to cbd
left_join(cbd_a.Carotene_Monogynella, top_positions_a.Carotene_Monogynella, by = "Tissue.code") -> cbd_a.Carotene_Monogynella

# calculate how much to nudge
data_Monogynella_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_Monogynella
cbd_a.Carotene_Monogynella$nudged <- max_a.Carotene_Monogynella$max * 1.05

# add CLDs to plot
a.Carotene_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_a.Carotene_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_Monogynella_boxplot


a.Carotene_Monogynella_boxplot


#### a.Carotene Cuscuta ####
plot_list_Cuscuta_a.Carotene[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "a.Car.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_a.Carotene <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_a.Carotene_Cuscuta <- with(data_Cuscuta_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_Cuscuta)
boxplot_positions_a.Carotene_Cuscuta <- as.data.frame(box.rslt_a.Carotene_Cuscuta$stats)

# what are these column tissue codes?
tissues_a.Carotene_Cuscuta <- levels(data_Cuscuta_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_Cuscuta) <- tissues_a.Carotene_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_Cuscuta <- boxplot_positions_a.Carotene_Cuscuta[5,]


# #add pairwise significance letter groups (compact letter display; CLD)
# cbd_a.Carotene_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__a.Car.Car"]])[["Letters"]])
# colnames(cbd_a.Carotene_Cuscuta)[1] <- "Letter"
# # turn rownames into first column for Tissue.code
# setDT(cbd_a.Carotene_Cuscuta, keep.rownames = "Tissue.code")


# Cuscuta a.Carotene was causing a problem in this pipeline because of the exact matches.. 
wilcox_list[["Cuscuta__a.Car.Car"]]
# QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__a.Car.Car"]])

# create letter groups manually 

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  data_loop_wilcox$FW.norm and data_loop_wilcox$Tissue.code 
# 
#     y     o
#   o 0.066 -
#   h 0.109 -
#   
#   P value adjustment method: BH

#  y and o the same
# y = a, o = a 
# h same as y = a

# make df
cbd_a.Carotene_Cuscuta_Tissue.code <- c("y", "o", "h")
cbd_a.Carotene_Cuscuta_Letters <- c("a", "a", "a")
cbd_a.Carotene_Cuscuta <- data.frame("Tissue.code" = cbd_a.Carotene_Cuscuta_Tissue.code, "Letter" = cbd_a.Carotene_Cuscuta_Letters)


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_Cuscuta
# now join these positions to cbd
left_join(cbd_a.Carotene_Cuscuta, top_positions_a.Carotene_Cuscuta, by = "Tissue.code") -> cbd_a.Carotene_Cuscuta

# calculate how much to nudge
data_Cuscuta_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_Cuscuta
cbd_a.Carotene_Cuscuta$nudged <- (max_a.Carotene_Cuscuta$max + 0.01) * 1.05


# add CLDs to plot
a.Carotene_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_a.Carotene_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_Cuscuta_boxplot


a.Carotene_Cuscuta_boxplot



#### a.Carotene Grammica ####
plot_list_Cuscuta_a.Carotene[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "a.Car.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_a.Carotene <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_a.Carotene_Grammica <- with(data_Grammica_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_Grammica)
boxplot_positions_a.Carotene_Grammica <- as.data.frame(box.rslt_a.Carotene_Grammica$stats)

# what are these column tissue codes?
tissues_a.Carotene_Grammica <- levels(data_Grammica_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_Grammica) <- tissues_a.Carotene_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_Grammica <- boxplot_positions_a.Carotene_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_Grammica
# now join these positions to cbd
left_join(cbd_a.Carotene_Grammica, top_positions_a.Carotene_Grammica, by = "Tissue.code") -> cbd_a.Carotene_Grammica

# calculate how much to nudge
data_Grammica_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_Grammica
cbd_a.Carotene_Grammica$nudged <- (max_a.Carotene_Grammica$max + 0.01) * 1.05

# add CLDs to plot
a.Carotene_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_a.Carotene_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_Grammica_boxplot


a.Carotene_Grammica_boxplot




#### a.Carotene C_purpurata alone ####

data_C_purpurata_a.Carotene <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == "C_purpurata")

a.Carotene_C_purpurata_boxplot <- ggplot(data_C_purpurata_a.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "a.Car.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

a.Carotene_C_purpurata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_a.Carotene <- dplyr::filter(data_a.Carotene_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_a.Carotene_C_purpurata <- with(data_C_purpurata_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_C_purpurata)
boxplot_positions_a.Carotene_C_purpurata <- as.data.frame(box.rslt_a.Carotene_C_purpurata$stats)

# what are these column tissue codes?
tissues_a.Carotene_C_purpurata <- levels(data_C_purpurata_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_C_purpurata) <- tissues_a.Carotene_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_C_purpurata <- boxplot_positions_a.Carotene_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_C_purpurata
# now join these positions to cbd
left_join(cbd_a.Carotene_C_purpurata, top_positions_a.Carotene_C_purpurata, by = "Tissue.code") -> cbd_a.Carotene_C_purpurata

# calculate how much to nudge
data_C_purpurata_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_C_purpurata
cbd_a.Carotene_C_purpurata$nudged <- max_a.Carotene_C_purpurata$max * 1.05

# add CLDs to plot
a.Carotene_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_a.Carotene_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_C_purpurata_boxplot


a.Carotene_C_purpurata_boxplot




#### b.Carotene loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for b.Carotene
data_b.Carotene_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "b.Car.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_b.Carotene_plots_cuscutasub$Tissue.code <- factor(data_b.Carotene_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))

plot_list_Cuscuta_b.Carotene = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_b.Carotene_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x= element_text(angle = -30, hjust = 0, size = 4), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Cuscuta_b.Carotene[[sub]] = p
  
}

#### b.Carotene. Monogynella ####
plot_list_Cuscuta_b.Carotene[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "b.Car.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_Monogynella_boxplot

# KW not sig so no post hoc
b.Carotene_Monogynella_boxplot


#### b.Carotene Cuscuta ####
plot_list_Cuscuta_b.Carotene[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "b.Car.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_Cuscuta_boxplot

# KW not sig so no post hoc

b.Carotene_Cuscuta_boxplot



#### b.Carotene Grammica ####
plot_list_Cuscuta_b.Carotene[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "b.Car.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_b.Carotene <- dplyr::filter(data_b.Carotene_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_b.Carotene_Grammica <- with(data_Grammica_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_Grammica)
boxplot_positions_b.Carotene_Grammica <- as.data.frame(box.rslt_b.Carotene_Grammica$stats)

# what are these column tissue codes?
tissues_b.Carotene_Grammica <- levels(data_Grammica_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_Grammica) <- tissues_b.Carotene_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_Grammica <- boxplot_positions_b.Carotene_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_Grammica
# now join these positions to cbd
left_join(cbd_b.Carotene_Grammica, top_positions_b.Carotene_Grammica, by = "Tissue.code") -> cbd_b.Carotene_Grammica

# calculate how much to nudge
data_Grammica_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_Grammica
cbd_b.Carotene_Grammica$nudged <- max_b.Carotene_Grammica$max * 1.05

# add CLDs to plot
b.Carotene_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_b.Carotene_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_Grammica_boxplot


b.Carotene_Grammica_boxplot




#### b.Carotene C_purpurata alone ####

data_C_purpurata_b.Carotene <- dplyr::filter(data_b.Carotene_plots_cuscutasub, Subgenus == "C_purpurata")

b.Carotene_C_purpurata_boxplot <- ggplot(data_C_purpurata_b.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x= element_text(angle = -30, hjust = 0, size = 4), 
        axis.ticks.x= element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "b.Car.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_b.Carotene <- dplyr::filter(data_b.Carotene_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_b.Carotene_C_purpurata <- with(data_C_purpurata_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_purpurata)
boxplot_positions_b.Carotene_C_purpurata <- as.data.frame(box.rslt_b.Carotene_C_purpurata$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_purpurata <- levels(data_C_purpurata_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_purpurata) <- tissues_b.Carotene_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_purpurata <- boxplot_positions_b.Carotene_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_purpurata
# now join these positions to cbd
left_join(cbd_b.Carotene_C_purpurata, top_positions_b.Carotene_C_purpurata, by = "Tissue.code") -> cbd_b.Carotene_C_purpurata

# calculate how much to nudge
data_C_purpurata_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_purpurata
cbd_b.Carotene_C_purpurata$nudged <- max_b.Carotene_C_purpurata$max * 1.05

# add CLDs to plot
b.Carotene_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_b.Carotene_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_purpurata_boxplot


b.Carotene_C_purpurata_boxplot



#### Tot.Car loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for Tot.Car
data_Tot.Car_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Tot.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Tot.Car_plots_cuscutasub$Tissue.code <- factor(data_Tot.Car_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_Tot.Car = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Tot.Car_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 8))
  plot_list_Cuscuta_Tot.Car[[sub]] = p
  
}

#### Tot.Car. Monogynella ####
plot_list_Cuscuta_Tot.Car[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Tot.Car <- dplyr::filter(data_Tot.Car_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Tot.Car_Monogynella <- with(data_Monogynella_Tot.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_Monogynella)
boxplot_positions_Tot.Car_Monogynella <- as.data.frame(box.rslt_Tot.Car_Monogynella$stats)

# what are these column tissue codes?
tissues_Tot.Car_Monogynella <- levels(data_Monogynella_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_Monogynella) <- tissues_Tot.Car_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_Monogynella <- boxplot_positions_Tot.Car_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_Monogynella
# now join these positions to cbd
left_join(cbd_Tot.Car_Monogynella, top_positions_Tot.Car_Monogynella, by = "Tissue.code") -> cbd_Tot.Car_Monogynella

# calculate how much to nudge
data_Monogynella_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_Monogynella
cbd_Tot.Car_Monogynella$nudged <- (max_Tot.Car_Monogynella$max + 0.05) * 1.05


# add CLDs to plot
Tot.Car_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Car_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_Monogynella_boxplot


Tot.Car_Monogynella_boxplot


#### Tot.Car Cuscuta ####
plot_list_Cuscuta_Tot.Car[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_Cuscuta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Cuscuta_Tot.Car <- dplyr::filter(data_Tot.Car_plots_cuscutasub, Subgenus == "Cuscuta")

box.rslt_Tot.Car_Cuscuta <- with(data_Cuscuta_Tot.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_Cuscuta)
boxplot_positions_Tot.Car_Cuscuta <- as.data.frame(box.rslt_Tot.Car_Cuscuta$stats)

# what are these column tissue codes?
tissues_Tot.Car_Cuscuta <- levels(data_Cuscuta_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_Cuscuta) <- tissues_Tot.Car_Cuscuta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_Cuscuta <- boxplot_positions_Tot.Car_Cuscuta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_Cuscuta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Cuscuta__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_Cuscuta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_Cuscuta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_Cuscuta %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_Cuscuta
# now join these positions to cbd
left_join(cbd_Tot.Car_Cuscuta, top_positions_Tot.Car_Cuscuta, by = "Tissue.code") -> cbd_Tot.Car_Cuscuta

# calculate how much to nudge
data_Cuscuta_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_Cuscuta
cbd_Tot.Car_Cuscuta$nudged <- (max_Tot.Car_Cuscuta$max + 0.05) * 1.05


# add CLDs to plot
Tot.Car_Cuscuta_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Car_Cuscuta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_Cuscuta_boxplot


Tot.Car_Cuscuta_boxplot



#### Tot.Car Grammica ####
plot_list_Cuscuta_Tot.Car[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Tot.Car <- dplyr::filter(data_Tot.Car_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Tot.Car_Grammica <- with(data_Grammica_Tot.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_Grammica)
boxplot_positions_Tot.Car_Grammica <- as.data.frame(box.rslt_Tot.Car_Grammica$stats)

# what are these column tissue codes?
tissues_Tot.Car_Grammica <- levels(data_Grammica_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_Grammica) <- tissues_Tot.Car_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_Grammica <- boxplot_positions_Tot.Car_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_Grammica
# now join these positions to cbd
left_join(cbd_Tot.Car_Grammica, top_positions_Tot.Car_Grammica, by = "Tissue.code") -> cbd_Tot.Car_Grammica

# calculate how much to nudge
data_Grammica_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_Grammica
cbd_Tot.Car_Grammica$nudged <- (max_Tot.Car_Grammica$max + 0.05) * 1.05


# add CLDs to plot
Tot.Car_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Tot.Car_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_Grammica_boxplot


Tot.Car_Grammica_boxplot




#### Tot.Car C. purpurata alone ####

data_C_purpurata_Tot.Car <- dplyr::filter(data_Tot.Car_plots_cuscutasub, Subgenus == "C_purpurata")

Tot.Car_C_purpurata_boxplot <- ggplot(data_C_purpurata_Tot.Car, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 8)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Tot.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

Tot.Car_C_purpurata_boxplot

# not sig







#### NVZ.Car loop through Monogynella, Cuscuta, and Grammica ####
loop_subgenera <- c("Monogynella","Cuscuta", "Grammica")
# dplyr::filter for NVZ.Car
data_NVZ.Car_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "NVZ.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_NVZ.Car_plots_cuscutasub$Tissue.code <- factor(data_NVZ.Car_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Cuscuta_NVZ.Car = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_NVZ.Car_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on Cuscuta plots
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Cuscuta_NVZ.Car[[sub]] = p
  
}

#### NVZ.Car Monogynella ####
plot_list_Cuscuta_NVZ.Car[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "NVZ.Car" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_Monogynella_boxplot

# KW not significant so no post-hoc

NVZ.Car_Monogynella_boxplot



#### NVZ.Car Cuscuta ####
plot_list_Cuscuta_NVZ.Car[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "NVZ.Car" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_Cuscuta_boxplot

# KW not significant so no post-hoc

NVZ.Car_Cuscuta_boxplot



#### NVZ.Car Grammica ####
plot_list_Cuscuta_NVZ.Car[["Grammica"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "NVZ.Car" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_Grammica_boxplot

# KW not significant so no post hoc

NVZ.Car_Grammica_boxplot




#### NVZ.Car C. purpurata alone ####
data_C_purpurata_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_cuscutasub, Subgenus == "C_purpurata")

NVZ.Car_C_purpurata_boxplot <- ggplot(data_C_purpurata_NVZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100))+ 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "NVZ.Car" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE)

NVZ.Car_C_purpurata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_purpurata_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_cuscutasub, Subgenus == "C_purpurata")

box.rslt_NVZ.Car_C_purpurata <- with(data_C_purpurata_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_purpurata)
boxplot_positions_NVZ.Car_C_purpurata <- as.data.frame(box.rslt_NVZ.Car_C_purpurata$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_purpurata <- levels(data_C_purpurata_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_purpurata) <- tissues_NVZ.Car_C_purpurata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_purpurata <- boxplot_positions_NVZ.Car_C_purpurata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_C_purpurata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["C_purpurata__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_purpurata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_purpurata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_purpurata %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_purpurata
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_purpurata, top_positions_NVZ.Car_C_purpurata, by = "Tissue.code") -> cbd_NVZ.Car_C_purpurata

# calculate how much to nudge
data_C_purpurata_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_purpurata
cbd_NVZ.Car_C_purpurata$nudged <- max_NVZ.Car_C_purpurata$max * 1.05

# add CLDs to plot
NVZ.Car_C_purpurata_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_NVZ.Car_C_purpurata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_purpurata_boxplot


NVZ.Car_C_purpurata_boxplot



#### add "absent" to those with mean FW.norm = 0 in loop ####

Subgenus_list_Subgenus <- unique(data_long_calcs$Subgenus)
Pigment_list_Subgenus <- unique(data_long_calcs$Pigment)
label <- "absent"
absent <- data.frame(label)

Subgenus_carotenoid_absent_list <- list()

for (subgenus in Subgenus_list_Subgenus) {
  for (pigment in Pigment_list_Subgenus) {
    data_loop <- data_long_calcs %>% dplyr::filter(Subgenus == subgenus & Pigment == pigment)
    if (nrow(data_loop) > 0) { if (mean(data_loop$FW.norm) == 0) {
      print("absent") -> Subgenus_carotenoid_absent_list[[pigment]][[subgenus]]
    }}}}

Subgenus_carotenoid_absent_list
# $a.Tocopherol
# $a.Tocopherol$Monogynella
# [1] "absent"
# 
# 
# $Neoxanthin
# $Neoxanthin$C_purpurata
# [1] "absent"
# 
# $Neoxanthin$Cuscuta
# [1] "absent"


# add the following geom_text to the above plots
# + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) 




#### COMBINED (faceted) carotenoid plot: subgenus ####
wrap_elements(gridtext::richtext_grob('Total carotenoids', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Tot.Car_ipomoea_boxplot + ggtitle('Ipomoea nil') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) +
  Tot.Car_Monogynella_boxplot + ggtitle('Monogynella') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) +
  Tot.Car_Cuscuta_boxplot + ggtitle('Cuscuta') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Tot.Car_Grammica_boxplot + ggtitle('Grammica') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) +
  Tot.Car_C_purpurata_boxplot + ggtitle('C. purpurata') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  wrap_elements(gridtext::richtext_grob('VAZ', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  VAZ.Car_ipomoea_boxplot +  
  VAZ.Car_Monogynella_boxplot +  
  VAZ.Car_Cuscuta_boxplot + 
  VAZ.Car_Grammica_boxplot + 
  VAZ.Car_C_purpurata_boxplot + 
  # wrap_elements(gridtext::richtext_grob('NVZ.Car', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  # NVZ.Car_ipomoea_boxplot + 
  # NVZ.Car_Monogynella_boxplot +  
  # NVZ.Car_Cuscuta_boxplot + 
  # NVZ.Car_Grammica_boxplot +
  # NVZ.Car_C_purpurata_boxplot +
  wrap_elements(gridtext::richtext_grob('Neoxanthin', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) +
  Neoxanthin_ipomoea_boxplot +
  Neoxanthin_Monogynella_boxplot +
  Neoxanthin_Cuscuta_boxplot + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_Grammica_boxplot +
  Neoxanthin_C_purpurata_boxplot + geom_text(size    = 2, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  wrap_elements(gridtext::richtext_grob('Lutein', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Lutein_ipomoea_boxplot + 
  Lutein_Monogynella_boxplot +  
  Lutein_Cuscuta_boxplot + 
  Lutein_Grammica_boxplot + 
  Lutein_C_purpurata_boxplot + 
  wrap_elements(gridtext::richtext_grob('Lutein epoxide', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Lutein.epoxide_ipomoea_boxplot + 
  Lutein.epoxide_Monogynella_boxplot +  
  Lutein.epoxide_Cuscuta_boxplot + 
  Lutein.epoxide_Grammica_boxplot + 
  Lutein.epoxide_C_purpurata_boxplot + 
  wrap_elements(gridtext::richtext_grob('*a*-Carotene', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  a.Carotene_ipomoea_boxplot +
  a.Carotene_Monogynella_boxplot +  
  a.Carotene_Cuscuta_boxplot + 
  a.Carotene_Grammica_boxplot + 
  a.Carotene_C_purpurata_boxplot + 
  wrap_elements(gridtext::richtext_grob('*b*-Carotene', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  b.Carotene_ipomoea_boxplot +
  b.Carotene_Monogynella_boxplot +  
  b.Carotene_Cuscuta_boxplot + 
  b.Carotene_Grammica_boxplot + 
  b.Carotene_C_purpurata_boxplot + plot_layout(nrow = 7, byrow = T) -> carotenoid_boxplot_new

carotenoid_boxplot_new 

pdf("../output/boxplots/carotenoid_boxplot_new.pdf", width=7,height=8) 
carotenoid_boxplot_new
dev.off()

# while (!is.null(dev.list())) dev.off()
# dev.set(dev.next())














#### Grammica only CHLOROPHYLL PLOTS  ####

#### Chl.a C_australis ####
data_C_australis <-  dplyr::filter(data_long_calcs_Grammica_plot, Species == "C_australis")
data_C_australis_Chl.a <- dplyr::filter(data_C_australis, Pigment == "Chl.a")

# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_C_australis$Tissue.code <- factor(data_C_australis$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))

Chl.a_C_australis_boxplot <- ggplot(data_C_australis_Chl.a, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 6)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

Chl.a_C_australis_boxplot 

# KW not sig so no post hoc
Chl.a_C_australis_boxplot



#### Chl.b C_australis ####
data_C_australis_Chl.b <- dplyr::filter(data_C_australis, Pigment == "Chl.b")


Chl.b_C_australis_boxplot <- ggplot(data_C_australis_Chl.b, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 5))+
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# KW not sig so no post hoc

Chl.b_C_australis_boxplot


#### Tot.Chl C_australis ####
data_C_australis_Tot.Chl <- dplyr::filter(data_C_australis, Pigment == "Tot.Chl")


Tot.Chl_C_australis_boxplot <- ggplot(data_C_australis_Tot.Chl, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 8)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

Tot.Chl_C_australis_boxplot


#### Chl.a.b C_australis plot with x axis ####
data_C_australis_Chl.a.b <- dplyr::filter(data_C_australis, Pigment == "Chl.a.b")


Chl.a.b_C_australis_boxplot <- ggplot(data_C_australis_Chl.a.b, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Mass ratio") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 9)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

Chl.a.b_C_australis_boxplot


#### Chl.a loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Chl.a
data_Chl.a_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Chl.a")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.a_plots_grammicaspe$Tissue.code <- factor(data_Chl.a_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Chl.a = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
    scale_y_continuous(position = "left", limits = c(0, 6))
  plot_list_Grammica_Chl.a[[sub]] = p
  
}


#### Chl.a C_polygonorum ####
plot_list_Grammica_Chl.a[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_polygonorum_boxplot

# KW not sig so no post hoc
Chl.a_C_polygonorum_boxplot


#### Chl.a C_sandwichiana ####
plot_list_Grammica_Chl.a[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_sandwichiana_boxplot


# KW not sig so no post hoc
Chl.a_C_sandwichiana_boxplot


#### Chl.a C_californica ####
plot_list_Grammica_Chl.a[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_californica_boxplot


# KW not sig so no post hoc
Chl.a_C_californica_boxplot

#### Chl.a C_compacta ####
plot_list_Grammica_Chl.a[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_compacta_boxplot


# KW not sig so no post hoc
Chl.a_C_compacta_boxplot


#### Chl.a C_cephalanthii ####
plot_list_Grammica_Chl.a[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Chl.a <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Chl.a_C_cephalanthii <- with(data_C_cephalanthii_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_C_cephalanthii)
boxplot_positions_Chl.a_C_cephalanthii <- as.data.frame(box.rslt_Chl.a_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Chl.a_C_cephalanthii <- levels(data_C_cephalanthii_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_C_cephalanthii) <- tissues_Chl.a_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_C_cephalanthii <- boxplot_positions_Chl.a_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Chl.a_C_cephalanthii, top_positions_Chl.a_C_cephalanthii, by = "Tissue.code") -> cbd_Chl.a_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_C_cephalanthii
cbd_Chl.a_C_cephalanthii$nudged <- max_Chl.a_C_cephalanthii$max * 1.05


# add CLDs to plot
Chl.a_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_C_cephalanthii_boxplot


Chl.a_C_cephalanthii_boxplot


#### Chl.a C_denticulata ####
plot_list_Grammica_Chl.a[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_denticulata_boxplot

# KW not sig so no post hoc 
Chl.a_C_denticulata_boxplot


#### Chl.a C_tasmanica ####
plot_list_Grammica_Chl.a[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Chl.a <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Chl.a_C_tasmanica <- with(data_C_tasmanica_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_C_tasmanica)
boxplot_positions_Chl.a_C_tasmanica <- as.data.frame(box.rslt_Chl.a_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Chl.a_C_tasmanica <- levels(data_C_tasmanica_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_C_tasmanica) <- tissues_Chl.a_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_C_tasmanica <- boxplot_positions_Chl.a_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_C_tasmanica
# now join these positions to cbd
left_join(cbd_Chl.a_C_tasmanica, top_positions_Chl.a_C_tasmanica, by = "Tissue.code") -> cbd_Chl.a_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_C_tasmanica
cbd_Chl.a_C_tasmanica$nudged <- max_Chl.a_C_tasmanica$max * 1.05


# add CLDs to plot
Chl.a_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_C_tasmanica_boxplot


Chl.a_C_tasmanica_boxplot


#### Chl.a C_costaricensis ####
plot_list_Grammica_Chl.a[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Chl.a <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Chl.a_C_costaricensis <- with(data_C_costaricensis_Chl.a, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_C_costaricensis)
boxplot_positions_Chl.a_C_costaricensis <- as.data.frame(box.rslt_Chl.a_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Chl.a_C_costaricensis <- levels(data_C_costaricensis_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_C_costaricensis) <- tissues_Chl.a_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_C_costaricensis <- boxplot_positions_Chl.a_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_C_costaricensis
# now join these positions to cbd
left_join(cbd_Chl.a_C_costaricensis, top_positions_Chl.a_C_costaricensis, by = "Tissue.code") -> cbd_Chl.a_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_C_costaricensis
cbd_Chl.a_C_costaricensis$nudged <- max_Chl.a_C_costaricensis$max * 1.05


# add CLDs to plot
Chl.a_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_C_costaricensis_boxplot


Chl.a_C_costaricensis_boxplot


#### Chl.a C_indecora species alone ####

data_Grammica_Chl.a <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == "C_indecora")

Chl.a_C_indecora_boxplot <- ggplot(data_Grammica_Chl.a, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 6)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Chl.a <- dplyr::filter(data_Chl.a_plots_grammicaspe, Species == "C_indecora")

box.rslt_Chl.a_C_indecora <- with(data_C_indecora_Chl.a, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a_C_indecora)
boxplot_positions_Chl.a_C_indecora <- as.data.frame(box.rslt_Chl.a_C_indecora$stats)

# what are these column tissue codes?
tissues_Chl.a_C_indecora <- levels(data_C_indecora_Chl.a$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a_C_indecora) <- tissues_Chl.a_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a_C_indecora <- boxplot_positions_Chl.a_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Chl.a"]])[["Letters"]])
colnames(cbd_Chl.a_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a_C_indecora
# now join these positions to cbd
left_join(cbd_Chl.a_C_indecora, top_positions_Chl.a_C_indecora, by = "Tissue.code") -> cbd_Chl.a_C_indecora

# calculate how much to nudge
data_C_indecora_Chl.a %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.a_C_indecora
cbd_Chl.a_C_indecora$nudged <- (max_Chl.a_C_indecora$max + 0.00) * 1.05

# add CLDs to plot
Chl.a_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a_C_indecora_boxplot


Chl.a_C_indecora_boxplot



#### Chl.b loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Chl.b
data_Chl.b_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Chl.b")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.b_plots_grammicaspe$Tissue.code <- factor(data_Chl.b_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Chl.b = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 5)) 
  plot_list_Grammica_Chl.b[[sub]] = p
  
}


#### Chl.b C_polygonorum ####
plot_list_Grammica_Chl.b[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_polygonorum_boxplot

# KW not sig so no post hoc
Chl.b_C_polygonorum_boxplot


#### Chl.b C_sandwichiana ####
plot_list_Grammica_Chl.b[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_sandwichiana_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_Chl.b_C_sandwichiana <- with(data_C_sandwichiana_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_sandwichiana)
boxplot_positions_Chl.b_C_sandwichiana <- as.data.frame(box.rslt_Chl.b_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_Chl.b_C_sandwichiana <- levels(data_C_sandwichiana_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_sandwichiana) <- tissues_Chl.b_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_sandwichiana <- boxplot_positions_Chl.b_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_sandwichiana, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_sandwichiana %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_sandwichiana
# now join these positions to cbd
left_join(cbd_Chl.b_C_sandwichiana, top_positions_Chl.b_C_sandwichiana, by = "Tissue.code") -> cbd_Chl.b_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_sandwichiana
cbd_Chl.b_C_sandwichiana$nudged <- max_Chl.b_C_sandwichiana$max * 1.05


# add CLDs to plot
Chl.b_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_sandwichiana_boxplot


Chl.b_C_sandwichiana_boxplot


#### Chl.b C_californica ####
plot_list_Grammica_Chl.b[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_californica")

box.rslt_Chl.b_C_californica <- with(data_C_californica_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_californica)
boxplot_positions_Chl.b_C_californica <- as.data.frame(box.rslt_Chl.b_C_californica$stats)

# what are these column tissue codes?
tissues_Chl.b_C_californica <- levels(data_C_californica_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_californica) <- tissues_Chl.b_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_californica <- boxplot_positions_Chl.b_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_californica
# now join these positions to cbd
left_join(cbd_Chl.b_C_californica, top_positions_Chl.b_C_californica, by = "Tissue.code") -> cbd_Chl.b_C_californica

# calculate how much to nudge
data_C_californica_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_californica
cbd_Chl.b_C_californica$nudged <- max_Chl.b_C_californica$max * 1.05


# add CLDs to plot
Chl.b_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_californica_boxplot


Chl.b_C_californica_boxplot


#### Chl.b C_compacta ####
plot_list_Grammica_Chl.b[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_compacta")

box.rslt_Chl.b_C_compacta <- with(data_C_compacta_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_compacta)
boxplot_positions_Chl.b_C_compacta <- as.data.frame(box.rslt_Chl.b_C_compacta$stats)

# what are these column tissue codes?
tissues_Chl.b_C_compacta <- levels(data_C_compacta_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_compacta) <- tissues_Chl.b_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_compacta <- boxplot_positions_Chl.b_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_compacta
# now join these positions to cbd
left_join(cbd_Chl.b_C_compacta, top_positions_Chl.b_C_compacta, by = "Tissue.code") -> cbd_Chl.b_C_compacta

# calculate how much to nudge
data_C_compacta_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_compacta
cbd_Chl.b_C_compacta$nudged <- max_Chl.b_C_compacta$max * 1.05


# add CLDs to plot
Chl.b_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_compacta_boxplot


Chl.b_C_compacta_boxplot


#### Chl.b C_cephalanthii ####
plot_list_Grammica_Chl.b[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Chl.b_C_cephalanthii <- with(data_C_cephalanthii_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_cephalanthii)
boxplot_positions_Chl.b_C_cephalanthii <- as.data.frame(box.rslt_Chl.b_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Chl.b_C_cephalanthii <- levels(data_C_cephalanthii_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_cephalanthii) <- tissues_Chl.b_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_cephalanthii <- boxplot_positions_Chl.b_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Chl.b_C_cephalanthii, top_positions_Chl.b_C_cephalanthii, by = "Tissue.code") -> cbd_Chl.b_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_cephalanthii
cbd_Chl.b_C_cephalanthii$nudged <- max_Chl.b_C_cephalanthii$max * 1.05


# add CLDs to plot
Chl.b_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_cephalanthii_boxplot


Chl.b_C_cephalanthii_boxplot


#### Chl.b C_denticulata ####
plot_list_Grammica_Chl.b[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_denticulata_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_denticulata")

box.rslt_Chl.b_C_denticulata <- with(data_C_denticulata_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_denticulata)
boxplot_positions_Chl.b_C_denticulata <- as.data.frame(box.rslt_Chl.b_C_denticulata$stats)

# what are these column tissue codes?
tissues_Chl.b_C_denticulata <- levels(data_C_denticulata_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_denticulata) <- tissues_Chl.b_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_denticulata <- boxplot_positions_Chl.b_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_denticulata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_denticulata %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_denticulata
# now join these positions to cbd
left_join(cbd_Chl.b_C_denticulata, top_positions_Chl.b_C_denticulata, by = "Tissue.code") -> cbd_Chl.b_C_denticulata

# calculate how much to nudge
data_C_denticulata_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_denticulata
cbd_Chl.b_C_denticulata$nudged <- max_Chl.b_C_denticulata$max * 1.05


# add CLDs to plot
Chl.b_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_denticulata_boxplot


Chl.b_C_denticulata_boxplot



#### Chl.b C_tasmanica ####
plot_list_Grammica_Chl.b[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Chl.b_C_tasmanica <- with(data_C_tasmanica_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_tasmanica)
boxplot_positions_Chl.b_C_tasmanica <- as.data.frame(box.rslt_Chl.b_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Chl.b_C_tasmanica <- levels(data_C_tasmanica_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_tasmanica) <- tissues_Chl.b_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_tasmanica <- boxplot_positions_Chl.b_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_tasmanica
# now join these positions to cbd
left_join(cbd_Chl.b_C_tasmanica, top_positions_Chl.b_C_tasmanica, by = "Tissue.code") -> cbd_Chl.b_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_tasmanica
cbd_Chl.b_C_tasmanica$nudged <- max_Chl.b_C_tasmanica$max * 1.05


# add CLDs to plot
Chl.b_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_tasmanica_boxplot


Chl.b_C_tasmanica_boxplot


#### Chl.b C_costaricensis ####
plot_list_Grammica_Chl.b[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.b_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Chl.b_C_costaricensis <- with(data_C_costaricensis_Chl.b, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_costaricensis)
boxplot_positions_Chl.b_C_costaricensis <- as.data.frame(box.rslt_Chl.b_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Chl.b_C_costaricensis <- levels(data_C_costaricensis_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_costaricensis) <- tissues_Chl.b_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_costaricensis <- boxplot_positions_Chl.b_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_costaricensis
# now join these positions to cbd
left_join(cbd_Chl.b_C_costaricensis, top_positions_Chl.b_C_costaricensis, by = "Tissue.code") -> cbd_Chl.b_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_costaricensis
cbd_Chl.b_C_costaricensis$nudged <- max_Chl.b_C_costaricensis$max * 1.05


# add CLDs to plot
Chl.b_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_costaricensis_boxplot


Chl.b_C_costaricensis_boxplot


#### Chl.b C_indecora species alone ####

data_Grammica_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_indecora")

Chl.b_C_indecora_boxplot <- ggplot(data_Grammica_Chl.b, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 5)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.b" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Chl.b <- dplyr::filter(data_Chl.b_plots_grammicaspe, Species == "C_indecora")

box.rslt_Chl.b_C_indecora <- with(data_C_indecora_Chl.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.b_C_indecora)
boxplot_positions_Chl.b_C_indecora <- as.data.frame(box.rslt_Chl.b_C_indecora$stats)

# what are these column tissue codes?
tissues_Chl.b_C_indecora <- levels(data_C_indecora_Chl.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.b_C_indecora) <- tissues_Chl.b_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.b_C_indecora <- boxplot_positions_Chl.b_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.b_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Chl.b"]])[["Letters"]])
colnames(cbd_Chl.b_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.b_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.b_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Chl.b_C_indecora
# now join these positions to cbd
left_join(cbd_Chl.b_C_indecora, top_positions_Chl.b_C_indecora, by = "Tissue.code") -> cbd_Chl.b_C_indecora

# calculate how much to nudge
data_C_indecora_Chl.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Chl.b_C_indecora
cbd_Chl.b_C_indecora$nudged <- (max_Chl.b_C_indecora$max * 1.05)

# add CLDs to plot
Chl.b_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.b_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.b_C_indecora_boxplot


Chl.b_C_indecora_boxplot


#### Tot.Chl loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Tot.Chl
data_Tot.Chl_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Tot.Chl")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Tot.Chl_plots_grammicaspe$Tissue.code <- factor(data_Tot.Chl_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Tot.Chl = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 8))
  plot_list_Grammica_Tot.Chl[[sub]] = p
  
}

#### Tot.Chl C_polygonorum ####
plot_list_Grammica_Tot.Chl[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_polygonorum_boxplot

# KW not sig so no post hoc
Tot.Chl_C_polygonorum_boxplot


#### Tot.Chl C_sandwichiana ####
plot_list_Grammica_Tot.Chl[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_sandwichiana_boxplot


# KW not sig so no post hoc
Tot.Chl_C_sandwichiana_boxplot


#### Tot.Chl C_californica ####
plot_list_Grammica_Tot.Chl[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_californica")

box.rslt_Tot.Chl_C_californica <- with(data_C_californica_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_californica)
boxplot_positions_Tot.Chl_C_californica <- as.data.frame(box.rslt_Tot.Chl_C_californica$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_californica <- levels(data_C_californica_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_californica) <- tissues_Tot.Chl_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_californica <- boxplot_positions_Tot.Chl_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_californica
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_californica, top_positions_Tot.Chl_C_californica, by = "Tissue.code") -> cbd_Tot.Chl_C_californica

# calculate how much to nudge
data_C_californica_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_californica
cbd_Tot.Chl_C_californica$nudged <- max_Tot.Chl_C_californica$max * 1.05


# add CLDs to plot
Tot.Chl_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_californica_boxplot


Tot.Chl_C_californica_boxplot


#### Tot.Chl C_compacta ####
plot_list_Grammica_Tot.Chl[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_compacta")

box.rslt_Tot.Chl_C_compacta <- with(data_C_compacta_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_compacta)
boxplot_positions_Tot.Chl_C_compacta <- as.data.frame(box.rslt_Tot.Chl_C_compacta$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_compacta <- levels(data_C_compacta_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_compacta) <- tissues_Tot.Chl_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_compacta <- boxplot_positions_Tot.Chl_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_compacta
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_compacta, top_positions_Tot.Chl_C_compacta, by = "Tissue.code") -> cbd_Tot.Chl_C_compacta

# calculate how much to nudge
data_C_compacta_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_compacta
cbd_Tot.Chl_C_compacta$nudged <- max_Tot.Chl_C_compacta$max * 1.05


# add CLDs to plot
Tot.Chl_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_compacta_boxplot


Tot.Chl_C_compacta_boxplot


#### Tot.Chl C_cephalanthii ####
plot_list_Grammica_Tot.Chl[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Tot.Chl_C_cephalanthii <- with(data_C_cephalanthii_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_cephalanthii)
boxplot_positions_Tot.Chl_C_cephalanthii <- as.data.frame(box.rslt_Tot.Chl_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_cephalanthii <- levels(data_C_cephalanthii_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_cephalanthii) <- tissues_Tot.Chl_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_cephalanthii <- boxplot_positions_Tot.Chl_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_cephalanthii, top_positions_Tot.Chl_C_cephalanthii, by = "Tissue.code") -> cbd_Tot.Chl_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_cephalanthii
cbd_Tot.Chl_C_cephalanthii$nudged <- max_Tot.Chl_C_cephalanthii$max * 1.05


# add CLDs to plot
Tot.Chl_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_cephalanthii_boxplot


Tot.Chl_C_cephalanthii_boxplot


#### Tot.Chl C_denticulata ####
plot_list_Grammica_Tot.Chl[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_denticulata_boxplot


# KW not sig so no post hoc

Tot.Chl_C_denticulata_boxplot



#### Tot.Chl C_tasmanica ####
plot_list_Grammica_Tot.Chl[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Tot.Chl_C_tasmanica <- with(data_C_tasmanica_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_tasmanica)
boxplot_positions_Tot.Chl_C_tasmanica <- as.data.frame(box.rslt_Tot.Chl_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_tasmanica <- levels(data_C_tasmanica_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_tasmanica) <- tissues_Tot.Chl_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_tasmanica <- boxplot_positions_Tot.Chl_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_tasmanica
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_tasmanica, top_positions_Tot.Chl_C_tasmanica, by = "Tissue.code") -> cbd_Tot.Chl_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_tasmanica
cbd_Tot.Chl_C_tasmanica$nudged <- max_Tot.Chl_C_tasmanica$max * 1.05


# add CLDs to plot
Tot.Chl_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_tasmanica_boxplot


Tot.Chl_C_tasmanica_boxplot


#### Tot.Chl C_costaricensis ####
plot_list_Grammica_Tot.Chl[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Chl_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Tot.Chl_C_costaricensis <- with(data_C_costaricensis_Tot.Chl, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_costaricensis)
boxplot_positions_Tot.Chl_C_costaricensis <- as.data.frame(box.rslt_Tot.Chl_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_costaricensis <- levels(data_C_costaricensis_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_costaricensis) <- tissues_Tot.Chl_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_costaricensis <- boxplot_positions_Tot.Chl_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_costaricensis
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_costaricensis, top_positions_Tot.Chl_C_costaricensis, by = "Tissue.code") -> cbd_Tot.Chl_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_costaricensis
cbd_Tot.Chl_C_costaricensis$nudged <- max_Tot.Chl_C_costaricensis$max * 1.05


# add CLDs to plot
Tot.Chl_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_costaricensis_boxplot


Tot.Chl_C_costaricensis_boxplot


#### Tot.Chl C_indecora species alone ####

data_Grammica_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_indecora")

Tot.Chl_C_indecora_boxplot <- ggplot(data_Grammica_Tot.Chl, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg Fresh Weight") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 8)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Chl" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Tot.Chl <- dplyr::filter(data_Tot.Chl_plots_grammicaspe, Species == "C_indecora")

box.rslt_Tot.Chl_C_indecora <- with(data_C_indecora_Tot.Chl, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Chl_C_indecora)
boxplot_positions_Tot.Chl_C_indecora <- as.data.frame(box.rslt_Tot.Chl_C_indecora$stats)

# what are these column tissue codes?
tissues_Tot.Chl_C_indecora <- levels(data_C_indecora_Tot.Chl$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Chl_C_indecora) <- tissues_Tot.Chl_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Chl_C_indecora <- boxplot_positions_Tot.Chl_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Chl_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Tot.Chl"]])[["Letters"]])
colnames(cbd_Tot.Chl_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Chl_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Chl_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Chl_C_indecora
# now join these positions to cbd
left_join(cbd_Tot.Chl_C_indecora, top_positions_Tot.Chl_C_indecora, by = "Tissue.code") -> cbd_Tot.Chl_C_indecora

# calculate how much to nudge
data_C_indecora_Tot.Chl %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Chl_C_indecora
cbd_Tot.Chl_C_indecora$nudged <- (max_Tot.Chl_C_indecora$max * 1.05)

# add CLDs to plot
Tot.Chl_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Chl_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Chl_C_indecora_boxplot


Tot.Chl_C_indecora_boxplot






#### Chl.a.b loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Chl.a.b
data_Chl.a.b_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Chl.a.b")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Chl.a.b_plots_grammicaspe$Tissue.code <- factor(data_Chl.a.b_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Chl.a.b = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 9))
  plot_list_Grammica_Chl.a.b[[sub]] = p
  
}

#### Chl.a.b C_polygonorum ####
plot_list_Grammica_Chl.a.b[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_polygonorum_boxplot

# KW not sig so no post hoc
Chl.a.b_C_polygonorum_boxplot


#### Chl.a.b C_sandwichiana ####
plot_list_Grammica_Chl.a.b[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_sandwichiana_boxplot


# KW not sig so no post hoc

Chl.a.b_C_sandwichiana_boxplot


#### Chl.a.b C_californica ####
plot_list_Grammica_Chl.a.b[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_californica_boxplot


# KW not sig so no post hoc
Chl.a.b_C_californica_boxplot


#### Chl.a.b C_compacta ####
plot_list_Grammica_Chl.a.b[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_compacta_boxplot



# KW not sig so no post hoc

Chl.a.b_C_compacta_boxplot


#### Chl.a.b C_cephalanthii ####
plot_list_Grammica_Chl.a.b[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Chl.a.b_C_cephalanthii <- with(data_C_cephalanthii_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_C_cephalanthii)
boxplot_positions_Chl.a.b_C_cephalanthii <- as.data.frame(box.rslt_Chl.a.b_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Chl.a.b_C_cephalanthii <- levels(data_C_cephalanthii_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_C_cephalanthii) <- tissues_Chl.a.b_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_C_cephalanthii <- boxplot_positions_Chl.a.b_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Chl.a.b_C_cephalanthii, top_positions_Chl.a.b_C_cephalanthii, by = "Tissue.code") -> cbd_Chl.a.b_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_C_cephalanthii
cbd_Chl.a.b_C_cephalanthii$nudged <- max_Chl.a.b_C_cephalanthii$max * 1.05


# add CLDs to plot
Chl.a.b_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a.b_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_C_cephalanthii_boxplot


Chl.a.b_C_cephalanthii_boxplot


#### Chl.a.b C_denticulata ####
plot_list_Grammica_Chl.a.b[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_denticulata_boxplot


# KW not sig so no post hoc

Chl.a.b_C_denticulata_boxplot



#### Chl.a.b C_tasmanica ####
plot_list_Grammica_Chl.a.b[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Chl.a.b_C_tasmanica <- with(data_C_tasmanica_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_C_tasmanica)
boxplot_positions_Chl.a.b_C_tasmanica <- as.data.frame(box.rslt_Chl.a.b_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Chl.a.b_C_tasmanica <- levels(data_C_tasmanica_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_C_tasmanica) <- tissues_Chl.a.b_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_C_tasmanica <- boxplot_positions_Chl.a.b_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_C_tasmanica
# now join these positions to cbd
left_join(cbd_Chl.a.b_C_tasmanica, top_positions_Chl.a.b_C_tasmanica, by = "Tissue.code") -> cbd_Chl.a.b_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_C_tasmanica
cbd_Chl.a.b_C_tasmanica$nudged <- max_Chl.a.b_C_tasmanica$max * 1.05


# add CLDs to plot
Chl.a.b_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a.b_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_C_tasmanica_boxplot


Chl.a.b_C_tasmanica_boxplot


#### Chl.a.b C_costaricensis ####
plot_list_Grammica_Chl.a.b[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Chl.a.b_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Chl.a.b_C_costaricensis <- with(data_C_costaricensis_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_C_costaricensis)
boxplot_positions_Chl.a.b_C_costaricensis <- as.data.frame(box.rslt_Chl.a.b_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Chl.a.b_C_costaricensis <- levels(data_C_costaricensis_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_C_costaricensis) <- tissues_Chl.a.b_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_C_costaricensis <- boxplot_positions_Chl.a.b_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_C_costaricensis
# now join these positions to cbd
left_join(cbd_Chl.a.b_C_costaricensis, top_positions_Chl.a.b_C_costaricensis, by = "Tissue.code") -> cbd_Chl.a.b_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_C_costaricensis
cbd_Chl.a.b_C_costaricensis$nudged <- max_Chl.a.b_C_costaricensis$max * 1.05


# add CLDs to plot
Chl.a.b_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a.b_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_C_costaricensis_boxplot


Chl.a.b_C_costaricensis_boxplot


#### Chl.a.b C_indecora species alone ####

data_Grammica_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == "C_indecora")

Chl.a.b_C_indecora_boxplot <- ggplot(data_Grammica_Chl.a.b, aes(x=Tissue.code, y=FW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Mass ratio") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 9)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Chl.a.b" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Chl.a.b <- dplyr::filter(data_Chl.a.b_plots_grammicaspe, Species == "C_indecora")

box.rslt_Chl.a.b_C_indecora <- with(data_C_indecora_Chl.a.b, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Chl.a.b_C_indecora)
boxplot_positions_Chl.a.b_C_indecora <- as.data.frame(box.rslt_Chl.a.b_C_indecora$stats)

# what are these column tissue codes?
tissues_Chl.a.b_C_indecora <- levels(data_C_indecora_Chl.a.b$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Chl.a.b_C_indecora) <- tissues_Chl.a.b_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Chl.a.b_C_indecora <- boxplot_positions_Chl.a.b_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Chl.a.b_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Chl.a.b"]])[["Letters"]])
colnames(cbd_Chl.a.b_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Chl.a.b_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Chl.a.b_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Chl.a.b_C_indecora
# now join these positions to cbd
left_join(cbd_Chl.a.b_C_indecora, top_positions_Chl.a.b_C_indecora, by = "Tissue.code") -> cbd_Chl.a.b_C_indecora

# calculate how much to nudge
data_C_indecora_Chl.a.b %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm)) -> max_Chl.a.b_C_indecora
cbd_Chl.a.b_C_indecora$nudged <- (max_Chl.a.b_C_indecora$max * 1.05)

# add CLDs to plot
Chl.a.b_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Chl.a.b_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Chl.a.b_C_indecora_boxplot


Chl.a.b_C_indecora_boxplot



#### add "absent" to those with mean FW.norm = 0 in loop ####

Species_list_Grammica <- unique(data_long_calcs_Grammica_plot$Species)
Pigment_list_Grammica <- unique(data_long_calcs_Grammica_plot$Pigment)
label <- "absent"
absent <- data.frame(label)

Grammica_carotenoid_absent_list <- list()

for (species in Species_list_Grammica) {
  for (pigment in Pigment_list_Grammica) {
    data_loop <- data_long_calcs_Grammica_plot %>% dplyr::filter(Species == species & Pigment == pigment)
    if (mean(data_loop$FW.norm) == 0) {
      print("absent") -> Grammica_carotenoid_absent_list[[pigment]][[species]]
    }}}

# 
# $Neoxanthin
# $Neoxanthin$C_sandwichiana
# [1] "absent"
# 
# $Neoxanthin$C_compacta
# [1] "absent"
# 
# $Neoxanthin$C_costaricensis
# [1] "absent"
# 
# $Neoxanthin$C_tasmanica
# [1] "absent"
# 
# $Neoxanthin$C_californica
# [1] "absent"
# 
# $Neoxanthin$C_indecora
# [1] "absent"
# 
# $Neoxanthin$C_denticulata
# [1] "absent"

# add the following geom_text to the above plots
# + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) 





#### COMBINED (faceted) chlorophyll plot: Grammica ONLY ####
wrap_elements(gridtext::richtext_grob('Chlorophyll *a*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.a_C_australis_boxplot + ggtitle('C. australis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_polygonorum_boxplot + ggtitle('C. polygonorum') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_sandwichiana_boxplot + ggtitle('C. sandwichiana') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_californica_boxplot + ggtitle('C. californica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_compacta_boxplot + ggtitle('C. compacta') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_cephalanthii_boxplot + ggtitle('C. cephalanthii') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_denticulata_boxplot + ggtitle('C. denticulata') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Chl.a_C_tasmanica_boxplot + ggtitle('C. tasmanica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Chl.a_C_costaricensis_boxplot + ggtitle('C. costaricensis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Chl.a_C_indecora_boxplot + ggtitle('C. indecora') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  wrap_elements(gridtext::richtext_grob('Chlorophyll *b*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.b_C_australis_boxplot +
  Chl.b_C_polygonorum_boxplot + 
  Chl.b_C_sandwichiana_boxplot +
  Chl.b_C_californica_boxplot + 
  Chl.b_C_compacta_boxplot +
  Chl.b_C_cephalanthii_boxplot +
  Chl.b_C_denticulata_boxplot +
  Chl.b_C_tasmanica_boxplot +
  Chl.b_C_costaricensis_boxplot +
  Chl.b_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('Total chlorophyll', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Tot.Chl_C_australis_boxplot +
  Tot.Chl_C_polygonorum_boxplot + 
  Tot.Chl_C_sandwichiana_boxplot +
  Tot.Chl_C_californica_boxplot + 
  Tot.Chl_C_compacta_boxplot +
  Tot.Chl_C_cephalanthii_boxplot +
  Tot.Chl_C_denticulata_boxplot +
  Tot.Chl_C_tasmanica_boxplot +
  Tot.Chl_C_costaricensis_boxplot +
  Tot.Chl_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('Chlorophyll *a*:*b*', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Chl.a.b_C_australis_boxplot +
  Chl.a.b_C_polygonorum_boxplot + 
  Chl.a.b_C_sandwichiana_boxplot +
  Chl.a.b_C_californica_boxplot + 
  Chl.a.b_C_compacta_boxplot +
  Chl.a.b_C_cephalanthii_boxplot +
  Chl.a.b_C_denticulata_boxplot +
  Chl.a.b_C_tasmanica_boxplot +
  Chl.a.b_C_costaricensis_boxplot +
  Chl.a.b_C_indecora_boxplot +
  plot_layout(nrow = 4, byrow = T) -> chlorophyll_boxplot_Grammica_new

chlorophyll_boxplot_Grammica_new 

pdf("../output/boxplots/chlorophyll_boxplot_Grammica_new.pdf", width=9,height=5.5) 
chlorophyll_boxplot_Grammica_new
dev.off()


#### ####



#### Grammica only CAROTENOID PLOTS  ####

#### VAZ.Car C_australis ####
data_C_australis <-  dplyr::filter(data_long_calcs_Grammica_plot, Species == "C_australis")
data_C_australis_VAZ.Car <- dplyr::filter(data_C_australis, Pigment == "VAZ.Car")

# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_C_australis$Tissue.code <- factor(data_C_australis$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))

VAZ.Car_C_australis_boxplot <- ggplot(data_C_australis_VAZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

VAZ.Car_C_australis_boxplot 

# KW not sig so no post hoc
VAZ.Car_C_australis_boxplot



#### Neoxanthin C_australis ####
data_C_australis_Neoxanthin <- dplyr::filter(data_C_australis, Pigment == "Neo.Car")


Neoxanthin_C_australis_boxplot <- ggplot(data_C_australis_Neoxanthin, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100))+
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Neo.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# KW not sig so no post hoc

Neoxanthin_C_australis_boxplot


#### Lutein.epoxide C_australis ####
data_C_australis_Lutein.epoxide <- dplyr::filter(data_C_australis, Pigment == "Lut.epo.Car")


Lutein.epoxide_C_australis_boxplot <- ggplot(data_C_australis_Lutein.epoxide, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

Lutein.epoxide_C_australis_boxplot



#### Lutein C_australis ####
data_C_australis_Lutein <- dplyr::filter(data_C_australis, Pigment == "Lut.Car")


Lutein_C_australis_boxplot <- ggplot(data_C_australis_Lutein, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

Lutein_C_australis_boxplot




#### a.Carotene C_australis ####
data_C_australis_a.Carotene <- dplyr::filter(data_C_australis, Pigment == "a.Car.Car")


a.Carotene_C_australis_boxplot <- ggplot(data_C_australis_a.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

a.Carotene_C_australis_boxplot




#### b.Carotene C_australis ####
data_C_australis_b.Carotene <- dplyr::filter(data_C_australis, Pigment == "b.Car.Car")


b.Carotene_C_australis_boxplot <- ggplot(data_C_australis_b.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

b.Carotene_C_australis_boxplot




#### Tot.Car C_australis ####
data_C_australis_Tot.Car <- dplyr::filter(data_C_australis, Pigment == "Tot.Car")


Tot.Car_C_australis_boxplot <- ggplot(data_C_australis_Tot.Car, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg FW") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 6)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc

Tot.Car_C_australis_boxplot


#### NVZ.Car C_australis plot with x axis ####
data_C_australis_NVZ.Car <- dplyr::filter(data_C_australis, Pigment == "NVZ.Car")


NVZ.Car_C_australis_boxplot <- ggplot(data_C_australis_NVZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 100)) +
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# KW not sig so no post hoc
NVZ.Car_C_australis_boxplot


#### VAZ.Car loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for VAZ.Car
data_VAZ.Car_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "VAZ.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_VAZ.Car_plots_grammicaspe$Tissue.code <- factor(data_VAZ.Car_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_VAZ.Car = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
    scale_y_continuous(position = "left", limits = c(0, 100))
  plot_list_Grammica_VAZ.Car[[sub]] = p
  
}


#### VAZ.Car C_polygonorum ####
plot_list_Grammica_VAZ.Car[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_polygonorum_boxplot

# not sig


#### VAZ.Car C_sandwichiana ####
plot_list_Grammica_VAZ.Car[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_sandwichiana_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_VAZ.Car_C_sandwichiana <- with(data_C_sandwichiana_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_sandwichiana)
boxplot_positions_VAZ.Car_C_sandwichiana <- as.data.frame(box.rslt_VAZ.Car_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_sandwichiana <- levels(data_C_sandwichiana_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_sandwichiana) <- tissues_VAZ.Car_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_sandwichiana <- boxplot_positions_VAZ.Car_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_sandwichiana, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_sandwichiana %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_sandwichiana
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_sandwichiana, top_positions_VAZ.Car_C_sandwichiana, by = "Tissue.code") -> cbd_VAZ.Car_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_sandwichiana
cbd_VAZ.Car_C_sandwichiana$nudged <- max_VAZ.Car_C_sandwichiana$max * 1.05


# add CLDs to plot
VAZ.Car_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_sandwichiana_boxplot


VAZ.Car_C_sandwichiana_boxplot


#### VAZ.Car C_californica ####
plot_list_Grammica_VAZ.Car[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_californica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_californica_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_californica")

box.rslt_VAZ.Car_C_californica <- with(data_C_californica_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_californica)
boxplot_positions_VAZ.Car_C_californica <- as.data.frame(box.rslt_VAZ.Car_C_californica$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_californica <- levels(data_C_californica_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_californica) <- tissues_VAZ.Car_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_californica <- boxplot_positions_VAZ.Car_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_californica
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_californica, top_positions_VAZ.Car_C_californica, by = "Tissue.code") -> cbd_VAZ.Car_C_californica

# calculate how much to nudge
data_C_californica_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_californica
cbd_VAZ.Car_C_californica$nudged <- max_VAZ.Car_C_californica$max * 1.05


# add CLDs to plot
VAZ.Car_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_californica_boxplot


VAZ.Car_C_californica_boxplot


#### VAZ.Car C_compacta ####
plot_list_Grammica_VAZ.Car[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_compacta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_compacta_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_compacta")

box.rslt_VAZ.Car_C_compacta <- with(data_C_compacta_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_compacta)
boxplot_positions_VAZ.Car_C_compacta <- as.data.frame(box.rslt_VAZ.Car_C_compacta$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_compacta <- levels(data_C_compacta_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_compacta) <- tissues_VAZ.Car_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_compacta <- boxplot_positions_VAZ.Car_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_compacta
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_compacta, top_positions_VAZ.Car_C_compacta, by = "Tissue.code") -> cbd_VAZ.Car_C_compacta

# calculate how much to nudge
data_C_compacta_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_compacta
cbd_VAZ.Car_C_compacta$nudged <- max_VAZ.Car_C_compacta$max * 1.05


# add CLDs to plot
VAZ.Car_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_compacta_boxplot


VAZ.Car_C_compacta_boxplot


#### VAZ.Car C_cephalanthii ####
plot_list_Grammica_VAZ.Car[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_VAZ.Car_C_cephalanthii <- with(data_C_cephalanthii_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_cephalanthii)
boxplot_positions_VAZ.Car_C_cephalanthii <- as.data.frame(box.rslt_VAZ.Car_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_cephalanthii <- levels(data_C_cephalanthii_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_cephalanthii) <- tissues_VAZ.Car_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_cephalanthii <- boxplot_positions_VAZ.Car_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_cephalanthii
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_cephalanthii, top_positions_VAZ.Car_C_cephalanthii, by = "Tissue.code") -> cbd_VAZ.Car_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_cephalanthii
cbd_VAZ.Car_C_cephalanthii$nudged <- max_VAZ.Car_C_cephalanthii$max * 1.05


# add CLDs to plot
VAZ.Car_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_cephalanthii_boxplot


VAZ.Car_C_cephalanthii_boxplot


#### VAZ.Car C_denticulata ####
plot_list_Grammica_VAZ.Car[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_denticulata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_denticulata")

box.rslt_VAZ.Car_C_denticulata <- with(data_C_denticulata_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_denticulata)
boxplot_positions_VAZ.Car_C_denticulata <- as.data.frame(box.rslt_VAZ.Car_C_denticulata$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_denticulata <- levels(data_C_denticulata_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_denticulata) <- tissues_VAZ.Car_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_denticulata <- boxplot_positions_VAZ.Car_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_denticulata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_denticulata %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_denticulata
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_denticulata, top_positions_VAZ.Car_C_denticulata, by = "Tissue.code") -> cbd_VAZ.Car_C_denticulata

# calculate how much to nudge
data_C_denticulata_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_denticulata
cbd_VAZ.Car_C_denticulata$nudged <- max_VAZ.Car_C_denticulata$max * 1.05


# add CLDs to plot
VAZ.Car_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_denticulata_boxplot


VAZ.Car_C_denticulata_boxplot


#### VAZ.Car C_tasmanica ####
plot_list_Grammica_VAZ.Car[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_VAZ.Car_C_tasmanica <- with(data_C_tasmanica_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_tasmanica)
boxplot_positions_VAZ.Car_C_tasmanica <- as.data.frame(box.rslt_VAZ.Car_C_tasmanica$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_tasmanica <- levels(data_C_tasmanica_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_tasmanica) <- tissues_VAZ.Car_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_tasmanica <- boxplot_positions_VAZ.Car_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_tasmanica
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_tasmanica, top_positions_VAZ.Car_C_tasmanica, by = "Tissue.code") -> cbd_VAZ.Car_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_tasmanica
cbd_VAZ.Car_C_tasmanica$nudged <- max_VAZ.Car_C_tasmanica$max * 1.05


# add CLDs to plot
VAZ.Car_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_tasmanica_boxplot


VAZ.Car_C_tasmanica_boxplot


#### VAZ.Car C_costaricensis ####
plot_list_Grammica_VAZ.Car[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> VAZ.Car_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_VAZ.Car_C_costaricensis <- with(data_C_costaricensis_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_costaricensis)
boxplot_positions_VAZ.Car_C_costaricensis <- as.data.frame(box.rslt_VAZ.Car_C_costaricensis$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_costaricensis <- levels(data_C_costaricensis_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_costaricensis) <- tissues_VAZ.Car_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_costaricensis <- boxplot_positions_VAZ.Car_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_costaricensis
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_costaricensis, top_positions_VAZ.Car_C_costaricensis, by = "Tissue.code") -> cbd_VAZ.Car_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_costaricensis
cbd_VAZ.Car_C_costaricensis$nudged <- max_VAZ.Car_C_costaricensis$max * 1.05


# add CLDs to plot
VAZ.Car_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_costaricensis_boxplot


VAZ.Car_C_costaricensis_boxplot


#### VAZ.Car C_indecora species alone ####

data_Grammica_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_indecora")

VAZ.Car_C_indecora_boxplot <- ggplot(data_Grammica_VAZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "VAZ.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_VAZ.Car <- dplyr::filter(data_VAZ.Car_plots_grammicaspe, Species == "C_indecora")

box.rslt_VAZ.Car_C_indecora <- with(data_C_indecora_VAZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_VAZ.Car_C_indecora)
boxplot_positions_VAZ.Car_C_indecora <- as.data.frame(box.rslt_VAZ.Car_C_indecora$stats)

# what are these column tissue codes?
tissues_VAZ.Car_C_indecora <- levels(data_C_indecora_VAZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_VAZ.Car_C_indecora) <- tissues_VAZ.Car_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_VAZ.Car_C_indecora <- boxplot_positions_VAZ.Car_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_VAZ.Car_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__VAZ.Car"]])[["Letters"]])
colnames(cbd_VAZ.Car_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_VAZ.Car_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_VAZ.Car_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_VAZ.Car_C_indecora
# now join these positions to cbd
left_join(cbd_VAZ.Car_C_indecora, top_positions_VAZ.Car_C_indecora, by = "Tissue.code") -> cbd_VAZ.Car_C_indecora

# calculate how much to nudge
data_C_indecora_VAZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_VAZ.Car_C_indecora
cbd_VAZ.Car_C_indecora$nudged <- max_VAZ.Car_C_indecora$max * 1.05


# add CLDs to plot
VAZ.Car_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_VAZ.Car_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> VAZ.Car_C_indecora_boxplot


VAZ.Car_C_indecora_boxplot



#### Neoxanthin loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Neoxanthin
data_Neoxanthin_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Neo.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Neoxanthin_plots_grammicaspe$Tissue.code <- factor(data_Neoxanthin_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Neoxanthin = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Neoxanthin_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100)) 
  plot_list_Grammica_Neoxanthin[[sub]] = p
  
}


#### Neoxanthin C_polygonorum ####
plot_list_Grammica_Neoxanthin[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Neo.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Neoxanthin_C_polygonorum_boxplot

# KW not sig so no post hoc
Neoxanthin_C_polygonorum_boxplot


#### Neoxanthin C_sandwichiana ####
plot_list_Grammica_Neoxanthin[["C_sandwichiana"]] -> Neoxanthin_C_sandwichiana_boxplot

# no KW and no post hoc
Neoxanthin_C_sandwichiana_boxplot


#### Neoxanthin C_californica ####
plot_list_Grammica_Neoxanthin[["C_californica"]] -> Neoxanthin_C_californica_boxplot

# no KW and no post hoc
Neoxanthin_C_californica_boxplot


#### Neoxanthin C_compacta ####
plot_list_Grammica_Neoxanthin[["C_compacta"]] -> Neoxanthin_C_compacta_boxplot

# no KW and no post hoc
Neoxanthin_C_compacta_boxplot


#### Neoxanthin C_cephalanthii ####
plot_list_Grammica_Neoxanthin[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Neo.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Neoxanthin_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Neoxanthin <- dplyr::filter(data_Neoxanthin_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Neoxanthin_C_cephalanthii <- with(data_C_cephalanthii_Neoxanthin, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Neoxanthin_C_cephalanthii)
boxplot_positions_Neoxanthin_C_cephalanthii <- as.data.frame(box.rslt_Neoxanthin_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Neoxanthin_C_cephalanthii <- levels(data_C_cephalanthii_Neoxanthin$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Neoxanthin_C_cephalanthii) <- tissues_Neoxanthin_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Neoxanthin_C_cephalanthii <- boxplot_positions_Neoxanthin_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Neoxanthin_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Neo.Car"]])[["Letters"]])
colnames(cbd_Neoxanthin_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Neoxanthin_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Neoxanthin_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Neoxanthin_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Neoxanthin_C_cephalanthii, top_positions_Neoxanthin_C_cephalanthii, by = "Tissue.code") -> cbd_Neoxanthin_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Neoxanthin %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Neoxanthin_C_cephalanthii
cbd_Neoxanthin_C_cephalanthii$nudged <- max_Neoxanthin_C_cephalanthii$max * 1.05


# add CLDs to plot
Neoxanthin_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Neoxanthin_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Neoxanthin_C_cephalanthii_boxplot


Neoxanthin_C_cephalanthii_boxplot


#### Neoxanthin C_denticulata ####
plot_list_Grammica_Neoxanthin[["C_denticulata"]] -> Neoxanthin_C_denticulata_boxplot

# no KW and no post hoc
Neoxanthin_C_denticulata_boxplot



#### Neoxanthin C_tasmanica ####
plot_list_Grammica_Neoxanthin[["C_tasmanica"]]  -> Neoxanthin_C_tasmanica_boxplot

# no KW and no post hoc

Neoxanthin_C_tasmanica_boxplot


#### Neoxanthin C_costaricensis ####
plot_list_Grammica_Neoxanthin[["C_costaricensis"]] -> Neoxanthin_C_costaricensis_boxplot

# no KW and no post hoc
Neoxanthin_C_costaricensis_boxplot


#### Neoxanthin C_indecora species alone ####

data_Grammica_Neoxanthin <- dplyr::filter(data_Neoxanthin_plots_grammicaspe, Species == "C_indecora")

Neoxanthin_C_indecora_boxplot <- ggplot(data_Grammica_Neoxanthin, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) 


# no KW and no post hoc
Neoxanthin_C_indecora_boxplot


#### Lutein.epoxide loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Lutein.epoxide
data_Lutein.epoxide_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Lut.epo.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Lutein.epoxide_plots_grammicaspe$Tissue.code <- factor(data_Lutein.epoxide_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Lutein.epoxide = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Lutein.epoxide_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Grammica_Lutein.epoxide[[sub]] = p
  
}

#### Lutein.epoxide C_polygonorum ####
plot_list_Grammica_Lutein.epoxide[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_polygonorum_boxplot

# KW not sig so no post hoc
Lutein.epoxide_C_polygonorum_boxplot


#### Lutein.epoxide C_sandwichiana ####
plot_list_Grammica_Lutein.epoxide[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_sandwichiana_boxplot


# KW not sig so no post hoc
Lutein.epoxide_C_sandwichiana_boxplot


#### Lutein.epoxide C_californica ####
plot_list_Grammica_Lutein.epoxide[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_californica_boxplot


# KW not sig no post hoc
Lutein.epoxide_C_californica_boxplot


#### Lutein.epoxide C_compacta ####
plot_list_Grammica_Lutein.epoxide[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_compacta_boxplot

# KW not sig no post hoc
Lutein.epoxide_C_compacta_boxplot


#### Lutein.epoxide C_cephalanthii ####
plot_list_Grammica_Lutein.epoxide[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Lutein.epoxide_C_cephalanthii <- with(data_C_cephalanthii_Lutein.epoxide, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein.epoxide_C_cephalanthii)
boxplot_positions_Lutein.epoxide_C_cephalanthii <- as.data.frame(box.rslt_Lutein.epoxide_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Lutein.epoxide_C_cephalanthii <- levels(data_C_cephalanthii_Lutein.epoxide$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein.epoxide_C_cephalanthii) <- tissues_Lutein.epoxide_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein.epoxide_C_cephalanthii <- boxplot_positions_Lutein.epoxide_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein.epoxide_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Lut.epo.Car"]])[["Letters"]])
colnames(cbd_Lutein.epoxide_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein.epoxide_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein.epoxide_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Lutein.epoxide_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Lutein.epoxide_C_cephalanthii, top_positions_Lutein.epoxide_C_cephalanthii, by = "Tissue.code") -> cbd_Lutein.epoxide_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Lutein.epoxide %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein.epoxide_C_cephalanthii
cbd_Lutein.epoxide_C_cephalanthii$nudged <- max_Lutein.epoxide_C_cephalanthii$max * 1.05


# add CLDs to plot
Lutein.epoxide_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein.epoxide_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein.epoxide_C_cephalanthii_boxplot


Lutein.epoxide_C_cephalanthii_boxplot


#### Lutein.epoxide C_denticulata ####
plot_list_Grammica_Lutein.epoxide[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_denticulata_boxplot

# not sig



#### Lutein.epoxide C_tasmanica ####
plot_list_Grammica_Lutein.epoxide[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Lutein.epoxide_C_tasmanica <- with(data_C_tasmanica_Lutein.epoxide, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein.epoxide_C_tasmanica)
boxplot_positions_Lutein.epoxide_C_tasmanica <- as.data.frame(box.rslt_Lutein.epoxide_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Lutein.epoxide_C_tasmanica <- levels(data_C_tasmanica_Lutein.epoxide$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein.epoxide_C_tasmanica) <- tissues_Lutein.epoxide_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein.epoxide_C_tasmanica <- boxplot_positions_Lutein.epoxide_C_tasmanica[5,]


# #add pairwise significance letter groups (compact letter display; CLD)
# cbd_Lutein.epoxide_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Lut.epo.Car"]])[["Letters"]])
# colnames(cbd_Lutein.epoxide_C_tasmanica)[1] <- "Letter"
# # turn rownames into first column for Tissue.code
# setDT(cbd_Lutein.epoxide_C_tasmanica, keep.rownames = "Tissue.code")

wilcox_list_Grammica[["C_tasmanica__Lutein.epoxide"]]

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  data_loop_wilcox$FW.norm and data_loop_wilcox$Tissue.code 
# 
#     y      o      h     
#   o -      -      -     
#   h 0.0028 0.0010 -     
#   f 0.0267 0.0125 0.0497
# 
# P value adjustment method: BH

# make df
cbd_Lutein.epoxide_C_tasmanica_Tissue.code <- c("y", "o", "h", "f")
cbd_Lutein.epoxide_C_tasmanica_Letters <- c("a", "b", "c", "c")
cbd_Lutein.epoxide_C_tasmanica <- data.frame("Tissue.code" = cbd_Lutein.epoxide_C_tasmanica_Tissue.code, "Letter" = cbd_Lutein.epoxide_C_tasmanica_Letters)





# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein.epoxide_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Lutein.epoxide_C_tasmanica
# now join these positions to cbd
left_join(cbd_Lutein.epoxide_C_tasmanica, top_positions_Lutein.epoxide_C_tasmanica, by = "Tissue.code") -> cbd_Lutein.epoxide_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Lutein.epoxide %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein.epoxide_C_tasmanica
cbd_Lutein.epoxide_C_tasmanica$nudged <- max_Lutein.epoxide_C_tasmanica$max * 1.05


# add CLDs to plot
Lutein.epoxide_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein.epoxide_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein.epoxide_C_tasmanica_boxplot


Lutein.epoxide_C_tasmanica_boxplot


#### Lutein.epoxide C_costaricensis ####
plot_list_Grammica_Lutein.epoxide[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein.epoxide_C_costaricensis_boxplot

# KW not sig no post hoc
Lutein.epoxide_C_costaricensis_boxplot


#### Lutein.epoxide C_indecora species alone ####

data_Grammica_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_grammicaspe, Species == "C_indecora")

Lutein.epoxide_C_indecora_boxplot <- ggplot(data_Grammica_Lutein.epoxide, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.epo.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Lutein.epoxide <- dplyr::filter(data_Lutein.epoxide_plots_grammicaspe, Species == "C_indecora")

box.rslt_Lutein.epoxide_C_indecora <- with(data_C_indecora_Lutein.epoxide, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein.epoxide_C_indecora)
boxplot_positions_Lutein.epoxide_C_indecora <- as.data.frame(box.rslt_Lutein.epoxide_C_indecora$stats)

# what are these column tissue codes?
tissues_Lutein.epoxide_C_indecora <- levels(data_C_indecora_Lutein.epoxide$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein.epoxide_C_indecora) <- tissues_Lutein.epoxide_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein.epoxide_C_indecora <- boxplot_positions_Lutein.epoxide_C_indecora[5,]


# #add pairwise significance letter groups (compact letter display; CLD)
# cbd_Lutein.epoxide_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Lut.epo.Car"]])[["Letters"]])
# colnames(cbd_Lutein.epoxide_C_indecora)[1] <- "Letter"
# # turn rownames into first column for Tissue.code
# setDT(cbd_Lutein.epoxide_C_indecora, keep.rownames = "Tissue.code")

wilcox_list_Grammica[["C_indecora__Lutein.epoxide"]]

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  data_loop_wilcox$FW.norm and data_loop_wilcox$Tissue.code 
# 
#     y    o    h   
#   o -    -    -   
#   h 0.14 0.14 -   
#   f 0.39 0.39 0.23
# 
# P value adjustment method: BH 

# make df
cbd_Lutein.epoxide_C_indecora_Tissue.code <- c("y", "o", "h", "f")
cbd_Lutein.epoxide_C_indecora_Letters <- c("a", "a", "a", "a")
cbd_Lutein.epoxide_C_indecora <- data.frame("Tissue.code" = cbd_Lutein.epoxide_C_indecora_Tissue.code, "Letter" = cbd_Lutein.epoxide_C_indecora_Letters)


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein.epoxide_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Lutein.epoxide_C_indecora
# now join these positions to cbd
left_join(cbd_Lutein.epoxide_C_indecora, top_positions_Lutein.epoxide_C_indecora, by = "Tissue.code") -> cbd_Lutein.epoxide_C_indecora

# calculate how much to nudge
data_C_indecora_Lutein.epoxide %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein.epoxide_C_indecora
cbd_Lutein.epoxide_C_indecora$nudged <- (max_Lutein.epoxide_C_indecora$max * 1.05)

# add CLDs to plot
Lutein.epoxide_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein.epoxide_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein.epoxide_C_indecora_boxplot


Lutein.epoxide_C_indecora_boxplot





#### Lutein loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Lutein
data_Lutein_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Lut.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Lutein_plots_grammicaspe$Tissue.code <- factor(data_Lutein_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Lutein = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Grammica_Lutein[[sub]] = p
  
}


#### Lutein C_polygonorum ####
plot_list_Grammica_Lutein[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_polygonorum_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_polygonorum")

box.rslt_Lutein_C_polygonorum <- with(data_C_polygonorum_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_polygonorum)
boxplot_positions_Lutein_C_polygonorum <- as.data.frame(box.rslt_Lutein_C_polygonorum$stats)

# what are these column tissue codes?
tissues_Lutein_C_polygonorum <- levels(data_C_polygonorum_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_polygonorum) <- tissues_Lutein_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_polygonorum <- boxplot_positions_Lutein_C_polygonorum[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_polygonorum, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_polygonorum %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_polygonorum
# now join these positions to cbd
left_join(cbd_Lutein_C_polygonorum, top_positions_Lutein_C_polygonorum, by = "Tissue.code") -> cbd_Lutein_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_polygonorum
cbd_Lutein_C_polygonorum$nudged <- max_Lutein_C_polygonorum$max * 1.05


# add CLDs to plot
Lutein_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_polygonorum_boxplot


Lutein_C_polygonorum_boxplot



#### Lutein C_sandwichiana ####
plot_list_Grammica_Lutein[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_sandwichiana_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_Lutein_C_sandwichiana <- with(data_C_sandwichiana_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_sandwichiana)
boxplot_positions_Lutein_C_sandwichiana <- as.data.frame(box.rslt_Lutein_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_Lutein_C_sandwichiana <- levels(data_C_sandwichiana_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_sandwichiana) <- tissues_Lutein_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_sandwichiana <- boxplot_positions_Lutein_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_sandwichiana, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_sandwichiana %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_sandwichiana
# now join these positions to cbd
left_join(cbd_Lutein_C_sandwichiana, top_positions_Lutein_C_sandwichiana, by = "Tissue.code") -> cbd_Lutein_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_sandwichiana
cbd_Lutein_C_sandwichiana$nudged <- max_Lutein_C_sandwichiana$max * 1.05


# add CLDs to plot
Lutein_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_sandwichiana_boxplot


Lutein_C_sandwichiana_boxplot



#### Lutein C_californica ####
plot_list_Grammica_Lutein[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_californica")

box.rslt_Lutein_C_californica <- with(data_C_californica_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_californica)
boxplot_positions_Lutein_C_californica <- as.data.frame(box.rslt_Lutein_C_californica$stats)

# what are these column tissue codes?
tissues_Lutein_C_californica <- levels(data_C_californica_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_californica) <- tissues_Lutein_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_californica <- boxplot_positions_Lutein_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_californica
# now join these positions to cbd
left_join(cbd_Lutein_C_californica, top_positions_Lutein_C_californica, by = "Tissue.code") -> cbd_Lutein_C_californica

# calculate how much to nudge
data_C_californica_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_californica
cbd_Lutein_C_californica$nudged <- max_Lutein_C_californica$max * 1.05


# add CLDs to plot
Lutein_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_californica_boxplot


Lutein_C_californica_boxplot


#### Lutein C_compacta ####
plot_list_Grammica_Lutein[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_compacta")

box.rslt_Lutein_C_compacta <- with(data_C_compacta_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_compacta)
boxplot_positions_Lutein_C_compacta <- as.data.frame(box.rslt_Lutein_C_compacta$stats)

# what are these column tissue codes?
tissues_Lutein_C_compacta <- levels(data_C_compacta_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_compacta) <- tissues_Lutein_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_compacta <- boxplot_positions_Lutein_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_compacta
# now join these positions to cbd
left_join(cbd_Lutein_C_compacta, top_positions_Lutein_C_compacta, by = "Tissue.code") -> cbd_Lutein_C_compacta

# calculate how much to nudge
data_C_compacta_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_compacta
cbd_Lutein_C_compacta$nudged <- max_Lutein_C_compacta$max * 1.05


# add CLDs to plot
Lutein_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_compacta_boxplot


Lutein_C_compacta_boxplot


#### Lutein C_cephalanthii ####
plot_list_Grammica_Lutein[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Lutein_C_cephalanthii <- with(data_C_cephalanthii_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_cephalanthii)
boxplot_positions_Lutein_C_cephalanthii <- as.data.frame(box.rslt_Lutein_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Lutein_C_cephalanthii <- levels(data_C_cephalanthii_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_cephalanthii) <- tissues_Lutein_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_cephalanthii <- boxplot_positions_Lutein_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Lutein_C_cephalanthii, top_positions_Lutein_C_cephalanthii, by = "Tissue.code") -> cbd_Lutein_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_cephalanthii
cbd_Lutein_C_cephalanthii$nudged <- max_Lutein_C_cephalanthii$max * 1.05


# add CLDs to plot
Lutein_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_cephalanthii_boxplot


Lutein_C_cephalanthii_boxplot


#### Lutein C_denticulata ####
plot_list_Grammica_Lutein[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_denticulata_boxplot


# KW not sig so no post hoc



#### Lutein C_tasmanica ####
plot_list_Grammica_Lutein[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Lutein_C_tasmanica <- with(data_C_tasmanica_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_tasmanica)
boxplot_positions_Lutein_C_tasmanica <- as.data.frame(box.rslt_Lutein_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Lutein_C_tasmanica <- levels(data_C_tasmanica_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_tasmanica) <- tissues_Lutein_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_tasmanica <- boxplot_positions_Lutein_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_tasmanica
# now join these positions to cbd
left_join(cbd_Lutein_C_tasmanica, top_positions_Lutein_C_tasmanica, by = "Tissue.code") -> cbd_Lutein_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_tasmanica
cbd_Lutein_C_tasmanica$nudged <- max_Lutein_C_tasmanica$max * 1.05


# add CLDs to plot
Lutein_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_tasmanica_boxplot


Lutein_C_tasmanica_boxplot


#### Lutein C_costaricensis ####
plot_list_Grammica_Lutein[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Lutein_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Lutein_C_costaricensis <- with(data_C_costaricensis_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_costaricensis)
boxplot_positions_Lutein_C_costaricensis <- as.data.frame(box.rslt_Lutein_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Lutein_C_costaricensis <- levels(data_C_costaricensis_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_costaricensis) <- tissues_Lutein_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_costaricensis <- boxplot_positions_Lutein_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_costaricensis
# now join these positions to cbd
left_join(cbd_Lutein_C_costaricensis, top_positions_Lutein_C_costaricensis, by = "Tissue.code") -> cbd_Lutein_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_costaricensis
cbd_Lutein_C_costaricensis$nudged <- max_Lutein_C_costaricensis$max * 1.05


# add CLDs to plot
Lutein_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_costaricensis_boxplot


Lutein_C_costaricensis_boxplot


#### Lutein C_indecora species alone ####

data_Grammica_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_indecora")

Lutein_C_indecora_boxplot <- ggplot(data_Grammica_Lutein, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Lut.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Lutein <- dplyr::filter(data_Lutein_plots_grammicaspe, Species == "C_indecora")

box.rslt_Lutein_C_indecora <- with(data_C_indecora_Lutein, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Lutein_C_indecora)
boxplot_positions_Lutein_C_indecora <- as.data.frame(box.rslt_Lutein_C_indecora$stats)

# what are these column tissue codes?
tissues_Lutein_C_indecora <- levels(data_C_indecora_Lutein$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Lutein_C_indecora) <- tissues_Lutein_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Lutein_C_indecora <- boxplot_positions_Lutein_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Lutein_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Lut.Car"]])[["Letters"]])
colnames(cbd_Lutein_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Lutein_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Lutein_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Lutein_C_indecora
# now join these positions to cbd
left_join(cbd_Lutein_C_indecora, top_positions_Lutein_C_indecora, by = "Tissue.code") -> cbd_Lutein_C_indecora

# calculate how much to nudge
data_C_indecora_Lutein %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_Lutein_C_indecora
cbd_Lutein_C_indecora$nudged <- (max_Lutein_C_indecora$max * 1.05)

# add CLDs to plot
Lutein_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Lutein_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Lutein_C_indecora_boxplot


Lutein_C_indecora_boxplot



#### a.Carotene loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for a.Carotene
data_a.Carotene_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "a.Car.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_a.Carotene_plots_grammicaspe$Tissue.code <- factor(data_a.Carotene_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_a.Carotene = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_a.Carotene_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Grammica_a.Carotene[[sub]] = p
  
}

#### a.Carotene C_polygonorum ####
plot_list_Grammica_a.Carotene[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_polygonorum_boxplot

# KW not sig so no post hoc
a.Carotene_C_polygonorum_boxplot


#### a.Carotene C_sandwichiana ####
plot_list_Grammica_a.Carotene[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_sandwichiana_boxplot


# KW not sig so no post hoc
a.Carotene_C_sandwichiana_boxplot


#### a.Carotene C_californica ####
plot_list_Grammica_a.Carotene[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_californica_boxplot


# KW not sig so no post hoc
a.Carotene_C_californica_boxplot


#### a.Carotene C_compacta ####
plot_list_Grammica_a.Carotene[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_a.Carotene <- dplyr::filter(data_a.Carotene_plots_grammicaspe, Species == "C_compacta")

box.rslt_a.Carotene_C_compacta <- with(data_C_compacta_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_C_compacta)
boxplot_positions_a.Carotene_C_compacta <- as.data.frame(box.rslt_a.Carotene_C_compacta$stats)

# what are these column tissue codes?
tissues_a.Carotene_C_compacta <- levels(data_C_compacta_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_C_compacta) <- tissues_a.Carotene_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_C_compacta <- boxplot_positions_a.Carotene_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_C_compacta
# now join these positions to cbd
left_join(cbd_a.Carotene_C_compacta, top_positions_a.Carotene_C_compacta, by = "Tissue.code") -> cbd_a.Carotene_C_compacta

# calculate how much to nudge
data_C_compacta_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_C_compacta
cbd_a.Carotene_C_compacta$nudged <- max_a.Carotene_C_compacta$max * 1.05


# add CLDs to plot
a.Carotene_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_a.Carotene_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_C_compacta_boxplot


a.Carotene_C_compacta_boxplot


#### a.Carotene C_cephalanthii ####
plot_list_Grammica_a.Carotene[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_a.Carotene <- dplyr::filter(data_a.Carotene_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_a.Carotene_C_cephalanthii <- with(data_C_cephalanthii_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_C_cephalanthii)
boxplot_positions_a.Carotene_C_cephalanthii <- as.data.frame(box.rslt_a.Carotene_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_a.Carotene_C_cephalanthii <- levels(data_C_cephalanthii_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_C_cephalanthii) <- tissues_a.Carotene_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_C_cephalanthii <- boxplot_positions_a.Carotene_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_C_cephalanthii
# now join these positions to cbd
left_join(cbd_a.Carotene_C_cephalanthii, top_positions_a.Carotene_C_cephalanthii, by = "Tissue.code") -> cbd_a.Carotene_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_C_cephalanthii
cbd_a.Carotene_C_cephalanthii$nudged <- max_a.Carotene_C_cephalanthii$max * 1.05


# add CLDs to plot
a.Carotene_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_a.Carotene_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_C_cephalanthii_boxplot


a.Carotene_C_cephalanthii_boxplot


#### a.Carotene C_denticulata ####
plot_list_Grammica_a.Carotene[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_denticulata_boxplot


# KW not sig so no post hoc
a.Carotene_C_denticulata_boxplot


#### a.Carotene C_tasmanica ####
plot_list_Grammica_a.Carotene[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_tasmanica_boxplot

# KW not sig so no post hoc


#### a.Carotene C_costaricensis ####
plot_list_Grammica_a.Carotene[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> a.Carotene_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_a.Carotene <- dplyr::filter(data_a.Carotene_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_a.Carotene_C_costaricensis <- with(data_C_costaricensis_a.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_a.Carotene_C_costaricensis)
boxplot_positions_a.Carotene_C_costaricensis <- as.data.frame(box.rslt_a.Carotene_C_costaricensis$stats)

# what are these column tissue codes?
tissues_a.Carotene_C_costaricensis <- levels(data_C_costaricensis_a.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_a.Carotene_C_costaricensis) <- tissues_a.Carotene_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_a.Carotene_C_costaricensis <- boxplot_positions_a.Carotene_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_a.Carotene_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__a.Car.Car"]])[["Letters"]])
colnames(cbd_a.Carotene_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_a.Carotene_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_a.Carotene_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_a.Carotene_C_costaricensis
# now join these positions to cbd
left_join(cbd_a.Carotene_C_costaricensis, top_positions_a.Carotene_C_costaricensis, by = "Tissue.code") -> cbd_a.Carotene_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_a.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_a.Carotene_C_costaricensis
cbd_a.Carotene_C_costaricensis$nudged <- max_a.Carotene_C_costaricensis$max * 1.05


# add CLDs to plot
a.Carotene_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_a.Carotene_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> a.Carotene_C_costaricensis_boxplot


a.Carotene_C_costaricensis_boxplot


#### a.Carotene C_indecora species alone ####

data_Grammica_a.Carotene <- dplyr::filter(data_a.Carotene_plots_grammicaspe, Species == "C_indecora")

a.Carotene_C_indecora_boxplot <- ggplot(data_Grammica_a.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "a.Car.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# KW not sig so no post hoc
a.Carotene_C_indecora_boxplot

#### b.Carotene loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for b.Carotene
data_b.Carotene_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "b.Car.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_b.Carotene_plots_grammicaspe$Tissue.code <- factor(data_b.Carotene_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_b.Carotene = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on C_sandwichiana plots
          axis.text.x= element_text(angle = -30, hjust = 0, size = 4), 
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Grammica_b.Carotene[[sub]] = p
  
}

#### b.Carotene C_polygonorum ####
plot_list_Grammica_b.Carotene[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_polygonorum_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_polygonorum")

box.rslt_b.Carotene_C_polygonorum <- with(data_C_polygonorum_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_polygonorum)
boxplot_positions_b.Carotene_C_polygonorum <- as.data.frame(box.rslt_b.Carotene_C_polygonorum$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_polygonorum <- levels(data_C_polygonorum_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_polygonorum) <- tissues_b.Carotene_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_polygonorum <- boxplot_positions_b.Carotene_C_polygonorum[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_polygonorum, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_polygonorum %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_polygonorum
# now join these positions to cbd
left_join(cbd_b.Carotene_C_polygonorum, top_positions_b.Carotene_C_polygonorum, by = "Tissue.code") -> cbd_b.Carotene_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_polygonorum
cbd_b.Carotene_C_polygonorum$nudged <- max_b.Carotene_C_polygonorum$max * 1.05


# add CLDs to plot
b.Carotene_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_b.Carotene_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_polygonorum_boxplot


b.Carotene_C_polygonorum_boxplot



#### b.Carotene C_sandwichiana ####
plot_list_Grammica_b.Carotene[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_sandwichiana_boxplot


# KW not sig so no post hoc
b.Carotene_C_sandwichiana_boxplot


#### b.Carotene C_californica ####
plot_list_Grammica_b.Carotene[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_californica_boxplot


# KW not sig so no post hoc
b.Carotene_C_californica_boxplot


#### b.Carotene C_compacta ####
plot_list_Grammica_b.Carotene[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_compacta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_compacta_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_compacta")

box.rslt_b.Carotene_C_compacta <- with(data_C_compacta_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_compacta)
boxplot_positions_b.Carotene_C_compacta <- as.data.frame(box.rslt_b.Carotene_C_compacta$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_compacta <- levels(data_C_compacta_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_compacta) <- tissues_b.Carotene_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_compacta <- boxplot_positions_b.Carotene_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_compacta
# now join these positions to cbd
left_join(cbd_b.Carotene_C_compacta, top_positions_b.Carotene_C_compacta, by = "Tissue.code") -> cbd_b.Carotene_C_compacta

# calculate how much to nudge
data_C_compacta_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_compacta
cbd_b.Carotene_C_compacta$nudged <- max_b.Carotene_C_compacta$max * 1.05


# add CLDs to plot
b.Carotene_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_b.Carotene_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_compacta_boxplot


b.Carotene_C_compacta_boxplot


#### b.Carotene C_cephalanthii ####
plot_list_Grammica_b.Carotene[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_b.Carotene_C_cephalanthii <- with(data_C_cephalanthii_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_cephalanthii)
boxplot_positions_b.Carotene_C_cephalanthii <- as.data.frame(box.rslt_b.Carotene_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_cephalanthii <- levels(data_C_cephalanthii_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_cephalanthii) <- tissues_b.Carotene_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_cephalanthii <- boxplot_positions_b.Carotene_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_cephalanthii
# now join these positions to cbd
left_join(cbd_b.Carotene_C_cephalanthii, top_positions_b.Carotene_C_cephalanthii, by = "Tissue.code") -> cbd_b.Carotene_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_cephalanthii
cbd_b.Carotene_C_cephalanthii$nudged <- max_b.Carotene_C_cephalanthii$max * 1.05


# add CLDs to plot
b.Carotene_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_b.Carotene_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_cephalanthii_boxplot


b.Carotene_C_cephalanthii_boxplot


#### b.Carotene C_denticulata ####
plot_list_Grammica_b.Carotene[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_denticulata_boxplot


# KW not sig so no post hoc
b.Carotene_C_denticulata_boxplot



#### b.Carotene C_tasmanica ####
plot_list_Grammica_b.Carotene[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_b.Carotene_C_tasmanica <- with(data_C_tasmanica_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_tasmanica)
boxplot_positions_b.Carotene_C_tasmanica <- as.data.frame(box.rslt_b.Carotene_C_tasmanica$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_tasmanica <- levels(data_C_tasmanica_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_tasmanica) <- tissues_b.Carotene_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_tasmanica <- boxplot_positions_b.Carotene_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_tasmanica
# now join these positions to cbd
left_join(cbd_b.Carotene_C_tasmanica, top_positions_b.Carotene_C_tasmanica, by = "Tissue.code") -> cbd_b.Carotene_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_tasmanica
cbd_b.Carotene_C_tasmanica$nudged <- max_b.Carotene_C_tasmanica$max * 1.05


# add CLDs to plot
b.Carotene_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_b.Carotene_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_tasmanica_boxplot


b.Carotene_C_tasmanica_boxplot


#### b.Carotene C_costaricensis ####
plot_list_Grammica_b.Carotene[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> b.Carotene_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_b.Carotene_C_costaricensis <- with(data_C_costaricensis_b.Carotene, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_b.Carotene_C_costaricensis)
boxplot_positions_b.Carotene_C_costaricensis <- as.data.frame(box.rslt_b.Carotene_C_costaricensis$stats)

# what are these column tissue codes?
tissues_b.Carotene_C_costaricensis <- levels(data_C_costaricensis_b.Carotene$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_b.Carotene_C_costaricensis) <- tissues_b.Carotene_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_b.Carotene_C_costaricensis <- boxplot_positions_b.Carotene_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_b.Carotene_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__b.Car.Car"]])[["Letters"]])
colnames(cbd_b.Carotene_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_b.Carotene_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_b.Carotene_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_b.Carotene_C_costaricensis
# now join these positions to cbd
left_join(cbd_b.Carotene_C_costaricensis, top_positions_b.Carotene_C_costaricensis, by = "Tissue.code") -> cbd_b.Carotene_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_b.Carotene %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_b.Carotene_C_costaricensis
cbd_b.Carotene_C_costaricensis$nudged <- max_b.Carotene_C_costaricensis$max * 1.05


# add CLDs to plot
b.Carotene_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_b.Carotene_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> b.Carotene_C_costaricensis_boxplot


b.Carotene_C_costaricensis_boxplot


#### b.Carotene C_indecora species alone ####

data_Grammica_b.Carotene <- dplyr::filter(data_b.Carotene_plots_grammicaspe, Species == "C_indecora")

b.Carotene_C_indecora_boxplot <- ggplot(data_Grammica_b.Carotene, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
  scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
  scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
  geom_boxplot(outlier.size = 0.1, lwd=0.2) +
  theme_minimal() +
  theme(text = element_text(size=10),
        strip.text.x = element_text(angle=0, face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"), 
        axis.text.y = element_text(size = 5), 
        axis.text.x= element_text(angle = -30, hjust = 0, size = 4), 
        axis.ticks.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "b.Car.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# KW not sig so no post hoc
b.Carotene_C_indecora_boxplot


#### Tot.Car loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Tot.Car
data_Tot.Car_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Tot.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Tot.Car_plots_grammicaspe$Tissue.code <- factor(data_Tot.Car_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Tot.Car = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 6))
  plot_list_Grammica_Tot.Car[[sub]] = p
  
}

#### Tot.Car C_polygonorum ####
plot_list_Grammica_Tot.Car[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_polygonorum_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_polygonorum")

box.rslt_Tot.Car_C_polygonorum <- with(data_C_polygonorum_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_polygonorum)
boxplot_positions_Tot.Car_C_polygonorum <- as.data.frame(box.rslt_Tot.Car_C_polygonorum$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_polygonorum <- levels(data_C_polygonorum_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_polygonorum) <- tissues_Tot.Car_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_polygonorum <- boxplot_positions_Tot.Car_C_polygonorum[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_polygonorum, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_polygonorum %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_polygonorum
# now join these positions to cbd
left_join(cbd_Tot.Car_C_polygonorum, top_positions_Tot.Car_C_polygonorum, by = "Tissue.code") -> cbd_Tot.Car_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_polygonorum
cbd_Tot.Car_C_polygonorum$nudged <- max_Tot.Car_C_polygonorum$max * 1.05


# add CLDs to plot
Tot.Car_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_polygonorum_boxplot


Tot.Car_C_polygonorum_boxplot



#### Tot.Car C_sandwichiana ####
plot_list_Grammica_Tot.Car[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_sandwichiana_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_Tot.Car_C_sandwichiana <- with(data_C_sandwichiana_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_sandwichiana)
boxplot_positions_Tot.Car_C_sandwichiana <- as.data.frame(box.rslt_Tot.Car_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_sandwichiana <- levels(data_C_sandwichiana_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_sandwichiana) <- tissues_Tot.Car_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_sandwichiana <- boxplot_positions_Tot.Car_C_sandwichiana[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_sandwichiana, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_sandwichiana %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_sandwichiana
# now join these positions to cbd
left_join(cbd_Tot.Car_C_sandwichiana, top_positions_Tot.Car_C_sandwichiana, by = "Tissue.code") -> cbd_Tot.Car_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_sandwichiana
cbd_Tot.Car_C_sandwichiana$nudged <- max_Tot.Car_C_sandwichiana$max * 1.05


# add CLDs to plot
Tot.Car_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_sandwichiana_boxplot


Tot.Car_C_sandwichiana_boxplot


#### Tot.Car C_californica ####
plot_list_Grammica_Tot.Car[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_californica")

box.rslt_Tot.Car_C_californica <- with(data_C_californica_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_californica)
boxplot_positions_Tot.Car_C_californica <- as.data.frame(box.rslt_Tot.Car_C_californica$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_californica <- levels(data_C_californica_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_californica) <- tissues_Tot.Car_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_californica <- boxplot_positions_Tot.Car_C_californica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_californica
# now join these positions to cbd
left_join(cbd_Tot.Car_C_californica, top_positions_Tot.Car_C_californica, by = "Tissue.code") -> cbd_Tot.Car_C_californica

# calculate how much to nudge
data_C_californica_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_californica
cbd_Tot.Car_C_californica$nudged <- max_Tot.Car_C_californica$max * 1.05


# add CLDs to plot
Tot.Car_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_californica_boxplot


Tot.Car_C_californica_boxplot


#### Tot.Car C_compacta ####
plot_list_Grammica_Tot.Car[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_compacta_boxplot



# use base R boxplot to get the coordinates of the boxes
data_C_compacta_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_compacta")

box.rslt_Tot.Car_C_compacta <- with(data_C_compacta_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_compacta)
boxplot_positions_Tot.Car_C_compacta <- as.data.frame(box.rslt_Tot.Car_C_compacta$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_compacta <- levels(data_C_compacta_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_compacta) <- tissues_Tot.Car_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_compacta <- boxplot_positions_Tot.Car_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_compacta
# now join these positions to cbd
left_join(cbd_Tot.Car_C_compacta, top_positions_Tot.Car_C_compacta, by = "Tissue.code") -> cbd_Tot.Car_C_compacta

# calculate how much to nudge
data_C_compacta_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_compacta
cbd_Tot.Car_C_compacta$nudged <- max_Tot.Car_C_compacta$max * 1.05


# add CLDs to plot
Tot.Car_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_compacta_boxplot


Tot.Car_C_compacta_boxplot


#### Tot.Car C_cephalanthii ####
plot_list_Grammica_Tot.Car[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Tot.Car_C_cephalanthii <- with(data_C_cephalanthii_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_cephalanthii)
boxplot_positions_Tot.Car_C_cephalanthii <- as.data.frame(box.rslt_Tot.Car_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_cephalanthii <- levels(data_C_cephalanthii_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_cephalanthii) <- tissues_Tot.Car_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_cephalanthii <- boxplot_positions_Tot.Car_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Tot.Car_C_cephalanthii, top_positions_Tot.Car_C_cephalanthii, by = "Tissue.code") -> cbd_Tot.Car_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_cephalanthii
cbd_Tot.Car_C_cephalanthii$nudged <- max_Tot.Car_C_cephalanthii$max * 1.05


# add CLDs to plot
Tot.Car_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_cephalanthii_boxplot


Tot.Car_C_cephalanthii_boxplot


#### Tot.Car C_denticulata ####
plot_list_Grammica_Tot.Car[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_denticulata_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_denticulata")

box.rslt_Tot.Car_C_denticulata <- with(data_C_denticulata_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_denticulata)
boxplot_positions_Tot.Car_C_denticulata <- as.data.frame(box.rslt_Tot.Car_C_denticulata$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_denticulata <- levels(data_C_denticulata_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_denticulata) <- tissues_Tot.Car_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_denticulata <- boxplot_positions_Tot.Car_C_denticulata[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_denticulata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_denticulata %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_denticulata
# now join these positions to cbd
left_join(cbd_Tot.Car_C_denticulata, top_positions_Tot.Car_C_denticulata, by = "Tissue.code") -> cbd_Tot.Car_C_denticulata

# calculate how much to nudge
data_C_denticulata_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_denticulata
cbd_Tot.Car_C_denticulata$nudged <- max_Tot.Car_C_denticulata$max * 1.05


# add CLDs to plot
Tot.Car_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_denticulata_boxplot


Tot.Car_C_denticulata_boxplot




#### Tot.Car C_tasmanica ####
plot_list_Grammica_Tot.Car[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Tot.Car_C_tasmanica <- with(data_C_tasmanica_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_tasmanica)
boxplot_positions_Tot.Car_C_tasmanica <- as.data.frame(box.rslt_Tot.Car_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_tasmanica <- levels(data_C_tasmanica_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_tasmanica) <- tissues_Tot.Car_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_tasmanica <- boxplot_positions_Tot.Car_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_tasmanica
# now join these positions to cbd
left_join(cbd_Tot.Car_C_tasmanica, top_positions_Tot.Car_C_tasmanica, by = "Tissue.code") -> cbd_Tot.Car_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_tasmanica
cbd_Tot.Car_C_tasmanica$nudged <- max_Tot.Car_C_tasmanica$max * 1.05


# add CLDs to plot
Tot.Car_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_tasmanica_boxplot


Tot.Car_C_tasmanica_boxplot


#### Tot.Car C_costaricensis ####
plot_list_Grammica_Tot.Car[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Tot.Car_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Tot.Car_C_costaricensis <- with(data_C_costaricensis_Tot.Car, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Tot.Car_C_costaricensis)
boxplot_positions_Tot.Car_C_costaricensis <- as.data.frame(box.rslt_Tot.Car_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Tot.Car_C_costaricensis <- levels(data_C_costaricensis_Tot.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Tot.Car_C_costaricensis) <- tissues_Tot.Car_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Tot.Car_C_costaricensis <- boxplot_positions_Tot.Car_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Tot.Car_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Tot.Car"]])[["Letters"]])
colnames(cbd_Tot.Car_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Tot.Car_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Tot.Car_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Tot.Car_C_costaricensis
# now join these positions to cbd
left_join(cbd_Tot.Car_C_costaricensis, top_positions_Tot.Car_C_costaricensis, by = "Tissue.code") -> cbd_Tot.Car_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Tot.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Tot.Car_C_costaricensis
cbd_Tot.Car_C_costaricensis$nudged <- max_Tot.Car_C_costaricensis$max * 1.05


# add CLDs to plot
Tot.Car_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Tot.Car_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Tot.Car_C_costaricensis_boxplot


Tot.Car_C_costaricensis_boxplot


#### Tot.Car C_indecora species alone ####

data_Grammica_Tot.Car <- dplyr::filter(data_Tot.Car_plots_grammicaspe, Species == "C_indecora")

Tot.Car_C_indecora_boxplot <- ggplot(data_Grammica_Tot.Car, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg FW") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 6)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Tot.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# KW not sig so no post hoc
Tot.Car_C_indecora_boxplot



#### NVZ.Car loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for NVZ.Car
data_NVZ.Car_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "NVZ.Car")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_NVZ.Car_plots_grammicaspe$Tissue.code <- factor(data_NVZ.Car_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_NVZ.Car = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) + 
    scale_fill_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_color_manual(name = "Tissue", labels = c("Seedling", "Young", "Old", "Haustorium", "Flower", "Seed"),values = c("sdlg" = seedling, "l" = leaf, "y" = young, "o" = old, "h" = haustorium, "f" = flower, "s" = seed)) +
    scale_x_discrete(name = "Tissue", labels = c("sdlg" = "Sg", "y" = "Y", "o" = "O", "h" = "H", "f" = "F", "s" = "Sd"), drop = FALSE) +
    geom_boxplot(outlier.size = 0.1, lwd=0.2) +
    theme_minimal() +
    theme(text = element_text(size=10),
          strip.text.x = element_text(angle=0, face = "bold"),
          strip.text.y.left = element_text(angle = 0, face = "bold"), 
          axis.text.y = element_blank(), # uncomment later but need for now to set limits on C_sandwichiana plots
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    guides(colour = guide_legend(nrow = 1)) +
    scale_y_continuous(position = "right", limits = c(0, 100))
  plot_list_Grammica_NVZ.Car[[sub]] = p
  
}

#### NVZ.Car C_polygonorum ####
plot_list_Grammica_NVZ.Car[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_polygonorum_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_polygonorum_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_polygonorum")

box.rslt_NVZ.Car_C_polygonorum <- with(data_C_polygonorum_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_polygonorum)
boxplot_positions_NVZ.Car_C_polygonorum <- as.data.frame(box.rslt_NVZ.Car_C_polygonorum$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_polygonorum <- levels(data_C_polygonorum_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_polygonorum) <- tissues_NVZ.Car_C_polygonorum

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_polygonorum <- boxplot_positions_NVZ.Car_C_polygonorum[5,]


#add pairwise significance letter groups (polygonorum letter display; CLD)
cbd_NVZ.Car_C_polygonorum <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_polygonorum__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_polygonorum)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_polygonorum, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_polygonorum %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_polygonorum
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_polygonorum, top_positions_NVZ.Car_C_polygonorum, by = "Tissue.code") -> cbd_NVZ.Car_C_polygonorum

# calculate how much to nudge
data_C_polygonorum_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_polygonorum
cbd_NVZ.Car_C_polygonorum$nudged <- max_NVZ.Car_C_polygonorum$max * 1.05


# add CLDs to plot
NVZ.Car_C_polygonorum_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_polygonorum,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_polygonorum_boxplot


NVZ.Car_C_polygonorum_boxplot


#### NVZ.Car C_sandwichiana ####
plot_list_Grammica_NVZ.Car[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_sandwichiana_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_sandwichiana_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_sandwichiana")

box.rslt_NVZ.Car_C_sandwichiana <- with(data_C_sandwichiana_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_sandwichiana)
boxplot_positions_NVZ.Car_C_sandwichiana <- as.data.frame(box.rslt_NVZ.Car_C_sandwichiana$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_sandwichiana <- levels(data_C_sandwichiana_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_sandwichiana) <- tissues_NVZ.Car_C_sandwichiana

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_sandwichiana <- boxplot_positions_NVZ.Car_C_sandwichiana[5,]


#add pairwise significance letter groups (polygonorum letter display; CLD)
cbd_NVZ.Car_C_sandwichiana <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_sandwichiana__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_sandwichiana)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_sandwichiana, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_sandwichiana %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_sandwichiana
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_sandwichiana, top_positions_NVZ.Car_C_sandwichiana, by = "Tissue.code") -> cbd_NVZ.Car_C_sandwichiana

# calculate how much to nudge
data_C_sandwichiana_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_sandwichiana
cbd_NVZ.Car_C_sandwichiana$nudged <- max_NVZ.Car_C_sandwichiana$max * 1.05


# add CLDs to plot
NVZ.Car_C_sandwichiana_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_sandwichiana,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_sandwichiana_boxplot


NVZ.Car_C_sandwichiana_boxplot


#### NVZ.Car C_californica ####
plot_list_Grammica_NVZ.Car[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_californica_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_californica_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_californica")

box.rslt_NVZ.Car_C_californica <- with(data_C_californica_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_californica)
boxplot_positions_NVZ.Car_C_californica <- as.data.frame(box.rslt_NVZ.Car_C_californica$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_californica <- levels(data_C_californica_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_californica) <- tissues_NVZ.Car_C_californica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_californica <- boxplot_positions_NVZ.Car_C_californica[5,]


#add pairwise significance letter groups (polygonorum letter display; CLD)
cbd_NVZ.Car_C_californica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_californica__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_californica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_californica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_californica %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_californica
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_californica, top_positions_NVZ.Car_C_californica, by = "Tissue.code") -> cbd_NVZ.Car_C_californica

# calculate how much to nudge
data_C_californica_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_californica
cbd_NVZ.Car_C_californica$nudged <- max_NVZ.Car_C_californica$max * 1.05


# add CLDs to plot
NVZ.Car_C_californica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_californica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_californica_boxplot


NVZ.Car_C_californica_boxplot


#### NVZ.Car C_compacta ####
plot_list_Grammica_NVZ.Car[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_compacta_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_compacta_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_compacta")

box.rslt_NVZ.Car_C_compacta <- with(data_C_compacta_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_compacta)
boxplot_positions_NVZ.Car_C_compacta <- as.data.frame(box.rslt_NVZ.Car_C_compacta$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_compacta <- levels(data_C_compacta_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_compacta) <- tissues_NVZ.Car_C_compacta

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_compacta <- boxplot_positions_NVZ.Car_C_compacta[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_C_compacta <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_compacta__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_compacta)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_compacta, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_compacta %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_compacta
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_compacta, top_positions_NVZ.Car_C_compacta, by = "Tissue.code") -> cbd_NVZ.Car_C_compacta

# calculate how much to nudge
data_C_compacta_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_compacta
cbd_NVZ.Car_C_compacta$nudged <- max_NVZ.Car_C_compacta$max * 1.05


# add CLDs to plot
NVZ.Car_C_compacta_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_compacta,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_compacta_boxplot


NVZ.Car_C_compacta_boxplot


#### NVZ.Car C_cephalanthii ####
plot_list_Grammica_NVZ.Car[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_NVZ.Car_C_cephalanthii <- with(data_C_cephalanthii_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_cephalanthii)
boxplot_positions_NVZ.Car_C_cephalanthii <- as.data.frame(box.rslt_NVZ.Car_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_cephalanthii <- levels(data_C_cephalanthii_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_cephalanthii) <- tissues_NVZ.Car_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_cephalanthii <- boxplot_positions_NVZ.Car_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_cephalanthii
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_cephalanthii, top_positions_NVZ.Car_C_cephalanthii, by = "Tissue.code") -> cbd_NVZ.Car_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_cephalanthii
cbd_NVZ.Car_C_cephalanthii$nudged <- max_NVZ.Car_C_cephalanthii$max * 1.05


# add CLDs to plot
NVZ.Car_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_cephalanthii_boxplot


NVZ.Car_C_cephalanthii_boxplot


#### NVZ.Car C_denticulata ####
plot_list_Grammica_NVZ.Car[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_denticulata_boxplot


# use base R boxplot to get the coordinates of the boxes
data_C_denticulata_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_denticulata")

box.rslt_NVZ.Car_C_denticulata <- with(data_C_denticulata_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_denticulata)
boxplot_positions_NVZ.Car_C_denticulata <- as.data.frame(box.rslt_NVZ.Car_C_denticulata$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_denticulata <- levels(data_C_denticulata_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_denticulata) <- tissues_NVZ.Car_C_denticulata

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_denticulata <- boxplot_positions_NVZ.Car_C_denticulata[5,]


#add pairwise significance letter groups (polygonorum letter display; CLD)
cbd_NVZ.Car_C_denticulata <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_denticulata__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_denticulata)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_denticulata, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_denticulata %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_denticulata
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_denticulata, top_positions_NVZ.Car_C_denticulata, by = "Tissue.code") -> cbd_NVZ.Car_C_denticulata

# calculate how much to nudge
data_C_denticulata_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_denticulata
cbd_NVZ.Car_C_denticulata$nudged <- max_NVZ.Car_C_denticulata$max * 1.05


# add CLDs to plot
NVZ.Car_C_denticulata_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_denticulata,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_denticulata_boxplot


NVZ.Car_C_denticulata_boxplot


#### NVZ.Car C_tasmanica ####
plot_list_Grammica_NVZ.Car[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_NVZ.Car_C_tasmanica <- with(data_C_tasmanica_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_tasmanica)
boxplot_positions_NVZ.Car_C_tasmanica <- as.data.frame(box.rslt_NVZ.Car_C_tasmanica$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_tasmanica <- levels(data_C_tasmanica_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_tasmanica) <- tissues_NVZ.Car_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_tasmanica <- boxplot_positions_NVZ.Car_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_tasmanica
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_tasmanica, top_positions_NVZ.Car_C_tasmanica, by = "Tissue.code") -> cbd_NVZ.Car_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_tasmanica
cbd_NVZ.Car_C_tasmanica$nudged <- max_NVZ.Car_C_tasmanica$max * 1.05


# add CLDs to plot
NVZ.Car_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_tasmanica_boxplot


NVZ.Car_C_tasmanica_boxplot


#### NVZ.Car C_costaricensis ####
plot_list_Grammica_NVZ.Car[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> NVZ.Car_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_NVZ.Car_C_costaricensis <- with(data_C_costaricensis_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_costaricensis)
boxplot_positions_NVZ.Car_C_costaricensis <- as.data.frame(box.rslt_NVZ.Car_C_costaricensis$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_costaricensis <- levels(data_C_costaricensis_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_costaricensis) <- tissues_NVZ.Car_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_costaricensis <- boxplot_positions_NVZ.Car_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_NVZ.Car_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_costaricensis
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_costaricensis, top_positions_NVZ.Car_C_costaricensis, by = "Tissue.code") -> cbd_NVZ.Car_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_costaricensis
cbd_NVZ.Car_C_costaricensis$nudged <- max_NVZ.Car_C_costaricensis$max * 1.05


# add CLDs to plot
NVZ.Car_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_costaricensis_boxplot


NVZ.Car_C_costaricensis_boxplot


#### NVZ.Car C_indecora species alone ####

data_Grammica_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_indecora")

NVZ.Car_C_indecora_boxplot <- ggplot(data_Grammica_NVZ.Car, aes(x=Tissue.code, y=(FW.norm*100), color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("% of Carotenoids") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 100)) +
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "NVZ.Car" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

# use base R boxplot to get the coordinates of the boxes
data_C_indecora_NVZ.Car <- dplyr::filter(data_NVZ.Car_plots_grammicaspe, Species == "C_indecora")

box.rslt_NVZ.Car_C_indecora <- with(data_C_indecora_NVZ.Car, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_NVZ.Car_C_indecora)
boxplot_positions_NVZ.Car_C_indecora <- as.data.frame(box.rslt_NVZ.Car_C_indecora$stats)

# what are these column tissue codes?
tissues_NVZ.Car_C_indecora <- levels(data_C_indecora_NVZ.Car$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_NVZ.Car_C_indecora) <- tissues_NVZ.Car_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_NVZ.Car_C_indecora <- boxplot_positions_NVZ.Car_C_indecora[5,]


#add pairwise significance letter groups (polygonorum letter display; CLD)
cbd_NVZ.Car_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__NVZ.Car"]])[["Letters"]])
colnames(cbd_NVZ.Car_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_NVZ.Car_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_NVZ.Car_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_NVZ.Car_C_indecora
# now join these positions to cbd
left_join(cbd_NVZ.Car_C_indecora, top_positions_NVZ.Car_C_indecora, by = "Tissue.code") -> cbd_NVZ.Car_C_indecora

# calculate how much to nudge
data_C_indecora_NVZ.Car %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(FW.norm*100)) -> max_NVZ.Car_C_indecora
cbd_NVZ.Car_C_indecora$nudged <- max_NVZ.Car_C_indecora$max * 1.05


# add CLDs to plot
NVZ.Car_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_NVZ.Car_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> NVZ.Car_C_indecora_boxplot


NVZ.Car_C_indecora_boxplot







#### COMBINED (faceted) carotenoid plot: Grammica ONLY ####
wrap_elements(gridtext::richtext_grob('Total carotenoids', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  Tot.Car_C_australis_boxplot + ggtitle('C. australis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Tot.Car_C_polygonorum_boxplot + ggtitle('C. polygonorum') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Tot.Car_C_sandwichiana_boxplot + ggtitle('C. sandwichiana') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Tot.Car_C_californica_boxplot + ggtitle('C. californica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Tot.Car_C_compacta_boxplot + ggtitle('C. compacta') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Tot.Car_C_cephalanthii_boxplot + ggtitle('C. cephalanthii') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Tot.Car_C_denticulata_boxplot + ggtitle('C. denticulata') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Tot.Car_C_tasmanica_boxplot + ggtitle('C. tasmanica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Tot.Car_C_costaricensis_boxplot + ggtitle('C. costaricensis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Tot.Car_C_indecora_boxplot + ggtitle('C. indecora') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  wrap_elements(gridtext::richtext_grob('VAZ', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  VAZ.Car_C_australis_boxplot + 
  VAZ.Car_C_polygonorum_boxplot + 
  VAZ.Car_C_sandwichiana_boxplot +  
  VAZ.Car_C_californica_boxplot +  
  VAZ.Car_C_compacta_boxplot + 
  VAZ.Car_C_cephalanthii_boxplot + 
  VAZ.Car_C_denticulata_boxplot + 
  VAZ.Car_C_tasmanica_boxplot + 
  VAZ.Car_C_costaricensis_boxplot + 
  VAZ.Car_C_indecora_boxplot + 
  # wrap_elements(gridtext::richtext_grob('NVZ.Car', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  # NVZ.Car_C_australis_boxplot +
  # NVZ.Car_C_polygonorum_boxplot + 
  # NVZ.Car_C_sandwichiana_boxplot +
  # NVZ.Car_C_californica_boxplot + 
  # NVZ.Car_C_compacta_boxplot +
  # NVZ.Car_C_cephalanthii_boxplot +
  # NVZ.Car_C_denticulata_boxplot +
  # NVZ.Car_C_tasmanica_boxplot +
  # NVZ.Car_C_costaricensis_boxplot +
  # NVZ.Car_C_indecora_boxplot +
wrap_elements(gridtext::richtext_grob('Neoxanthin', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) +
  Neoxanthin_C_australis_boxplot +
  Neoxanthin_C_polygonorum_boxplot +
  Neoxanthin_C_sandwichiana_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_californica_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_compacta_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_cephalanthii_boxplot +
  Neoxanthin_C_denticulata_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_tasmanica_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_costaricensis_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  Neoxanthin_C_indecora_boxplot + geom_text(size    = 1.8, color = "black", data    = absent, inherit.aes = T, mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"), nudge_y = -.5, parse = TRUE) +
  wrap_elements(gridtext::richtext_grob('Lutein', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  Lutein_C_australis_boxplot +
  Lutein_C_polygonorum_boxplot + 
  Lutein_C_sandwichiana_boxplot +
  Lutein_C_californica_boxplot + 
  Lutein_C_compacta_boxplot +
  Lutein_C_cephalanthii_boxplot +
  Lutein_C_denticulata_boxplot +
  Lutein_C_tasmanica_boxplot +
  Lutein_C_costaricensis_boxplot +
  Lutein_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('Lutein epoxide', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  Lutein.epoxide_C_australis_boxplot +
  Lutein.epoxide_C_polygonorum_boxplot + 
  Lutein.epoxide_C_sandwichiana_boxplot +
  Lutein.epoxide_C_californica_boxplot + 
  Lutein.epoxide_C_compacta_boxplot +
  Lutein.epoxide_C_cephalanthii_boxplot +
  Lutein.epoxide_C_denticulata_boxplot +
  Lutein.epoxide_C_tasmanica_boxplot +
  Lutein.epoxide_C_costaricensis_boxplot +
  Lutein.epoxide_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('*a*-Carotene', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  a.Carotene_C_australis_boxplot +
  a.Carotene_C_polygonorum_boxplot + 
  a.Carotene_C_sandwichiana_boxplot +
  a.Carotene_C_californica_boxplot + 
  a.Carotene_C_compacta_boxplot +
  a.Carotene_C_cephalanthii_boxplot +
  a.Carotene_C_denticulata_boxplot +
  a.Carotene_C_tasmanica_boxplot +
  a.Carotene_C_costaricensis_boxplot +
  a.Carotene_C_indecora_boxplot +
  wrap_elements(gridtext::richtext_grob('*b*-Carotene', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  b.Carotene_C_australis_boxplot +
  b.Carotene_C_polygonorum_boxplot + 
  b.Carotene_C_sandwichiana_boxplot +
  b.Carotene_C_californica_boxplot + 
  b.Carotene_C_compacta_boxplot +
  b.Carotene_C_cephalanthii_boxplot +
  b.Carotene_C_denticulata_boxplot +
  b.Carotene_C_tasmanica_boxplot +
  b.Carotene_C_costaricensis_boxplot +
  b.Carotene_C_indecora_boxplot +
  plot_layout(nrow = 7, byrow = T) -> carotenoid_boxplot_Grammica_new

carotenoid_boxplot_Grammica_new 

pdf("../output/boxplots/carotenoid_boxplot_Grammica_new.pdf", width=9,height=5.5) 
carotenoid_boxplot_Grammica_new
dev.off()


#### Car.Chl PLOTS! ####
#### Car.Chl31 Ipomoea_nil ####
data_ipomoea_Car.Chl31 <- dplyr::filter(data_ipomoea, Pigment == "Car.Chl31")

Car.Chl31_ipomoea_boxplot <- ggplot(data_ipomoea_Car.Chl31, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
  scale_fill_manual(name = "Tissue", labels = c("Leaf", "Young", "Old", "Flower", "Seed"),values = c("l" = leaf, "y" = young, "o" = old, "f" = flower, "s" = seed)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 mmol/mol") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0,15)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Car.Chl31" & Subgenus == "Ipomoea nil"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 



# use base R boxplot to get the coordinates of the boxes
box.rslt_Car.Chl31_ipomoea <- with(data_ipomoea_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_ipomoea)
boxplot_positions_Car.Chl31_ipomoea <- as.data.frame(box.rslt_Car.Chl31_ipomoea$stats)

# what are these column tissue codes?
tissues_Car.Chl31_ipomoea <- levels(data_ipomoea_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_ipomoea) <- tissues_Car.Chl31_ipomoea

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_ipomoea <- boxplot_positions_Car.Chl31_ipomoea[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_ipomoea <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Ipomoea_nil__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_ipomoea)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_ipomoea, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_ipomoea %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_ipomoea
# now join these positions to cbd
left_join(cbd_Car.Chl31_ipomoea, top_positions_Car.Chl31_ipomoea, by = "Tissue.code") -> cbd_Car.Chl31_ipomoea

# calculate how much to nudge
data_ipomoea_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_ipomoea
cbd_Car.Chl31_ipomoea$nudged <- max_Car.Chl31_ipomoea$max* 1.05


# add CLDs to plot
Car.Chl31_ipomoea_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Car.Chl31_ipomoea,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_ipomoea_boxplot


Car.Chl31_ipomoea_boxplot


#### Car.Chl31 loop through Monogynella, Cuscuta, and C. purpurata ####
loop_subgenera <- c("Monogynella","Cuscuta", "C_purpurata")
# dplyr::filter for Car.Chl31
data_Car.Chl31_plots_cuscutasub <- dplyr::filter(data_long_calcs, Pigment == "Car.Chl31")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Car.Chl31_plots_cuscutasub$Tissue.code <- factor(data_Car.Chl31_plots_cuscutasub$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))

plot_list_Cuscuta_Car.Chl31 = list()

for (sub in (loop_subgenera)) {
  
  data_loop <- dplyr::filter(data_Car.Chl31_plots_cuscutasub, Subgenus == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0,15))
  plot_list_Cuscuta_Car.Chl31[[sub]] = p
  
}

#### Car.Chl31. Monogynella ####
plot_list_Cuscuta_Car.Chl31[["Monogynella"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Car.Chl31" & Subgenus == "Monogynella"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_Monogynella_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Monogynella_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_cuscutasub, Subgenus == "Monogynella")

box.rslt_Car.Chl31_Monogynella <- with(data_Monogynella_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_Monogynella)
boxplot_positions_Car.Chl31_Monogynella <- as.data.frame(box.rslt_Car.Chl31_Monogynella$stats)

# what are these column tissue codes?
tissues_Car.Chl31_Monogynella <- levels(data_Monogynella_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_Monogynella) <- tissues_Car.Chl31_Monogynella

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_Monogynella <- boxplot_positions_Car.Chl31_Monogynella[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_Monogynella <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Monogynella__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_Monogynella)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_Monogynella, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_Monogynella %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_Monogynella
# now join these positions to cbd
left_join(cbd_Car.Chl31_Monogynella, top_positions_Car.Chl31_Monogynella, by = "Tissue.code") -> cbd_Car.Chl31_Monogynella

# calculate how much to nudge
data_Monogynella_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_Monogynella
cbd_Car.Chl31_Monogynella$nudged <- max_Car.Chl31_Monogynella$max * 1.05


# add CLDs to plot
Car.Chl31_Monogynella_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Car.Chl31_Monogynella,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_Monogynella_boxplot


Car.Chl31_Monogynella_boxplot


#### Car.Chl31 Cuscuta ####
plot_list_Cuscuta_Car.Chl31[["Cuscuta"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Car.Chl31" & Subgenus == "Cuscuta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_Cuscuta_boxplot

# KW not sig so no post hoc

Car.Chl31_Cuscuta_boxplot



#### Car.Chl31 C_purpurata ####
plot_list_Cuscuta_Car.Chl31[["C_purpurata"]] + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Car.Chl31" & Subgenus == "C. purpurata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_purpurata_boxplot

# KW not sig so no post hoc

Car.Chl31_C_purpurata_boxplot




#### Car.Chl31 Grammica subgenus alone ####

data_Grammica_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_cuscutasub, Subgenus == "Grammica")

Car.Chl31_Grammica_boxplot <- ggplot(data_Grammica_Car.Chl31, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 mmol/mol") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0,15)) + 
  geom_text(
    size    = 2,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal, Pigment == "Car.Chl31" & Subgenus == "Grammica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 

Car.Chl31_Grammica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_Grammica_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_cuscutasub, Subgenus == "Grammica")

box.rslt_Car.Chl31_Grammica <- with(data_Grammica_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_Grammica)
boxplot_positions_Car.Chl31_Grammica <- as.data.frame(box.rslt_Car.Chl31_Grammica$stats)

# what are these column tissue codes?
tissues_Car.Chl31_Grammica <- levels(data_Grammica_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_Grammica) <- tissues_Car.Chl31_Grammica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_Grammica <- boxplot_positions_Car.Chl31_Grammica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_Grammica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list[["Grammica__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_Grammica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_Grammica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_Grammica %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_Grammica
# now join these positions to cbd
left_join(cbd_Car.Chl31_Grammica, top_positions_Car.Chl31_Grammica, by = "Tissue.code") -> cbd_Car.Chl31_Grammica

# calculate how much to nudge
data_Grammica_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_Grammica
cbd_Car.Chl31_Grammica$nudged <- max_Car.Chl31_Grammica$max * 1.05

# add CLDs to plot
Car.Chl31_Grammica_boxplot + 
  geom_text(
    size    = 2,
    color = "black",
    data    = cbd_Car.Chl31_Grammica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_Grammica_boxplot


Car.Chl31_Grammica_boxplot


#### COMBINED (faceted) car:chl plot: subgenus ####
wrap_elements(gridtext::richtext_grob('Carotenoid:Chlorophyll', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 8, fontface = 'bold'))) + 
  Car.Chl31_ipomoea_boxplot + ggtitle('Ipomoea nil') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Car.Chl31_Monogynella_boxplot + ggtitle('Monogynella') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Car.Chl31_Cuscuta_boxplot + ggtitle('Cuscuta') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Car.Chl31_C_purpurata_boxplot + ggtitle('C. purpurata') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + 
  Car.Chl31_Grammica_boxplot + ggtitle('Grammica') + theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold.italic")) + plot_layout(nrow = 1, byrow = T) -> car.chl_boxplot_new

car.chl_boxplot_new 

pdf("../output/boxplots/car.chl_boxplot_new.pdf", width=7,height=2) 
car.chl_boxplot_new
dev.off()

#### car:chl plot Grammica ONLY ####

#### Car.Chl31 C_australis ####
data_C_australis <-  dplyr::filter(data_long_calcs_Grammica_plot, Species == "C_australis")
data_C_australis_Car.Chl31 <- dplyr::filter(data_C_australis, Pigment == "Car.Chl31")

# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_C_australis$Tissue.code <- factor(data_C_australis$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))

Car.Chl31_C_australis_boxplot <- ggplot(data_C_australis_Car.Chl31, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 ng/mg FW") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "left", limits = c(0, 15)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_australis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.6, parse = TRUE) 

Car.Chl31_C_australis_boxplot 

# KW not sig so no post hoc
Car.Chl31_C_australis_boxplot

#### Car.Chl31 loop through C_polygonorum, C_sandwichiana, C_californica, C_compacta, C_cephalanthii, C_denticulata, C_tasmanica, C_costaricensis, and C. indecora ####
loop_species <- c("C_polygonorum", "C_sandwichiana", "C_californica", "C_compacta", "C_cephalanthii", "C_denticulata", "C_tasmanica", "C_costaricensis", "C_indecora")
# dplyr::filter for Car.Chl31
data_Car.Chl31_plots_grammicaspe <- dplyr::filter(data_long_calcs_Grammica_plot, Pigment == "Car.Chl31")
# drop unused factor levels from tissues (e.g. haustorium from Ipomoea)
data_Car.Chl31_plots_grammicaspe$Tissue.code <- factor(data_Car.Chl31_plots_grammicaspe$Tissue.code, levels = c("sdlg", "y", "o", "h", "f", "s"))


plot_list_Grammica_Car.Chl31 = list()

for (sub in (loop_species)) {
  
  data_loop <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == sub)
  p <- ggplot(data_loop, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) + 
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
    scale_y_continuous(position = "right", limits = c(0, 15))
  plot_list_Grammica_Car.Chl31[[sub]] = p
  
}

#### Car.Chl31 C_polygonorum ####
plot_list_Grammica_Car.Chl31[["C_polygonorum"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_polygonorum"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_polygonorum_boxplot

# KW not sig so no post hoc
Car.Chl31_C_polygonorum_boxplot


#### Car.Chl31 C_sandwichiana ####
plot_list_Grammica_Car.Chl31[["C_sandwichiana"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_sandwichiana"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_sandwichiana_boxplot


# KW not sig so no post hoc
Car.Chl31_C_sandwichiana_boxplot


#### Car.Chl31 C_californica ####
plot_list_Grammica_Car.Chl31[["C_californica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_californica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_californica_boxplot


# KW not sig so no post hoc

Car.Chl31_C_californica_boxplot


#### Car.Chl31 C_compacta ####
plot_list_Grammica_Car.Chl31[["C_compacta"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_compacta"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_compacta_boxplot



# KW not sig so no post hoc

Car.Chl31_C_compacta_boxplot


#### Car.Chl31 C_cephalanthii ####
plot_list_Grammica_Car.Chl31[["C_cephalanthii"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_cephalanthii"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_cephalanthii_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_cephalanthii_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == "C_cephalanthii")

box.rslt_Car.Chl31_C_cephalanthii <- with(data_C_cephalanthii_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_C_cephalanthii)
boxplot_positions_Car.Chl31_C_cephalanthii <- as.data.frame(box.rslt_Car.Chl31_C_cephalanthii$stats)

# what are these column tissue codes?
tissues_Car.Chl31_C_cephalanthii <- levels(data_C_cephalanthii_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_C_cephalanthii) <- tissues_Car.Chl31_C_cephalanthii

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_C_cephalanthii <- boxplot_positions_Car.Chl31_C_cephalanthii[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_C_cephalanthii <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_cephalanthii__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_C_cephalanthii)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_C_cephalanthii, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_C_cephalanthii %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_C_cephalanthii
# now join these positions to cbd
left_join(cbd_Car.Chl31_C_cephalanthii, top_positions_Car.Chl31_C_cephalanthii, by = "Tissue.code") -> cbd_Car.Chl31_C_cephalanthii

# calculate how much to nudge
data_C_cephalanthii_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_C_cephalanthii
cbd_Car.Chl31_C_cephalanthii$nudged <- max_Car.Chl31_C_cephalanthii$max * 1.05


# add CLDs to plot
Car.Chl31_C_cephalanthii_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Car.Chl31_C_cephalanthii,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_C_cephalanthii_boxplot


Car.Chl31_C_cephalanthii_boxplot


#### Car.Chl31 C_denticulata ####
plot_list_Grammica_Car.Chl31[["C_denticulata"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_denticulata"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_denticulata_boxplot


# KW not sig so no post hoc

Car.Chl31_C_denticulata_boxplot



#### Car.Chl31 C_tasmanica ####
plot_list_Grammica_Car.Chl31[["C_tasmanica"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_tasmanica"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_tasmanica_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_tasmanica_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == "C_tasmanica")

box.rslt_Car.Chl31_C_tasmanica <- with(data_C_tasmanica_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_C_tasmanica)
boxplot_positions_Car.Chl31_C_tasmanica <- as.data.frame(box.rslt_Car.Chl31_C_tasmanica$stats)

# what are these column tissue codes?
tissues_Car.Chl31_C_tasmanica <- levels(data_C_tasmanica_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_C_tasmanica) <- tissues_Car.Chl31_C_tasmanica

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_C_tasmanica <- boxplot_positions_Car.Chl31_C_tasmanica[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_C_tasmanica <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_tasmanica__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_C_tasmanica)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_C_tasmanica, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_C_tasmanica %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_C_tasmanica
# now join these positions to cbd
left_join(cbd_Car.Chl31_C_tasmanica, top_positions_Car.Chl31_C_tasmanica, by = "Tissue.code") -> cbd_Car.Chl31_C_tasmanica

# calculate how much to nudge
data_C_tasmanica_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_C_tasmanica
cbd_Car.Chl31_C_tasmanica$nudged <- max_Car.Chl31_C_tasmanica$max * 1.05


# add CLDs to plot
Car.Chl31_C_tasmanica_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Car.Chl31_C_tasmanica,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_C_tasmanica_boxplot


Car.Chl31_C_tasmanica_boxplot


#### Car.Chl31 C_costaricensis ####
plot_list_Grammica_Car.Chl31[["C_costaricensis"]] + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_costaricensis"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) -> Car.Chl31_C_costaricensis_boxplot

# use base R boxplot to get the coordinates of the boxes
data_C_costaricensis_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == "C_costaricensis")

box.rslt_Car.Chl31_C_costaricensis <- with(data_C_costaricensis_Car.Chl31, graphics::boxplot(logFW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_C_costaricensis)
boxplot_positions_Car.Chl31_C_costaricensis <- as.data.frame(box.rslt_Car.Chl31_C_costaricensis$stats)

# what are these column tissue codes?
tissues_Car.Chl31_C_costaricensis <- levels(data_C_costaricensis_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_C_costaricensis) <- tissues_Car.Chl31_C_costaricensis

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_C_costaricensis <- boxplot_positions_Car.Chl31_C_costaricensis[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_C_costaricensis <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_costaricensis__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_C_costaricensis)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_C_costaricensis, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_C_costaricensis %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_C_costaricensis
# now join these positions to cbd
left_join(cbd_Car.Chl31_C_costaricensis, top_positions_Car.Chl31_C_costaricensis, by = "Tissue.code") -> cbd_Car.Chl31_C_costaricensis

# calculate how much to nudge
data_C_costaricensis_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_C_costaricensis
cbd_Car.Chl31_C_costaricensis$nudged <- max_Car.Chl31_C_costaricensis$max * 1.05


# add CLDs to plot
Car.Chl31_C_costaricensis_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Car.Chl31_C_costaricensis,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_C_costaricensis_boxplot


Car.Chl31_C_costaricensis_boxplot


#### Car.Chl31 C_indecora species alone ####

data_Grammica_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == "C_indecora")

Car.Chl31_C_indecora_boxplot <- ggplot(data_Grammica_Car.Chl31, aes(x=Tissue.code, y=logFW.norm, color=Tissue.code)) +
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
        axis.title.y = element_text(size = 5),
        legend.position = "none") +
  ylab("Log10 mmol/mol") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(position = "right", limits = c(0, 15)) + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = dplyr::filter(dat_text_plot_kruskal_Grammica, Pigment == "Car.Chl31" & Species == "C_indecora"),
    inherit.aes = T,
    mapping = aes(x = Inf, y = Inf, label = label, vjust = "top", hjust = "right"),
    nudge_y = -.5, parse = TRUE) 


# use base R boxplot to get the coordinates of the boxes
data_C_indecora_Car.Chl31 <- dplyr::filter(data_Car.Chl31_plots_grammicaspe, Species == "C_indecora")

box.rslt_Car.Chl31_C_indecora <- with(data_C_indecora_Car.Chl31, graphics::boxplot(FW.norm ~ Tissue.code, plot = FALSE))
str(box.rslt_Car.Chl31_C_indecora)
boxplot_positions_Car.Chl31_C_indecora <- as.data.frame(box.rslt_Car.Chl31_C_indecora$stats)

# what are these column tissue codes?
tissues_Car.Chl31_C_indecora <- levels(data_C_indecora_Car.Chl31$Tissue.code)
# add appropriate tissues to position df
colnames(boxplot_positions_Car.Chl31_C_indecora) <- tissues_Car.Chl31_C_indecora

# fifth row of boxplot_positions gives the y coordinates for the tops of the whiskers
top_positions_Car.Chl31_C_indecora <- boxplot_positions_Car.Chl31_C_indecora[5,]


#add pairwise significance letter groups (compact letter display; CLD)
cbd_Car.Chl31_C_indecora <- as.data.frame(QsRutils::make_letter_assignments(wilcox_list_Grammica[["C_indecora__Car.Chl31"]])[["Letters"]])
colnames(cbd_Car.Chl31_C_indecora)[1] <- "Letter"
# turn rownames into first column for Tissue.code
setDT(cbd_Car.Chl31_C_indecora, keep.rownames = "Tissue.code")


# add a column y.position taken from top_positions based on mtaching up Tissue.code
# first reshape top_positions so that colnames are a column called Tissue.code
top_positions_Car.Chl31_C_indecora %>% gather(., Tissue.code, y.position) -> top_positions_Car.Chl31_C_indecora
# now join these positions to cbd
left_join(cbd_Car.Chl31_C_indecora, top_positions_Car.Chl31_C_indecora, by = "Tissue.code") -> cbd_Car.Chl31_C_indecora

# calculate how much to nudge
data_C_indecora_Car.Chl31 %>% group_by(Tissue.code) %>% dplyr::summarize(., max = max(logFW.norm)) -> max_Car.Chl31_C_indecora
cbd_Car.Chl31_C_indecora$nudged <- (max_Car.Chl31_C_indecora$max * 1.05)

# add CLDs to plot
Car.Chl31_C_indecora_boxplot + 
  geom_text(
    size    = 1.8,
    color = "black",
    data    = cbd_Car.Chl31_C_indecora,
    inherit.aes = T,
    mapping = aes(x = Tissue.code, y = nudged, label = Letter, vjust = 0)) -> Car.Chl31_C_indecora_boxplot


Car.Chl31_C_indecora_boxplot



#### COMBINED (faceted) car:chl plot: Grammica ONLY ####
wrap_elements(gridtext::richtext_grob('Carotenoid: Chorophyll', rot = 90, hjust = 0.5, vjust = 1, padding = unit(c(0, 0, 0, 0), "pt"), gp = gpar(fontsize = 6, fontface = 'bold'))) + 
  Car.Chl31_C_australis_boxplot + ggtitle('C. australis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_polygonorum_boxplot + ggtitle('C. polygonorum') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_sandwichiana_boxplot + ggtitle('C. sandwichiana') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_californica_boxplot + ggtitle('C. californica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_compacta_boxplot + ggtitle('C. compacta') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_cephalanthii_boxplot + ggtitle('C. cephalanthii') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_denticulata_boxplot + ggtitle('C. denticulata') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) + 
  Car.Chl31_C_tasmanica_boxplot + ggtitle('C. tasmanica') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Car.Chl31_C_costaricensis_boxplot + ggtitle('C. costaricensis') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  Car.Chl31_C_indecora_boxplot + ggtitle('C. indecora') + theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold.italic")) +
  plot_layout(nrow = 1, byrow = T) -> car.chl_boxplot_Grammica_new

car.chl_boxplot_Grammica_new 

pdf("../output/boxplots/car.chl_boxplot_Grammica_new.pdf", width=9,height=2) 
car.chl_boxplot_Grammica_new
dev.off()
