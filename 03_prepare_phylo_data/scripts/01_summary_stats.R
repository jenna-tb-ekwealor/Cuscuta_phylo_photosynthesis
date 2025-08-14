library(openxlsx)
library(tidyverse)
library(plyr)
library(rstatix)
library(gtools)
library(rstudioapi)
library(FSA)
library(data.table)

##### set working directory ##### 
# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )


##### load data #####
# skip the first row which has units info (check raw csv for units)
samples <- read.csv("../../01_pigments/data/pigment_organized_data_tab_AS_edit.csv", 
                    skip = 1, header = T, 
                    stringsAsFactors = FALSE, 
                    na.strings = c("NA","NaN","","#DIV/0!")) 


# DATA CLEANING -----------------------------------------------------------
##### sample cleanup #####

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


#### subset data ####
# select ng pigments, remove peak.11.9 column, omit d-Tocopherol because no standard and it's in area not ng
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
write.csv(ipo_le, file = "../output/ipomoea_le.csv", row.names = F)

# create table of just ipomoea beta Carotene for correlation analysis
ipo_b.Car <- data_long_calcs %>% filter(Species == "Ipomoea_nil") %>% filter(Pigment == "b.Carotene")
# save a copy of this dataset
write.csv(ipo_b.Car, file = "../output/ipomoea_b.Car.csv", row.names = F)


# STATISTICAL TESTS -------------------------------------------------------
##### overall patterns stats ##### 
# export only correlation pigments
correlation_pigs <- c("Neoxanthin", "Lutein.epoxide", "b.Carotene", "Lutein")

##### summary stats #####
# summary stats per Accession.No (summarizes those that have replicates), which is also per species: mean, n, and std dev  

summary_accession <- data_long_calcs %>%
  dplyr::filter(Pigment %in% correlation_pigs) %>%
  group_by(Subgenus, Species, Accession.No, Tissue.code, Pigment)  %>%
  dplyr::summarize(Mean = mean(FW.norm, na.rm = TRUE), n = sum(!is.na(FW.norm)), sd = sd(FW.norm, na.rm = TRUE))
# replace NaN with NA (for 0/0 ratios)
is.na(summary_accession$Mean) <- sapply(summary_accession$Mean, is.nan)

# write to csv
write_csv(summary_accession, file = "../output/pigments_species_summary_corrs.csv")
