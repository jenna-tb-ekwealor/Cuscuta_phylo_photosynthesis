library(dplyr)
library(plyr)
library(ape)
library(data.table)
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

pigments <- read.csv("../../01_pigments/output/stat_results/pigments_species_summary.csv")

le <- pigments %>% filter(Pigment == "Lutein.epoxide")
neo <- pigments %>% filter(Pigment == "Neoxanthin")
car.chl <- pigments %>% filter(Pigment == "Car.Chl31")

pam <- read.csv("../../02_pam/output/stat_results/fluorescence_species_summary.csv")

fvfm <- pam %>% filter(Metric == "Fv.Fm") 	
phipsii <-  pam %>% filter(Metric == "φPSII")
npq <- pam %>% filter(Metric == "NPQ")

# omit C. sandwichiana
fvfm <- fvfm %>% filter(Species != "C_sandwichiana")
phipsii <- phipsii %>% filter(Species != "C_sandwichiana")
npq <- npq %>% filter(Species != "C_sandwichiana")



# FVFM LE ------------------------------------------------------------
#### for seedling ####
le_seedling <- le %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seedling)[colnames(le_seedling) == 'Mean'] <- 'le'

fvfm_seedling <- fvfm %>% filter(Tissue.edit == "sdlg") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_seedling)[colnames(fvfm_seedling) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_seedling <- full_join(le_seedling, fvfm_seedling, by = "Species")

# # make a duplicate of C_sandwichiana
# le_fvfm_seedling_C_sandwichiana <- le_fvfm_seedling %>% filter(Species == "C_sandwichiana") 
# le_fvfm_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_seedling_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# le_fvfm_seedling <- rbind(le_fvfm_seedling, le_fvfm_seedling_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_seedling %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_seedling

# omit Ipomoea 
le_fvfm_seedling <- le_fvfm_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_seedling <- le_fvfm_seedling %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_seedling <- na.omit(le_fvfm_seedling)

# transpose df 
le_fvfm_seedling <- data.table::transpose(le_fvfm_seedling)
names(le_fvfm_seedling) <- le_fvfm_seedling[1,] # move row 1 to col names
le_fvfm_seedling <- le_fvfm_seedling[-1,] # delete first row
# add rownames
rownames(le_fvfm_seedling) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_seedling, file = "../output/le_fvfm_seedling.nex", format = "continuous", datablock = TRUE)






#### for young stem ####
le_young <- le %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_young)[colnames(le_young) == 'Mean'] <- 'le'

fvfm_young <- fvfm %>% filter(Tissue.edit == "y") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_young)[colnames(fvfm_young) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_young <- full_join(le_young, fvfm_young, by = "Species")


# make a duplicate of C_sandwichiana
# le_fvfm_young_C_sandwichiana <- le_fvfm_young %>% filter(Species == "C_sandwichiana") 
# le_fvfm_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_young_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_fvfm_young <- rbind(le_fvfm_young, le_fvfm_young_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_young %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_young

# omit Ipomoea 
le_fvfm_young <- le_fvfm_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_young <- le_fvfm_young %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_young <- na.omit(le_fvfm_young)

# transpose df 
le_fvfm_young <- data.table::transpose(le_fvfm_young)
names(le_fvfm_young) <- le_fvfm_young[1,] # move row 1 to col names
le_fvfm_young <- le_fvfm_young[-1,] # delete first row
# add rownames
rownames(le_fvfm_young) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_young, file = "../output/le_fvfm_young.nex", format = "continuous", datablock = TRUE)





#### for old stem ####
le_old <- le %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_old)[colnames(le_old) == 'Mean'] <- 'le'

fvfm_old <- fvfm %>% filter(Tissue.edit == "o") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_old)[colnames(fvfm_old) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_old <- full_join(le_old, fvfm_old, by = "Species")


# make a duplicate of C_sandwichiana
# le_fvfm_old_C_sandwichiana <- le_fvfm_old %>% filter(Species == "C_sandwichiana") 
# le_fvfm_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_old_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_fvfm_old <- rbind(le_fvfm_old, le_fvfm_old_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_old %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_old

# omit Ipomoea 
le_fvfm_old <- le_fvfm_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_old <- le_fvfm_old %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_old <- na.omit(le_fvfm_old)

# transpose df 
le_fvfm_old <- data.table::transpose(le_fvfm_old)
names(le_fvfm_old) <- le_fvfm_old[1,] # move row 1 to col names
le_fvfm_old <- le_fvfm_old[-1,] # delete first row
# add rownames
rownames(le_fvfm_old) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_old, file = "../output/le_fvfm_old.nex", format = "continuous", datablock = TRUE)




#### for haustorium ####
le_haustorium <- le %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_haustorium)[colnames(le_haustorium) == 'Mean'] <- 'le'

fvfm_haustorium <- fvfm %>% filter(Tissue.edit == "h") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_haustorium)[colnames(fvfm_haustorium) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_haustorium <- full_join(le_haustorium, fvfm_haustorium, by = "Species")


# make a duplicate of C_sandwichiana
# le_fvfm_haustorium_C_sandwichiana <- le_fvfm_haustorium %>% filter(Species == "C_sandwichiana") 
# le_fvfm_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_haustorium_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_fvfm_haustorium <- rbind(le_fvfm_haustorium, le_fvfm_haustorium_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_haustorium %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_haustorium

# omit Ipomoea 
le_fvfm_haustorium <- le_fvfm_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_haustorium <- le_fvfm_haustorium %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_haustorium <- na.omit(le_fvfm_haustorium)

# transpose df 
le_fvfm_haustorium <- data.table::transpose(le_fvfm_haustorium)
names(le_fvfm_haustorium) <- le_fvfm_haustorium[1,] # move row 1 to col names
le_fvfm_haustorium <- le_fvfm_haustorium[-1,] # delete first row
# add rownames
rownames(le_fvfm_haustorium) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_haustorium, file = "../output/le_fvfm_haustorium.nex", format = "continuous", datablock = TRUE)




#### for flower ####
le_flower <- le %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_flower)[colnames(le_flower) == 'Mean'] <- 'le'

fvfm_flower <- fvfm %>% filter(Tissue.edit == "f") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_flower)[colnames(fvfm_flower) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_flower <- full_join(le_flower, fvfm_flower, by = "Species")


# make a duplicate of C_sandwichiana
# le_fvfm_flower_C_sandwichiana <- le_fvfm_flower %>% filter(Species == "C_sandwichiana") 
# le_fvfm_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_flower_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_fvfm_flower <- rbind(le_fvfm_flower, le_fvfm_flower_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_flower %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_flower

# omit Ipomoea 
le_fvfm_flower <- le_fvfm_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_flower <- le_fvfm_flower %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_flower <- na.omit(le_fvfm_flower)

# transpose df 
le_fvfm_flower <- data.table::transpose(le_fvfm_flower)
names(le_fvfm_flower) <- le_fvfm_flower[1,] # move row 1 to col names
le_fvfm_flower <- le_fvfm_flower[-1,] # delete first row
# add rownames
rownames(le_fvfm_flower) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_flower, file = "../output/le_fvfm_flower.nex", format = "continuous", datablock = TRUE)




#### for seed ####
le_seed <- le %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seed)[colnames(le_seed) == 'Mean'] <- 'le'

fvfm_seed <- fvfm %>% filter(Tissue.edit == "s") %>% dplyr::select(Species, mean)
# change "mean" to fvfm
colnames(fvfm_seed)[colnames(fvfm_seed) == 'mean'] <- 'fvfm'


# join le and fvfm 
le_fvfm_seed <- full_join(le_seed, fvfm_seed, by = "Species")


# make a duplicate of C_sandwichiana
# le_fvfm_seed_C_sandwichiana <- le_fvfm_seed %>% filter(Species == "C_sandwichiana") 
# le_fvfm_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_fvfm_seed_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_fvfm_seed <- rbind(le_fvfm_seed, le_fvfm_seed_C_sandwichiana)

# replace species names with how they are listed in the tree
le_fvfm_seed %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_fvfm_seed

# omit Ipomoea 
le_fvfm_seed <- le_fvfm_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_fvfm_seed <- le_fvfm_seed %>% dplyr::select(Species_tree, le, fvfm)

le_fvfm_seed <- na.omit(le_fvfm_seed)

# transpose df 
le_fvfm_seed <- data.table::transpose(le_fvfm_seed)
names(le_fvfm_seed) <- le_fvfm_seed[1,] # move row 1 to col names
le_fvfm_seed <- le_fvfm_seed[-1,] # delete first row
# add rownames
rownames(le_fvfm_seed) <- c("le", "fvfm") # did not need perhaps 

# write to nexus 
write.nexus.data(le_fvfm_seed, file = "../output/le_fvfm_seed.nex", format = "continuous", datablock = TRUE)







# PHIPSII LE ------------------------------------------------------------
##### for seedling #####
le_seedling <- le %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seedling)[colnames(le_seedling) == 'Mean'] <- 'le'

phipsii_seedling <- phipsii %>% filter(Tissue.edit == "sdlg") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_seedling)[colnames(phipsii_seedling) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_seedling <- full_join(le_seedling, phipsii_seedling, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_seedling_C_sandwichiana <- le_phipsii_seedling %>% filter(Species == "C_sandwichiana") 
# le_phipsii_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_seedling_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_seedling <- rbind(le_phipsii_seedling, le_phipsii_seedling_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_seedling %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_seedling

# omit Ipomoea 
le_phipsii_seedling <- le_phipsii_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_seedling <- le_phipsii_seedling %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_seedling <- na.omit(le_phipsii_seedling)

# transpose df 
le_phipsii_seedling <- data.table::transpose(le_phipsii_seedling)
names(le_phipsii_seedling) <- le_phipsii_seedling[1,] # move row 1 to col names
le_phipsii_seedling <- le_phipsii_seedling[-1,] # delete first row
# add rownames
rownames(le_phipsii_seedling) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_seedling, file = "../output/le_phipsii_seedling.nex", format = "continuous", datablock = TRUE)






##### for young stem #####
le_young <- le %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_young)[colnames(le_young) == 'Mean'] <- 'le'

phipsii_young <- phipsii %>% filter(Tissue.edit == "y") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_young)[colnames(phipsii_young) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_young <- full_join(le_young, phipsii_young, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_young_C_sandwichiana <- le_phipsii_young %>% filter(Species == "C_sandwichiana") 
# le_phipsii_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_young_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_young <- rbind(le_phipsii_young, le_phipsii_young_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_young %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_young

# omit Ipomoea 
le_phipsii_young <- le_phipsii_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_young <- le_phipsii_young %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_young <- na.omit(le_phipsii_young)

# transpose df 
le_phipsii_young <- data.table::transpose(le_phipsii_young)
names(le_phipsii_young) <- le_phipsii_young[1,] # move row 1 to col names
le_phipsii_young <- le_phipsii_young[-1,] # delete first row
# add rownames
rownames(le_phipsii_young) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_young, file = "../output/le_phipsii_young.nex", format = "continuous", datablock = TRUE)





##### for old stem #####
le_old <- le %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_old)[colnames(le_old) == 'Mean'] <- 'le'

phipsii_old <- phipsii %>% filter(Tissue.edit == "o") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_old)[colnames(phipsii_old) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_old <- full_join(le_old, phipsii_old, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_old_C_sandwichiana <- le_phipsii_old %>% filter(Species == "C_sandwichiana") 
# le_phipsii_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_old_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_old <- rbind(le_phipsii_old, le_phipsii_old_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_old %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_old

# omit Ipomoea 
le_phipsii_old <- le_phipsii_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_old <- le_phipsii_old %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_old <- na.omit(le_phipsii_old)

# transpose df 
le_phipsii_old <- data.table::transpose(le_phipsii_old)
names(le_phipsii_old) <- le_phipsii_old[1,] # move row 1 to col names
le_phipsii_old <- le_phipsii_old[-1,] # delete first row
# add rownames
rownames(le_phipsii_old) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_old, file = "../output/le_phipsii_old.nex", format = "continuous", datablock = TRUE)




##### for haustorium #####
le_haustorium <- le %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_haustorium)[colnames(le_haustorium) == 'Mean'] <- 'le'

phipsii_haustorium <- phipsii %>% filter(Tissue.edit == "h") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_haustorium)[colnames(phipsii_haustorium) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_haustorium <- full_join(le_haustorium, phipsii_haustorium, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_haustorium_C_sandwichiana <- le_phipsii_haustorium %>% filter(Species == "C_sandwichiana") 
# le_phipsii_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_haustorium_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_haustorium <- rbind(le_phipsii_haustorium, le_phipsii_haustorium_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_haustorium %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_haustorium

# omit Ipomoea 
le_phipsii_haustorium <- le_phipsii_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_haustorium <- le_phipsii_haustorium %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_haustorium <- na.omit(le_phipsii_haustorium)

# transpose df 
le_phipsii_haustorium <- data.table::transpose(le_phipsii_haustorium)
names(le_phipsii_haustorium) <- le_phipsii_haustorium[1,] # move row 1 to col names
le_phipsii_haustorium <- le_phipsii_haustorium[-1,] # delete first row
# add rownames
rownames(le_phipsii_haustorium) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_haustorium, file = "../output/le_phipsii_haustorium.nex", format = "continuous", datablock = TRUE)




##### for flower #####
le_flower <- le %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_flower)[colnames(le_flower) == 'Mean'] <- 'le'

phipsii_flower <- phipsii %>% filter(Tissue.edit == "f") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_flower)[colnames(phipsii_flower) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_flower <- full_join(le_flower, phipsii_flower, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_flower_C_sandwichiana <- le_phipsii_flower %>% filter(Species == "C_sandwichiana") 
# le_phipsii_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_flower_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_flower <- rbind(le_phipsii_flower, le_phipsii_flower_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_flower %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_flower

# omit Ipomoea 
le_phipsii_flower <- le_phipsii_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_flower <- le_phipsii_flower %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_flower <- na.omit(le_phipsii_flower)

# transpose df 
le_phipsii_flower <- data.table::transpose(le_phipsii_flower)
names(le_phipsii_flower) <- le_phipsii_flower[1,] # move row 1 to col names
le_phipsii_flower <- le_phipsii_flower[-1,] # delete first row
# add rownames
rownames(le_phipsii_flower) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_flower, file = "../output/le_phipsii_flower.nex", format = "continuous", datablock = TRUE)




##### for seed #####
le_seed <- le %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seed)[colnames(le_seed) == 'Mean'] <- 'le'

phipsii_seed <- phipsii %>% filter(Tissue.edit == "s") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_seed)[colnames(phipsii_seed) == 'mean'] <- 'phipsii'


# join le and phipsii 
le_phipsii_seed <- full_join(le_seed, phipsii_seed, by = "Species")


# make a duplicate of C_sandwichiana
# le_phipsii_seed_C_sandwichiana <- le_phipsii_seed %>% filter(Species == "C_sandwichiana") 
# le_phipsii_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_phipsii_seed_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_phipsii_seed <- rbind(le_phipsii_seed, le_phipsii_seed_C_sandwichiana)

# replace species names with how they are listed in the tree
le_phipsii_seed %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_phipsii_seed

# omit Ipomoea 
le_phipsii_seed <- le_phipsii_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_phipsii_seed <- le_phipsii_seed %>% dplyr::select(Species_tree, le, phipsii)

le_phipsii_seed <- na.omit(le_phipsii_seed)

# transpose df 
le_phipsii_seed <- data.table::transpose(le_phipsii_seed)
names(le_phipsii_seed) <- le_phipsii_seed[1,] # move row 1 to col names
le_phipsii_seed <- le_phipsii_seed[-1,] # delete first row
# add rownames
rownames(le_phipsii_seed) <- c("le", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(le_phipsii_seed, file = "../output/le_phipsii_seed.nex", format = "continuous", datablock = TRUE)






# CAR/CHL LE ------------------------------------------------------------
##### for seedling #####
le_seedling <- le %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seedling)[colnames(le_seedling) == 'Mean'] <- 'le'

car.chl_seedling <- car.chl %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_seedling)[colnames(car.chl_seedling) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_seedling <- full_join(le_seedling, car.chl_seedling, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_seedling_C_sandwichiana <- le_car.chl_seedling %>% filter(Species == "C_sandwichiana") 
# le_car.chl_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_seedling_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_seedling <- rbind(le_car.chl_seedling, le_car.chl_seedling_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_seedling %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_seedling

# omit Ipomoea 
le_car.chl_seedling <- le_car.chl_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_seedling <- le_car.chl_seedling %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_seedling <- na.omit(le_car.chl_seedling)

# transpose df 
le_car.chl_seedling <- data.table::transpose(le_car.chl_seedling)
names(le_car.chl_seedling) <- le_car.chl_seedling[1,] # move row 1 to col names
le_car.chl_seedling <- le_car.chl_seedling[-1,] # delete first row
# add rownames
rownames(le_car.chl_seedling) <- c("le", "car.chl") # did not need perhaps 

# write to nexus 
write.nexus.data(le_car.chl_seedling, file = "../output/le_car.chl_seedling.nex", format = "continuous", datablock = TRUE)






##### for young stem #####
le_young <- le %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_young)[colnames(le_young) == 'Mean'] <- 'le'

car.chl_young <- car.chl %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_young)[colnames(car.chl_young) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_young <- full_join(le_young, car.chl_young, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_young_C_sandwichiana <- le_car.chl_young %>% filter(Species == "C_sandwichiana") 
# le_car.chl_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_young_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_young <- rbind(le_car.chl_young, le_car.chl_young_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_young %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_young

# omit Ipomoea 
le_car.chl_young <- le_car.chl_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_young <- le_car.chl_young %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_young <- na.omit(le_car.chl_young)

# transpose df 
le_car.chl_young <- data.table::transpose(le_car.chl_young)
names(le_car.chl_young) <- le_car.chl_young[1,] # move row 1 to col names
le_car.chl_young <- le_car.chl_young[-1,] # delete first row
# add rownames
rownames(le_car.chl_young) <- c("le", "car.chl") # did not need perhaps 

# write to nexus 
write.nexus.data(le_car.chl_young, file = "../output/le_car.chl_young.nex", format = "continuous", datablock = TRUE)





##### for old stem #####
le_old <- le %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_old)[colnames(le_old) == 'Mean'] <- 'le'

car.chl_old <- car.chl %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_old)[colnames(car.chl_old) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_old <- full_join(le_old, car.chl_old, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_old_C_sandwichiana <- le_car.chl_old %>% filter(Species == "C_sandwichiana") 
# le_car.chl_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_old_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_old <- rbind(le_car.chl_old, le_car.chl_old_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_old %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_old

# omit Ipomoea 
le_car.chl_old <- le_car.chl_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_old <- le_car.chl_old %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_old <- na.omit(le_car.chl_old)

# transpose df 
le_car.chl_old <- data.table::transpose(le_car.chl_old)
names(le_car.chl_old) <- le_car.chl_old[1,] # move row 1 to col names
le_car.chl_old <- le_car.chl_old[-1,] # delete first row
# add rownames
rownames(le_car.chl_old) <- c("le", "car.chl") # did not need perhaps 

# write to nexus 
write.nexus.data(le_car.chl_old, file = "../output/le_car.chl_old.nex", format = "continuous", datablock = TRUE)




##### for haustorium #####
le_haustorium <- le %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_haustorium)[colnames(le_haustorium) == 'Mean'] <- 'le'

car.chl_haustorium <- car.chl %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_haustorium)[colnames(car.chl_haustorium) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_haustorium <- full_join(le_haustorium, car.chl_haustorium, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_haustorium_C_sandwichiana <- le_car.chl_haustorium %>% filter(Species == "C_sandwichiana") 
# le_car.chl_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_haustorium_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_haustorium <- rbind(le_car.chl_haustorium, le_car.chl_haustorium_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_haustorium %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_haustorium

# omit Ipomoea 
le_car.chl_haustorium <- le_car.chl_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_haustorium <- le_car.chl_haustorium %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_haustorium <- na.omit(le_car.chl_haustorium)

# transpose df 
le_car.chl_haustorium <- data.table::transpose(le_car.chl_haustorium)
names(le_car.chl_haustorium) <- le_car.chl_haustorium[1,] # move row 1 to col names
le_car.chl_haustorium <- le_car.chl_haustorium[-1,] # delete first row
# add rownames
rownames(le_car.chl_haustorium) <- c("le", "car.chl") # did not need perhaps 

# write to nexus 
write.nexus.data(le_car.chl_haustorium, file = "../output/le_car.chl_haustorium.nex", format = "continuous", datablock = TRUE)




##### for flower #####
le_flower <- le %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_flower)[colnames(le_flower) == 'Mean'] <- 'le'

car.chl_flower <- car.chl %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_flower)[colnames(car.chl_flower) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_flower <- full_join(le_flower, car.chl_flower, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_flower_C_sandwichiana <- le_car.chl_flower %>% filter(Species == "C_sandwichiana") 
# le_car.chl_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_flower_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_flower <- rbind(le_car.chl_flower, le_car.chl_flower_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_flower %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_flower

# omit Ipomoea 
le_car.chl_flower <- le_car.chl_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_flower <- le_car.chl_flower %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_flower <- na.omit(le_car.chl_flower)

# transpose df 
le_car.chl_flower <- data.table::transpose(le_car.chl_flower)
names(le_car.chl_flower) <- le_car.chl_flower[1,] # move row 1 to col names
le_car.chl_flower <- le_car.chl_flower[-1,] # delete first row
# add rownames
rownames(le_car.chl_flower) <- c("le", "car.chl") # did not need perhaps 

# write to nexus 
write.nexus.data(le_car.chl_flower, file = "../output/le_car.chl_flower.nex", format = "continuous", datablock = TRUE)




##### for seed #####
le_seed <- le %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seed)[colnames(le_seed) == 'Mean'] <- 'le'

car.chl_seed <- car.chl %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "mean" to car.chl
colnames(car.chl_seed)[colnames(car.chl_seed) == 'Mean'] <- 'car.chl'


# join le and car.chl 
le_car.chl_seed <- full_join(le_seed, car.chl_seed, by = "Species")


# make a duplicate of C_sandwichiana
# le_car.chl_seed_C_sandwichiana <- le_car.chl_seed %>% filter(Species == "C_sandwichiana") 
# le_car.chl_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_car.chl_seed_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_car.chl_seed <- rbind(le_car.chl_seed, le_car.chl_seed_C_sandwichiana)

# replace species names with how they are listed in the tree
le_car.chl_seed %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_car.chl_seed

# omit Ipomoea 
le_car.chl_seed <- le_car.chl_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_car.chl_seed <- le_car.chl_seed %>% dplyr::select(Species_tree, le, car.chl)

le_car.chl_seed <- na.omit(le_car.chl_seed)

# transpose df 
le_car.chl_seed <- data.table::transpose(le_car.chl_seed)
names(le_car.chl_seed) <- le_car.chl_seed[1,] # move row 1 to col names
le_car.chl_seed <- le_car.chl_seed[-1,] # delete first row
# add rownames
rownames(le_car.chl_seed) <- c("le", "car.chl") # did not need perhaps 


# write to nexus 
write.nexus.data(le_car.chl_seed, file = "../output/le_car.chl_seed.nex", format = "continuous", datablock = TRUE)



# PHIPSII NEOXANTHIN ------------------------------------------------------------
##### for seedling #####
neo_seedling <- neo %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(neo_seedling)[colnames(neo_seedling) == 'Mean'] <- 'neo'

phipsii_seedling <- phipsii %>% filter(Tissue.edit == "sdlg") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_seedling)[colnames(phipsii_seedling) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_seedling <- full_join(neo_seedling, phipsii_seedling, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_seedling_C_sandwichiana <- neo_phipsii_seedling %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_seedling_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_seedling <- rbind(neo_phipsii_seedling, neo_phipsii_seedling_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_seedling %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_seedling

# omit Ipomoea 
neo_phipsii_seedling <- neo_phipsii_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_seedling <- neo_phipsii_seedling %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_seedling <- na.omit(neo_phipsii_seedling)

# transpose df 
neo_phipsii_seedling <- data.table::transpose(neo_phipsii_seedling)
names(neo_phipsii_seedling) <- neo_phipsii_seedling[1,] # move row 1 to col names
neo_phipsii_seedling <- neo_phipsii_seedling[-1,] # delete first row
# add rownames
rownames(neo_phipsii_seedling) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_seedling, file = "../output/neo_phipsii_seedling.nex", format = "continuous", datablock = TRUE)






##### for young stem #####
neo_young <- neo %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "Mean" to neo
colnames(neo_young)[colnames(neo_young) == 'Mean'] <- 'neo'

phipsii_young <- phipsii %>% filter(Tissue.edit == "y") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_young)[colnames(phipsii_young) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_young <- full_join(neo_young, phipsii_young, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_young_C_sandwichiana <- neo_phipsii_young %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_young_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_young <- rbind(neo_phipsii_young, neo_phipsii_young_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_young %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_young

# omit Ipomoea 
neo_phipsii_young <- neo_phipsii_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_young <- neo_phipsii_young %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_young <- na.omit(neo_phipsii_young)

# transpose df 
neo_phipsii_young <- data.table::transpose(neo_phipsii_young)
names(neo_phipsii_young) <- neo_phipsii_young[1,] # move row 1 to col names
neo_phipsii_young <- neo_phipsii_young[-1,] # delete first row
# add rownames
rownames(neo_phipsii_young) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_young, file = "../output/neo_phipsii_young.nex", format = "continuous", datablock = TRUE)





##### for old stem #####
neo_old <- neo %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(neo_old)[colnames(neo_old) == 'Mean'] <- 'neo'

phipsii_old <- phipsii %>% filter(Tissue.edit == "o") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_old)[colnames(phipsii_old) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_old <- full_join(neo_old, phipsii_old, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_old_C_sandwichiana <- neo_phipsii_old %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_old_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_old <- rbind(neo_phipsii_old, neo_phipsii_old_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_old %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_old

# omit Ipomoea 
neo_phipsii_old <- neo_phipsii_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_old <- neo_phipsii_old %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_old <- na.omit(neo_phipsii_old)

# transpose df 
neo_phipsii_old <- data.table::transpose(neo_phipsii_old)
names(neo_phipsii_old) <- neo_phipsii_old[1,] # move row 1 to col names
neo_phipsii_old <- neo_phipsii_old[-1,] # delete first row
# add rownames
rownames(neo_phipsii_old) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_old, file = "../output/neo_phipsii_old.nex", format = "continuous", datablock = TRUE)




##### for haustorium #####
neo_haustorium <- neo %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(neo_haustorium)[colnames(neo_haustorium) == 'Mean'] <- 'neo'

phipsii_haustorium <- phipsii %>% filter(Tissue.edit == "h") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_haustorium)[colnames(phipsii_haustorium) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_haustorium <- full_join(neo_haustorium, phipsii_haustorium, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_haustorium_C_sandwichiana <- neo_phipsii_haustorium %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_haustorium_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_haustorium <- rbind(neo_phipsii_haustorium, neo_phipsii_haustorium_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_haustorium %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_haustorium

# omit Ipomoea 
neo_phipsii_haustorium <- neo_phipsii_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_haustorium <- neo_phipsii_haustorium %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_haustorium <- na.omit(neo_phipsii_haustorium)

# transpose df 
neo_phipsii_haustorium <- data.table::transpose(neo_phipsii_haustorium)
names(neo_phipsii_haustorium) <- neo_phipsii_haustorium[1,] # move row 1 to col names
neo_phipsii_haustorium <- neo_phipsii_haustorium[-1,] # delete first row
# add rownames
rownames(neo_phipsii_haustorium) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_haustorium, file = "../output/neo_phipsii_haustorium.nex", format = "continuous", datablock = TRUE)




##### for flower #####
neo_flower <- neo %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(neo_flower)[colnames(neo_flower) == 'Mean'] <- 'neo'

phipsii_flower <- phipsii %>% filter(Tissue.edit == "f") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_flower)[colnames(phipsii_flower) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_flower <- full_join(neo_flower, phipsii_flower, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_flower_C_sandwichiana <- neo_phipsii_flower %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_flower_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_flower <- rbind(neo_phipsii_flower, neo_phipsii_flower_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_flower %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_flower

# omit Ipomoea 
neo_phipsii_flower <- neo_phipsii_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_flower <- neo_phipsii_flower %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_flower <- na.omit(neo_phipsii_flower)

# transpose df 
neo_phipsii_flower <- data.table::transpose(neo_phipsii_flower)
names(neo_phipsii_flower) <- neo_phipsii_flower[1,] # move row 1 to col names
neo_phipsii_flower <- neo_phipsii_flower[-1,] # delete first row
# add rownames
rownames(neo_phipsii_flower) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_flower, file = "../output/neo_phipsii_flower.nex", format = "continuous", datablock = TRUE)




##### for seed #####
neo_seed <- neo %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(neo_seed)[colnames(neo_seed) == 'Mean'] <- 'neo'

phipsii_seed <- phipsii %>% filter(Tissue.edit == "s") %>% dplyr::select(Species, mean)
# change "mean" to phipsii
colnames(phipsii_seed)[colnames(phipsii_seed) == 'mean'] <- 'phipsii'


# join neo and phipsii 
neo_phipsii_seed <- full_join(neo_seed, phipsii_seed, by = "Species")


# make a duplicate of C_sandwichiana
# neo_phipsii_seed_C_sandwichiana <- neo_phipsii_seed %>% filter(Species == "C_sandwichiana") 
# neo_phipsii_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_seed_C_sandwichiana$Species)

# join second sandwichiana to full 
# neo_phipsii_seed <- rbind(neo_phipsii_seed, neo_phipsii_seed_C_sandwichiana)

# replace species names with how they are listed in the tree
neo_phipsii_seed %>%
  dplyr::mutate(Species_tree = case_when(
    Species == "C_australis"  ~  "Cuscuta_australis", 
    Species == "C_californica"  ~  "Cuscuta_californica",
    Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
    Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
    Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
    Species == "C_compacta"	~ "Cuscuta_compacta",
    Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
    Species == "C_denticulata"	~ "Cuscuta_denticulata",
    Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
    Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
    Species == "C_gracillima"	~ "Cuscuta_gracillima",
    Species == "C_indecora"	~ "Cuscuta_indecora",
    Species == "C_purpurata"	~ "Cuscuta_purpurata", 
    Species == "C_africana"	~ "Cuscuta_africana",
    Species == "C_epithymum"	~ "Cuscuta_epithymum",
    Species == "C_monogyna"	~ "Cuscuta_monogyna",
    Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_seed

# omit Ipomoea 
neo_phipsii_seed <- neo_phipsii_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

neo_phipsii_seed <- neo_phipsii_seed %>% dplyr::select(Species_tree, neo, phipsii)

neo_phipsii_seed <- na.omit(neo_phipsii_seed)

# transpose df 
neo_phipsii_seed <- data.table::transpose(neo_phipsii_seed)
names(neo_phipsii_seed) <- neo_phipsii_seed[1,] # move row 1 to col names
neo_phipsii_seed <- neo_phipsii_seed[-1,] # delete first row
# add rownames
rownames(neo_phipsii_seed) <- c("neo", "phipsii") # did not need perhaps 

# write to nexus 
write.nexus.data(neo_phipsii_seed, file = "../output/neo_phipsii_seed.nex", format = "continuous", datablock = TRUE)



