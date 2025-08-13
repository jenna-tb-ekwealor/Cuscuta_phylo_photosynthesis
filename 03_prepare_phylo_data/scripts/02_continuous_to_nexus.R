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

pigments <- read.csv("../output/pigments_species_summary_corrs.csv")

le <- pigments %>% filter(Pigment == "Lutein.epoxide")
neo <- pigments %>% filter(Pigment == "Neoxanthin")
b.Car <- pigments %>% filter(Pigment == "b.Carotene")

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






# # tot CAR LE ------------------------------------------------------------
# ##### for seedling #####
# le_seedling <- le %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_seedling)[colnames(le_seedling) == 'Mean'] <- 'le'
# 
# Tot.Car_seedling <- Tot.Car %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_seedling)[colnames(Tot.Car_seedling) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_seedling <- full_join(le_seedling, Tot.Car_seedling, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_seedling_C_sandwichiana <- le_Tot.Car_seedling %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_seedling_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_seedling <- rbind(le_Tot.Car_seedling, le_Tot.Car_seedling_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_seedling %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_seedling
# 
# # omit Ipomoea 
# le_Tot.Car_seedling <- le_Tot.Car_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_seedling <- le_Tot.Car_seedling %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_seedling <- na.omit(le_Tot.Car_seedling)
# 
# # transpose df 
# le_Tot.Car_seedling <- data.table::transpose(le_Tot.Car_seedling)
# names(le_Tot.Car_seedling) <- le_Tot.Car_seedling[1,] # move row 1 to col names
# le_Tot.Car_seedling <- le_Tot.Car_seedling[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_seedling) <- c("le", "Tot.Car") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_seedling, file = "../output/le_Tot.Car_seedling.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# 
# 
# ##### for young stem #####
# le_young <- le %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_young)[colnames(le_young) == 'Mean'] <- 'le'
# 
# Tot.Car_young <- Tot.Car %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_young)[colnames(Tot.Car_young) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_young <- full_join(le_young, Tot.Car_young, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_young_C_sandwichiana <- le_Tot.Car_young %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_young_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_young <- rbind(le_Tot.Car_young, le_Tot.Car_young_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_young %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_young
# 
# # omit Ipomoea 
# le_Tot.Car_young <- le_Tot.Car_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_young <- le_Tot.Car_young %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_young <- na.omit(le_Tot.Car_young)
# 
# # transpose df 
# le_Tot.Car_young <- data.table::transpose(le_Tot.Car_young)
# names(le_Tot.Car_young) <- le_Tot.Car_young[1,] # move row 1 to col names
# le_Tot.Car_young <- le_Tot.Car_young[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_young) <- c("le", "Tot.Car") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_young, file = "../output/le_Tot.Car_young.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# 
# ##### for old stem #####
# le_old <- le %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_old)[colnames(le_old) == 'Mean'] <- 'le'
# 
# Tot.Car_old <- Tot.Car %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_old)[colnames(Tot.Car_old) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_old <- full_join(le_old, Tot.Car_old, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_old_C_sandwichiana <- le_Tot.Car_old %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_old_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_old <- rbind(le_Tot.Car_old, le_Tot.Car_old_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_old %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_old
# 
# # omit Ipomoea 
# le_Tot.Car_old <- le_Tot.Car_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_old <- le_Tot.Car_old %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_old <- na.omit(le_Tot.Car_old)
# 
# # transpose df 
# le_Tot.Car_old <- data.table::transpose(le_Tot.Car_old)
# names(le_Tot.Car_old) <- le_Tot.Car_old[1,] # move row 1 to col names
# le_Tot.Car_old <- le_Tot.Car_old[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_old) <- c("le", "Tot.Car") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_old, file = "../output/le_Tot.Car_old.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for haustorium #####
# le_haustorium <- le %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_haustorium)[colnames(le_haustorium) == 'Mean'] <- 'le'
# 
# Tot.Car_haustorium <- Tot.Car %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_haustorium)[colnames(Tot.Car_haustorium) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_haustorium <- full_join(le_haustorium, Tot.Car_haustorium, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_haustorium_C_sandwichiana <- le_Tot.Car_haustorium %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_haustorium_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_haustorium <- rbind(le_Tot.Car_haustorium, le_Tot.Car_haustorium_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_haustorium %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_haustorium
# 
# # omit Ipomoea 
# le_Tot.Car_haustorium <- le_Tot.Car_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_haustorium <- le_Tot.Car_haustorium %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_haustorium <- na.omit(le_Tot.Car_haustorium)
# 
# # transpose df 
# le_Tot.Car_haustorium <- data.table::transpose(le_Tot.Car_haustorium)
# names(le_Tot.Car_haustorium) <- le_Tot.Car_haustorium[1,] # move row 1 to col names
# le_Tot.Car_haustorium <- le_Tot.Car_haustorium[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_haustorium) <- c("le", "Tot.Car") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_haustorium, file = "../output/le_Tot.Car_haustorium.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for flower #####
# le_flower <- le %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_flower)[colnames(le_flower) == 'Mean'] <- 'le'
# 
# Tot.Car_flower <- Tot.Car %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_flower)[colnames(Tot.Car_flower) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_flower <- full_join(le_flower, Tot.Car_flower, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_flower_C_sandwichiana <- le_Tot.Car_flower %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_flower_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_flower <- rbind(le_Tot.Car_flower, le_Tot.Car_flower_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_flower %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_flower
# 
# # omit Ipomoea 
# le_Tot.Car_flower <- le_Tot.Car_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_flower <- le_Tot.Car_flower %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_flower <- na.omit(le_Tot.Car_flower)
# 
# # transpose df 
# le_Tot.Car_flower <- data.table::transpose(le_Tot.Car_flower)
# names(le_Tot.Car_flower) <- le_Tot.Car_flower[1,] # move row 1 to col names
# le_Tot.Car_flower <- le_Tot.Car_flower[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_flower) <- c("le", "Tot.Car") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_flower, file = "../output/le_Tot.Car_flower.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for seed #####
# le_seed <- le %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(le_seed)[colnames(le_seed) == 'Mean'] <- 'le'
# 
# Tot.Car_seed <- Tot.Car %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# # change "mean" to Tot.Car
# colnames(Tot.Car_seed)[colnames(Tot.Car_seed) == 'Mean'] <- 'Tot.Car'
# 
# 
# # join le and Tot.Car 
# le_Tot.Car_seed <- full_join(le_seed, Tot.Car_seed, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # le_Tot.Car_seed_C_sandwichiana <- le_Tot.Car_seed %>% filter(Species == "C_sandwichiana") 
# # le_Tot.Car_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_Tot.Car_seed_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # le_Tot.Car_seed <- rbind(le_Tot.Car_seed, le_Tot.Car_seed_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# le_Tot.Car_seed %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_Tot.Car_seed
# 
# # omit Ipomoea 
# le_Tot.Car_seed <- le_Tot.Car_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# le_Tot.Car_seed <- le_Tot.Car_seed %>% dplyr::select(Species_tree, le, Tot.Car)
# 
# le_Tot.Car_seed <- na.omit(le_Tot.Car_seed)
# 
# # transpose df 
# le_Tot.Car_seed <- data.table::transpose(le_Tot.Car_seed)
# names(le_Tot.Car_seed) <- le_Tot.Car_seed[1,] # move row 1 to col names
# le_Tot.Car_seed <- le_Tot.Car_seed[-1,] # delete first row
# # add rownames
# rownames(le_Tot.Car_seed) <- c("le", "Tot.Car") # did not need perhaps 
# 
# 
# # write to nexus 
# write.nexus.data(le_Tot.Car_seed, file = "../output/le_Tot.Car_seed.nex", format = "continuous", datablock = TRUE)
# 



# BETA CAROTENE LE ------------------------------------------------------------
##### for seedling #####
le_seedling <- le %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seedling)[colnames(le_seedling) == 'Mean'] <- 'le'

b.Car_seedling <- b.Car %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_seedling)[colnames(b.Car_seedling) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_seedling <- full_join(le_seedling, b.Car_seedling, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_seedling_C_sandwichiana <- le_b.Car_seedling %>% filter(Species == "C_sandwichiana") 
# le_b.Car_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_seedling_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_seedling <- rbind(le_b.Car_seedling, le_b.Car_seedling_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_seedling %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_seedling

# omit Ipomoea 
le_b.Car_seedling <- le_b.Car_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_seedling <- le_b.Car_seedling %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_seedling <- na.omit(le_b.Car_seedling)

# transpose df 
le_b.Car_seedling <- data.table::transpose(le_b.Car_seedling)
names(le_b.Car_seedling) <- le_b.Car_seedling[1,] # move row 1 to col names
le_b.Car_seedling <- le_b.Car_seedling[-1,] # delete first row
# add rownames
rownames(le_b.Car_seedling) <- c("le", "b.Car") # did not need perhaps 

# write to nexus 
write.nexus.data(le_b.Car_seedling, file = "../output/le_b.Car_seedling.nex", format = "continuous", datablock = TRUE)






##### for young stem #####
le_young <- le %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_young)[colnames(le_young) == 'Mean'] <- 'le'

b.Car_young <- b.Car %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_young)[colnames(b.Car_young) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_young <- full_join(le_young, b.Car_young, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_young_C_sandwichiana <- le_b.Car_young %>% filter(Species == "C_sandwichiana") 
# le_b.Car_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_young_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_young <- rbind(le_b.Car_young, le_b.Car_young_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_young %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_young

# omit Ipomoea 
le_b.Car_young <- le_b.Car_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_young <- le_b.Car_young %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_young <- na.omit(le_b.Car_young)

# transpose df 
le_b.Car_young <- data.table::transpose(le_b.Car_young)
names(le_b.Car_young) <- le_b.Car_young[1,] # move row 1 to col names
le_b.Car_young <- le_b.Car_young[-1,] # delete first row
# add rownames
rownames(le_b.Car_young) <- c("le", "b.Car") # did not need perhaps 

# write to nexus 
write.nexus.data(le_b.Car_young, file = "../output/le_b.Car_young.nex", format = "continuous", datablock = TRUE)





##### for old stem #####
le_old <- le %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_old)[colnames(le_old) == 'Mean'] <- 'le'

b.Car_old <- b.Car %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_old)[colnames(b.Car_old) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_old <- full_join(le_old, b.Car_old, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_old_C_sandwichiana <- le_b.Car_old %>% filter(Species == "C_sandwichiana") 
# le_b.Car_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_old_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_old <- rbind(le_b.Car_old, le_b.Car_old_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_old %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_old

# omit Ipomoea 
le_b.Car_old <- le_b.Car_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_old <- le_b.Car_old %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_old <- na.omit(le_b.Car_old)

# transpose df 
le_b.Car_old <- data.table::transpose(le_b.Car_old)
names(le_b.Car_old) <- le_b.Car_old[1,] # move row 1 to col names
le_b.Car_old <- le_b.Car_old[-1,] # delete first row
# add rownames
rownames(le_b.Car_old) <- c("le", "b.Car") # did not need perhaps 

# write to nexus 
write.nexus.data(le_b.Car_old, file = "../output/le_b.Car_old.nex", format = "continuous", datablock = TRUE)




##### for haustorium #####
le_haustorium <- le %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_haustorium)[colnames(le_haustorium) == 'Mean'] <- 'le'

b.Car_haustorium <- b.Car %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_haustorium)[colnames(b.Car_haustorium) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_haustorium <- full_join(le_haustorium, b.Car_haustorium, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_haustorium_C_sandwichiana <- le_b.Car_haustorium %>% filter(Species == "C_sandwichiana") 
# le_b.Car_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_haustorium_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_haustorium <- rbind(le_b.Car_haustorium, le_b.Car_haustorium_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_haustorium %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_haustorium

# omit Ipomoea 
le_b.Car_haustorium <- le_b.Car_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_haustorium <- le_b.Car_haustorium %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_haustorium <- na.omit(le_b.Car_haustorium)

# transpose df 
le_b.Car_haustorium <- data.table::transpose(le_b.Car_haustorium)
names(le_b.Car_haustorium) <- le_b.Car_haustorium[1,] # move row 1 to col names
le_b.Car_haustorium <- le_b.Car_haustorium[-1,] # delete first row
# add rownames
rownames(le_b.Car_haustorium) <- c("le", "b.Car") # did not need perhaps 

# write to nexus 
write.nexus.data(le_b.Car_haustorium, file = "../output/le_b.Car_haustorium.nex", format = "continuous", datablock = TRUE)




##### for flower #####
le_flower <- le %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_flower)[colnames(le_flower) == 'Mean'] <- 'le'

b.Car_flower <- b.Car %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_flower)[colnames(b.Car_flower) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_flower <- full_join(le_flower, b.Car_flower, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_flower_C_sandwichiana <- le_b.Car_flower %>% filter(Species == "C_sandwichiana") 
# le_b.Car_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_flower_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_flower <- rbind(le_b.Car_flower, le_b.Car_flower_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_flower %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_flower

# omit Ipomoea 
le_b.Car_flower <- le_b.Car_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_flower <- le_b.Car_flower %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_flower <- na.omit(le_b.Car_flower)

# transpose df 
le_b.Car_flower <- data.table::transpose(le_b.Car_flower)
names(le_b.Car_flower) <- le_b.Car_flower[1,] # move row 1 to col names
le_b.Car_flower <- le_b.Car_flower[-1,] # delete first row
# add rownames
rownames(le_b.Car_flower) <- c("le", "b.Car") # did not need perhaps 

# write to nexus 
write.nexus.data(le_b.Car_flower, file = "../output/le_b.Car_flower.nex", format = "continuous", datablock = TRUE)




##### for seed #####
le_seed <- le %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "Mean" to le
colnames(le_seed)[colnames(le_seed) == 'Mean'] <- 'le'

b.Car_seed <- b.Car %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# change "mean" to b.Car
colnames(b.Car_seed)[colnames(b.Car_seed) == 'Mean'] <- 'b.Car'


# join le and b.Car 
le_b.Car_seed <- full_join(le_seed, b.Car_seed, by = "Species")


# make a duplicate of C_sandwichiana
# le_b.Car_seed_C_sandwichiana <- le_b.Car_seed %>% filter(Species == "C_sandwichiana") 
# le_b.Car_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", le_b.Car_seed_C_sandwichiana$Species)

# join second sandwichiana to full 
# le_b.Car_seed <- rbind(le_b.Car_seed, le_b.Car_seed_C_sandwichiana)

# replace species names with how they are listed in the tree
le_b.Car_seed %>%
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
    Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> le_b.Car_seed

# omit Ipomoea 
le_b.Car_seed <- le_b.Car_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")

le_b.Car_seed <- le_b.Car_seed %>% dplyr::select(Species_tree, le, b.Car)

le_b.Car_seed <- na.omit(le_b.Car_seed)

# transpose df 
le_b.Car_seed <- data.table::transpose(le_b.Car_seed)
names(le_b.Car_seed) <- le_b.Car_seed[1,] # move row 1 to col names
le_b.Car_seed <- le_b.Car_seed[-1,] # delete first row
# add rownames
rownames(le_b.Car_seed) <- c("le", "b.Car") # did not need perhaps 


# write to nexus 
write.nexus.data(le_b.Car_seed, file = "../output/le_b.Car_seed.nex", format = "continuous", datablock = TRUE)


# 
# 
# # PHIPSII NEOXANTHIN ------------------------------------------------------------
# ##### for seedling #####
# neo_seedling <- neo %>% filter(Tissue.code == "sdlg") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(neo_seedling)[colnames(neo_seedling) == 'Mean'] <- 'neo'
# 
# phipsii_seedling <- phipsii %>% filter(Tissue.edit == "sdlg") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_seedling)[colnames(phipsii_seedling) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_seedling <- full_join(neo_seedling, phipsii_seedling, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_seedling_C_sandwichiana <- neo_phipsii_seedling %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_seedling_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_seedling_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_seedling <- rbind(neo_phipsii_seedling, neo_phipsii_seedling_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_seedling %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_seedling
# 
# # omit Ipomoea 
# neo_phipsii_seedling <- neo_phipsii_seedling %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_seedling <- neo_phipsii_seedling %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_seedling <- na.omit(neo_phipsii_seedling)
# 
# # transpose df 
# neo_phipsii_seedling <- data.table::transpose(neo_phipsii_seedling)
# names(neo_phipsii_seedling) <- neo_phipsii_seedling[1,] # move row 1 to col names
# neo_phipsii_seedling <- neo_phipsii_seedling[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_seedling) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_seedling, file = "../output/neo_phipsii_seedling.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# 
# 
# ##### for young stem #####
# neo_young <- neo %>% filter(Tissue.code == "y") %>% dplyr::select(Species, Mean)
# # change "Mean" to neo
# colnames(neo_young)[colnames(neo_young) == 'Mean'] <- 'neo'
# 
# phipsii_young <- phipsii %>% filter(Tissue.edit == "y") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_young)[colnames(phipsii_young) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_young <- full_join(neo_young, phipsii_young, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_young_C_sandwichiana <- neo_phipsii_young %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_young_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_young_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_young <- rbind(neo_phipsii_young, neo_phipsii_young_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_young %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_young
# 
# # omit Ipomoea 
# neo_phipsii_young <- neo_phipsii_young %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_young <- neo_phipsii_young %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_young <- na.omit(neo_phipsii_young)
# 
# # transpose df 
# neo_phipsii_young <- data.table::transpose(neo_phipsii_young)
# names(neo_phipsii_young) <- neo_phipsii_young[1,] # move row 1 to col names
# neo_phipsii_young <- neo_phipsii_young[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_young) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_young, file = "../output/neo_phipsii_young.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# 
# ##### for old stem #####
# neo_old <- neo %>% filter(Tissue.code == "o") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(neo_old)[colnames(neo_old) == 'Mean'] <- 'neo'
# 
# phipsii_old <- phipsii %>% filter(Tissue.edit == "o") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_old)[colnames(phipsii_old) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_old <- full_join(neo_old, phipsii_old, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_old_C_sandwichiana <- neo_phipsii_old %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_old_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_old_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_old <- rbind(neo_phipsii_old, neo_phipsii_old_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_old %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_old
# 
# # omit Ipomoea 
# neo_phipsii_old <- neo_phipsii_old %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_old <- neo_phipsii_old %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_old <- na.omit(neo_phipsii_old)
# 
# # transpose df 
# neo_phipsii_old <- data.table::transpose(neo_phipsii_old)
# names(neo_phipsii_old) <- neo_phipsii_old[1,] # move row 1 to col names
# neo_phipsii_old <- neo_phipsii_old[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_old) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_old, file = "../output/neo_phipsii_old.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for haustorium #####
# neo_haustorium <- neo %>% filter(Tissue.code == "h") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(neo_haustorium)[colnames(neo_haustorium) == 'Mean'] <- 'neo'
# 
# phipsii_haustorium <- phipsii %>% filter(Tissue.edit == "h") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_haustorium)[colnames(phipsii_haustorium) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_haustorium <- full_join(neo_haustorium, phipsii_haustorium, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_haustorium_C_sandwichiana <- neo_phipsii_haustorium %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_haustorium_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_haustorium_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_haustorium <- rbind(neo_phipsii_haustorium, neo_phipsii_haustorium_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_haustorium %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_haustorium
# 
# # omit Ipomoea 
# neo_phipsii_haustorium <- neo_phipsii_haustorium %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_haustorium <- neo_phipsii_haustorium %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_haustorium <- na.omit(neo_phipsii_haustorium)
# 
# # transpose df 
# neo_phipsii_haustorium <- data.table::transpose(neo_phipsii_haustorium)
# names(neo_phipsii_haustorium) <- neo_phipsii_haustorium[1,] # move row 1 to col names
# neo_phipsii_haustorium <- neo_phipsii_haustorium[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_haustorium) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_haustorium, file = "../output/neo_phipsii_haustorium.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for flower #####
# neo_flower <- neo %>% filter(Tissue.code == "f") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(neo_flower)[colnames(neo_flower) == 'Mean'] <- 'neo'
# 
# phipsii_flower <- phipsii %>% filter(Tissue.edit == "f") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_flower)[colnames(phipsii_flower) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_flower <- full_join(neo_flower, phipsii_flower, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_flower_C_sandwichiana <- neo_phipsii_flower %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_flower_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_flower_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_flower <- rbind(neo_phipsii_flower, neo_phipsii_flower_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_flower %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_flower
# 
# # omit Ipomoea 
# neo_phipsii_flower <- neo_phipsii_flower %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_flower <- neo_phipsii_flower %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_flower <- na.omit(neo_phipsii_flower)
# 
# # transpose df 
# neo_phipsii_flower <- data.table::transpose(neo_phipsii_flower)
# names(neo_phipsii_flower) <- neo_phipsii_flower[1,] # move row 1 to col names
# neo_phipsii_flower <- neo_phipsii_flower[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_flower) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_flower, file = "../output/neo_phipsii_flower.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
# 
# ##### for seed #####
# neo_seed <- neo %>% filter(Tissue.code == "s") %>% dplyr::select(Species, Mean)
# # change "Mean" to le
# colnames(neo_seed)[colnames(neo_seed) == 'Mean'] <- 'neo'
# 
# phipsii_seed <- phipsii %>% filter(Tissue.edit == "s") %>% dplyr::select(Species, mean)
# # change "mean" to phipsii
# colnames(phipsii_seed)[colnames(phipsii_seed) == 'mean'] <- 'phipsii'
# 
# 
# # join neo and phipsii 
# neo_phipsii_seed <- full_join(neo_seed, phipsii_seed, by = "Species")
# 
# 
# # make a duplicate of C_sandwichiana
# # neo_phipsii_seed_C_sandwichiana <- neo_phipsii_seed %>% filter(Species == "C_sandwichiana") 
# # neo_phipsii_seed_C_sandwichiana$Species <- gsub("C_sandwichiana", "C_sandwichiana_2", neo_phipsii_seed_C_sandwichiana$Species)
# 
# # join second sandwichiana to full 
# # neo_phipsii_seed <- rbind(neo_phipsii_seed, neo_phipsii_seed_C_sandwichiana)
# 
# # replace species names with how they are listed in the tree
# neo_phipsii_seed %>%
#   dplyr::mutate(Species_tree = case_when(
#     Species == "C_australis"  ~  "Cuscuta_australis", 
#     Species == "C_californica"  ~  "Cuscuta_californica",
#     Species == "C_sandwichiana" ~  "Cuscuta_sandwichiana",
#     Species == "C_sandwichiana_2"	~ "Cuscuta_sandwichiana",
#     Species == "C_polygonorum"	~ "Cuscuta_polygonorum",
#     Species == "C_compacta"	~ "Cuscuta_compacta",
#     Species == "C_cephalanthii"	~ "Cuscuta_cephalanthi",
#     Species == "C_denticulata"	~ "Cuscuta_denticulata",
#     Species == "C_tasmanica"	~ "Cuscuta_tasmanica",
#     Species == "C_costaricensis"	~ "Cuscuta_costaricensis", 
#     Species == "C_gracillima"	~ "Cuscuta_gracillima",
#     Species == "C_indecora"	~ "Cuscuta_indecora",
#     Species == "C_purpurata"	~ "Cuscuta_purpurata", 
#     Species == "C_africana"	~ "Cuscuta_africana",
#     Species == "C_epithymum"	~ "Cuscuta_epithymum",
#     Species == "C_monogyna"	~ "Cuscuta_monogyna",
#     Species == "C_lupuliformis" ~ "Cuscuta_lupuliformis",
#     Species == "Ipomoea_nil" ~ "Ipomoea_spp_AF146016_MG973745"), .before = Species) -> neo_phipsii_seed
# 
# # omit Ipomoea 
# neo_phipsii_seed <- neo_phipsii_seed %>% filter(Species_tree != "Ipomoea_spp_AF146016_MG973745")
# 
# neo_phipsii_seed <- neo_phipsii_seed %>% dplyr::select(Species_tree, neo, phipsii)
# 
# neo_phipsii_seed <- na.omit(neo_phipsii_seed)
# 
# # transpose df 
# neo_phipsii_seed <- data.table::transpose(neo_phipsii_seed)
# names(neo_phipsii_seed) <- neo_phipsii_seed[1,] # move row 1 to col names
# neo_phipsii_seed <- neo_phipsii_seed[-1,] # delete first row
# # add rownames
# rownames(neo_phipsii_seed) <- c("neo", "phipsii") # did not need perhaps 
# 
# # write to nexus 
# write.nexus.data(neo_phipsii_seed, file = "../output/neo_phipsii_seed.nex", format = "continuous", datablock = TRUE)
# 
# 
# 
