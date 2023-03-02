library(phytools)
# devtools::install_github('emmanuelparadis/ape')
library(ape) # cran version has some bugs in chronos()
# https://www.cynkra.com/blog/2021-03-16-gfortran-macos/ gcc with fortran installed by homebrew has problems
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

ctr <- chronos.control()
ctr$dual.iter.max <- 100
ctr$iter.max <- 1e7

# version with both rbcl and its for sandwichiana--not using
# # load nexus format phylogram tree
# # use ape::drop.tip to remove tips if needed
# input_tree <- read.nexus("Cuscuta_rbcL_26S.raxml.support")
# # input_tree <- ape::drop.tip(input_tree, "Cuscuta_sandwichiana_155_rbcL")
# outgroup = "Ipomoea_spp_AF146016_MG973745"
# input_tree_root_outgroup <- root(input_tree, outgroup, resolve.root=TRUE)
# calib_tree <- chronos(input_tree_root_outgroup, lambda = 1, model = "relaxed", quiet = FALSE,
#                       calibration = makeChronosCalib(input_tree_root_outgroup),
#                       control = ctr)
# write.tree(calib_tree, file = "Cuscuta_rbcL_26S.raxml.support_chrom.tre")


# 
# make a tree with only nuclear for neoxanthin analysis and with relabeled tips
input_tree <- read.nexus("../data/Cuscuta_rbcL_26S_2021_08_06-relabeled-taxa.nexus")
outgroup = "Ipomoea_spp_AF146016_MG973745"
input_tree_root_outgroup <- root(input_tree, outgroup, resolve.root=TRUE)
calib_tree <- chronos(input_tree_root_outgroup, lambda = 1, model = "correlated", quiet = FALSE,
                      calibration = makeChronosCalib(input_tree_root_outgroup),
                      control = ctr)
write.tree(calib_tree, file = "../output/Cuscuta_rbcL_26S.raxml.support_chrom_sandwichiaina_nr.tre")

# 
# make a tree with only chloroplast for neoxanthin analysis and with relabeled tips
input_tree <- read.nexus("../data/Cuscuta_rbcL_26S_2021_08_11-sandwichiana_cp.raxml.support.nexus")
outgroup = "Ipomoea_spp_AF146016_MG973745"
input_tree_root_outgroup <- root(input_tree, outgroup, resolve.root=TRUE)
calib_tree <- chronos(input_tree_root_outgroup, lambda = 1, model = "correlated", quiet = FALSE,
                      calibration = makeChronosCalib(input_tree_root_outgroup),
                      control = ctr)
write.tree(calib_tree, file = "../output/Cuscuta_rbcL_26S_2021_08_11-sandwichiana_cp.raxml.support_chrom_sandwichiaina_cp.tre")


# make a new version with sandwichiana omitted and with relabeled tips
input_tree <- read.nexus("../data/Cuscuta_rbcL_26S_2021_08_06-relabeled-taxa.nexus")
outgroup = "Ipomoea_spp_AF146016_MG973745"
input_tree_root_outgroup <- root(input_tree, outgroup, resolve.root=TRUE)
calib_tree <- chronos(input_tree_root_outgroup, lambda = 1, model = "correlated", quiet = FALSE,
                      calibration = makeChronosCalib(input_tree_root_outgroup),
                      control = ctr)
calib_tree <- ape::drop.tip(calib_tree, "Cuscuta_sandwichiana")

write.tree(calib_tree, file = "../output/Cuscuta_rbcL_26S.raxml.support_chrom_sandwichiaina_omitted.tre")

