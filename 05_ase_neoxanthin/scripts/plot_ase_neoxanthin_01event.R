# library(devtools)
# devtools::install_github('cmt2/RevGadgets')
# devtools::install_github("YuLab-SMU/ggtree") 
# devtools::install_github("GuangchuangYu/treeio")
# install_github("lfabreti/convenience")
# install.packages("MCMCtrace")
library(RevGadgets)
library(coda)
library(ggplot2)
# library(convenience) 
library(rstudioapi)

# Getting the path of your current open file
# if not using rstudio, simply set your working directory to the scripts/ location of this script
# setwd(<location of scripts dir>)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
# print( getwd() )

#### summary of 4 runs combined ####
trace_quant <- RevGadgets::readTrace(path = '../output/../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj.log')
colnames(trace_quant[[1]]) 

RevGadgets::summarizeTrace(trace_quant, c('rates[1]', 'rates[2]', 'Posterior', 'Likelihood', 'Prior', 'is_irreverisible'))
# is_irreversible == 1 means it's irreversible
# so this actually means every samples was is_irreversible == 0
# (the posterior probability of the reversible model is 1)
# posterior of reversible model = 1

trace_quant_MCMC <- coda::as.mcmc(trace_quant[[1]])

coda::effectiveSize(trace_quant_MCMC)[c('rates[1]', 'rates[2]', 'Posterior', 'Likelihood', 'Prior', 'is_irreverisible')]

coda::traceplot(trace_quant_MCMC)

# ?processAncStates()
ase <- processAncStates('../output/output_ase_rjmcmc_01/ancestral_states_cuscuta_neoxanthin_rj.tree', state_labels = c("?" = "Missing", "0" = "Neoxanthin absent","1" = "Neoxanthin present"), labels_as_numbers = FALSE)

# ?plotAncStatesPIE()
# ?plotAncStatesMAP()


# plot
dir.create("../output/plots")
plotAncStatesPie(t = ase,
                 tip_labels_offset = 0.02, 
                 tip_labels_size = 2.25,
                 tip_pie_size = 1.0,
                 node_pie_size =1.5,
                 node_labels_size = 2.25,
                 tip_labels_italics = T,
                 node_labels_as = "state_posterior"
                 ) +
                 theme(legend.position = c(0.2,0.8),
                       legend.text=element_text(size=6)) -> ase_rj_pie


pdf("../output/plots/ASE_neoxanthin_rj_pie_01event_w_PP.pdf", width=6,height=3) 
ase_rj_pie
dev.off()






#### loop through 4 individual runs ####
runs <- 1:4
for (i in runs) {
  trace_quant[i] <- RevGadgets::readTrace(path = paste0('../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_', runs[i], '.log'))
  
  colnames(trace_quant[[1]]) 

  RevGadgets::summarizeTrace(trace_quant, c('rates[1]', 'rates[2]', 'Posterior', 'Likelihood', 'Prior', 'is_irreverisible'))
 

  trace_quant_MCMC <- coda::as.mcmc(trace_quant[[1]])
  
  coda::effectiveSize(trace_quant_MCMC)[c('rates[1]', 'rates[2]', 'Posterior', 'Likelihood', 'Prior', 'is_irreverisible')]
  
  coda::traceplot(trace_quant_MCMC)
}

# check parameter monitor for each run
print(RevGadgets::summarizeTrace(trace_quant[1], 'is_irreverisible'))
print(RevGadgets::summarizeTrace(trace_quant[2], 'is_irreverisible'))
print(RevGadgets::summarizeTrace(trace_quant[3], 'is_irreverisible'))
print(RevGadgets::summarizeTrace(trace_quant[4], 'is_irreverisible'))


# #### check convergence on 4 runs ####
# check_con <- convenience::checkConvergence(list_files = c("../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_1.log", "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_1.trees",
#                                  "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_2.log", "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_2.trees",
#                                  "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_3.log", "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_3.trees",
#                                  "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_4.log", "../output/output_ase_rjmcmc_01/cuscuta_neoxanthin_rj_run_4.trees") )
# 
# printTableSplits(check_con)
# 



