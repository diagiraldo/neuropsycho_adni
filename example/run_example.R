setwd("~/neuropsycho_adni")

source("scripts/example_functions.R")
load("example/tests_lists.RData")

sel_tests <- c("ADAS", "MMSE", "MOCA", "CLOCK", "TMT", "LM", "CATFL", "AVLT", "BNTMINT")
it <- unname(unlist(selitems[sel_tests]))

A <- read.csv("example/example_data.csv")
A <- reverse_scores(A, selitems)

# Standardization 
srbcoef <- read.csv("results/srb_parameters.csv")
A_srb <- calculate_srbz(A, it, srbcoef)
rm(srbcoef, selitems, seltotals)

# Calculate composite scores
load("results/cfa_estimates.RData")
namesfac <- colnames(LW$L)
S <- composite_scores(newdata = A_srb, weights = LW$W, centerit = meansub) 
S <- cbind(dplyr::select(A_srb, -all_of(it)), S)

# Assign MCI subgroup
grmeds <- read.csv("results/medoids_domainscores_MCIsubgroups_k4.csv")
MCIGR <- distance2meds(S, LW$Cz, grmeds, namesfac)
S <- merge(S, select(MCIGR, ID, GR))

# Predict progression
library(randomForest)
load("example/pretrained_RF_MCIprogression_prediction.RData")
PRED <- predict_MCIprogression(S, RFlist, namesfac)
S <- merge(S, select(PRED, ID, ends_with("prediction")))



