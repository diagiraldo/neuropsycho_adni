##########################################
### Confirmatory Factor Analysis (CFA) ###
### With SRB subscores ###################
##########################################
### D. Giraldo, Sept 2020
### Modified Jun 2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(reshape2)
library(lavaan)

source("scripts/plot_functions.R")
source("scripts/cfa_functions.R")

# Load data
load("processed_data/srb_subscores_firstvisit.RData")

# Filter data for CFA
B <- A_srb %>%
    filter(set == "FA")

# Specify factor model

famodel1 <- '
MEMORY =~ Q1SCORE + Q4SCORE + MOCADLREC + RAVLT.IMMED + AVTOT6 + AVTOTB + AVDEL30MIN + AVDELTOT + LIMMTOTAL + LDELTOTAL + MMRECALL
LANGUAGE =~ Q5SCORE + MOCANAM + BMNOCUE + CATANIMSC
EXECUTIVE =~ Q13SCORE + TRAASCOR + TRABSCOR + MOCASERIAL + TRAILS
VISUOSPATIAL =~ CLOCKSCOR + COPYSCOR + MOCACLOCK + Q3SCORE + CUBE + MMDRAW
ORIENTATION =~ Q7SCORE + MMORITIME + MMORISPACE + MOCAORI
ATTENTION =~ Q9SCORE + Q10SCORE + Q11SCORE + Q12SCORE + Q2SCORE
'
itmod <- its_in_model(famodel1)

# Run CFA
meansub <- colMeans(B[, itmod]) 
fitcfa1 <- cfa(famodel1, data = B[, itmod], sample.cov = cov(B[, itmod]), estimator = "ULS")
fitMeasures(fitcfa1)
s <- lavInspect(fitcfa1, what = "est")
m <- modindices(fitcfa1, sort = TRUE)
mloads <- m[m$op == "=~",]

# Extract CFA estimates: factor loadings and factor covariance, calculate weights
LW <- loadings_weights_fcov(fitcfa1)
namesfac <- colnames(LW$L)

# Calculate domain composite scores
# Same result as lavPredict(fitcfa1, newdata = A_srb[,itmod], method = "Bartlett")
S <- composite_scores(newdata = A_srb, weights = LW$W, centerit = meansub) 
S <- cbind(dplyr::select(A_srb, c(1:6, DIAGNOSIS:set)), S)

# Compare domain scores between all CN and MCI participants
grcomp <- compare_domain_scores(S, vars = namesfac, grvar = "DIAGNOSIS")

# Save CFA fit and weights
save(fitcfa1, file = "results/cfa_modelfit.RData")
save(meansub, LW, file = "results/cfa_estimates.RData")
# Save Weights in csv
save_weights_csv(LW$W, "results/weights_for_domainscores.csv")
# Save domain composite scores
save(S, namesfac, file = "processed_data/domain_scores_srb_firstvisit.RData")

# Plot Weights
pl <- weight.plot(LW$W, nam = c("Domain score", "Sub-score", "Weight"),
                  valrange = c(0,0.5), withvalues = TRUE)
ggsave("plots/DS_weights_cfa_revjun2021.png", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/DS_weights_cfa_revjun2021.eps", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")

# Plot Loadings
pl <- weight.plot(t(LW$L), nam = c("Factor", "Sub-score", "Loading"),
                  valrange = c(0,3.1), withvalues = TRUE)
pl
ggsave("plots/Factorloadings_cfa_revjun2021.png", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/Factorloadings_cfa_revjun2021.eps", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")

# Save Weights in csv
L <- LW$L
L <- L[rev(rownames(L)), rev(colnames(L))]
save_weights_csv(t(L), "results/factor_loadings.csv")

# Save comparisons between groups
write.table(grcomp, file = "results/comparison_domainscores_CN_MCI.csv", sep = ",", row.names = FALSE)

