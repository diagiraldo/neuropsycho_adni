####################################################
### MCI conversion prediction with domain scores ###
####################################################
### D. Giraldo, SEP2020
### Modified Jun 2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(randomForest)
library(ROCR)
library(reshape2)

source("scripts/plot_functions.R")
source("scripts/MCIpred_functions.R")

# Load data
load("processed_data/domain_scores_srb_firstvisit.RData")

# Filter data for EVALUATION
B <- S %>%
    filter(set == "EVAL")
rm(S)

# Include Raw test totals for baseline
load("processed_data/neuropsycho_seltests.RData")
sel_tests <- c("ADAS", "MMSE", "MOCA", "AVLT")
tot_tests <- unname(unlist(seltotals[sel_tests]))
B <- merge(B, dplyr::select(DF, c(1:5, tot_tests)), all.x = TRUE)
rm(DF, selitems, seltotals)
# Select totals and RAVLT-Immediate
tot <- tot_tests[c(1:4)]

# Vars and params for classification 
MCIlabs <- c("stableMCI", "converterMCI")
seltw <- seq(12,60,12)

###############################################
# Classification experiments in the manuscript 
###############################################
# List of classification experiments in manuscript
# 1. composite scores + edu + sex + age
# 2. tests scores + edu + sex + age
# 3. ADAS-Cog + edu + sex + age

expelist <- data.frame(fmlas = c(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c(tot, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("TOTAL13", "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ ")),
                       varsufix = c("dom_cov", "raw_cov", "adas_cov"),
                       expname = c("Domain scores", "Tests scores", "ADAS-Cog"),
                       stringsAsFactors = FALSE)

# Cross-validation params
trsplit <- 0.7

# data.frame with sample sizes 
sampsize <- samplesizes_train_test(B, seltw, MCIlabs, trsplit)

# RF params
ntr <- 200

# Cross-validation iterations
k_iter <- 1000

# Run cross-validation
results <- run_crossval_RF_expelist(B, MCIlabs, seltw, expelist, k_iter, rf.ntrees)

# Plot results
results <- results %>%
    mutate(time = factor(time, levels = seltw, labels = paste(seltw, "months")),
           input = factor(expname,
                          levels = expelist$expname))

perf_metric <- "auc"
ylab <- "AUC value"
pl <- MCIpred.boxplot(results, perf_metric, ylab)
pl
ggsave("plots/RFauc_boxplot_adas_revjun2021.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_adas_revjun2021.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Save results
save(results, sampsize, file = "results/results_MCIprogression_prediction_vJun2021.RData")

# -----
# Compare AUC 
seltw <- seq(12,60,12)
pval <- vector("numeric", length = length(seltw))
meandom <- vector("numeric", length = length(seltw))
meanother <- vector("numeric", length = length(seltw))
cohend <- vector("numeric", length = length(seltw))
for(t in 1:length(seltw)){
    tw <- seltw[t]
    tmpcomp <- filter(results, time == paste(tw, "months"))
    x <- tmpcomp$auc[tmpcomp$expname == "Domain scores"]
    y <- tmpcomp$auc[tmpcomp$expname == "Tests scores"]
    tt <- t.test(x,y, paired = TRUE)
    pval[t] <- tt$p.value
    meandom[t] <- round(mean(x), digits = 3)
    meanother[t] <- round(mean(y), digits = 3)
    cohend[t] <- tt$estimate/sd(x-y)
}
pval
cohend

