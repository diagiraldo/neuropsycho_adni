####################################################
### MCI conversion prediction compatison with other results ###
####################################################
### D. Giraldo,  Jun 2021

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

# Load Neuropsychiatric data
load("processed_data/neuropsychiatric_MCIinEval.RData")
B <- merge(B, PVIS, all = TRUE)
rm(PVIS)

# Vars and params for classification 
MCIlabs <- c("stableMCI", "converterMCI")
seltw <- seq(12,60,12)

#######################################################
# Classification experiments with neuropsychiatric data 
#######################################################
# List of classification experiments with neuropsychiatric scores
# 4. composite scores + edu + sex + age + neuropsychiatric scores
# 5. tests scores + edu + sex + age + neuropsychiatric scores
# 6. ADAS-Cog + edu + sex + age + neuropsychiatric scores

expelist <- data.frame(fmlas = c(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c(tot, "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("TOTAL13", "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ ")),
                       varsufix = c("dom_cov_np", "raw_cov_np", "adas_cov_np"),
                       expname = c("Domain scores + NP", "Tests scores + NP", "ADAS-Cog + NP"),
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
NPresults <- run_crossval_RF_expelist(B, MCIlabs, seltw, expelist, k_iter, rf.ntrees = ntr)

# Plot results
NPresults <- NPresults %>%
    mutate(time = factor(time, levels = seltw, labels = paste(seltw, "months")),
           input = factor(expname,
                          levels = expelist$expname))

perf_metric <- "auc"
ylab <- "AUC value"
pl <- MCIpred.boxplot(NPresults, perf_metric, ylab)
pl
ggsave("plots/RFauc_boxplot_withNP_revjun2021.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_withNP_revjun2021.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Save results
save(NPresults, file = "results/results_withpsychiatric_MCIprogression_prediction_vJun2021.RData")

# Compare with first results (without neuropsychiatric info)
load("results/results_MCIprogression_prediction_vJun2021.RData")
NPresults <- rbind(NPresults, results)
input0lev <- levels(results$input)
rm(results)

NPresults <- NPresults %>%
    mutate(labA = gsub(" \\+ NP", "", expname),
           labB = ifelse(grepl("NP", expname), "with psychiatric symptoms", "without psychiatric symptoms"),
           labA = factor(labA, levels = input0lev))

perf_metric <- "auc"
ylab <- "AUC value"
pl <- MCIpred.boxplot(filter(NPresults, labA == "Domain scores"), perf_metric, ylab, "labB") 
pl
ggsave("plots/RFauc_boxplot_psychiatric_revjun2021.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_psychiatric_revjun2021.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Compare AUC 
pval <- vector("numeric", length = length(seltw))
meanx <- vector("numeric", length = length(seltw))
meany <- vector("numeric", length = length(seltw))
cohend <- vector("numeric", length = length(seltw))
for(t in 1:length(seltw)){
    tw <- seltw[t]
    tmpcomp <- filter(NPresults, labA == "Domain scores" & time == paste(tw, "months"))
    x <- tmpcomp$auc[tmpcomp$labB == "with psychiatric symptoms"]
    y <- tmpcomp$auc[tmpcomp$labB == "without psychiatric symptoms"]
    tt <- t.test(x,y, paired = TRUE)
    pval[t] <- tt$p.value
    meanx[t] <- round(mean(x), digits = 3)
    meany[t] <- round(mean(y), digits = 3)
    cohend[t] <- tt$estimate/sd(x-y)
}
pval
meandif
cohend


rm(NPresults)

##################################################################
# Classification experiments to compare with the state-of-the-art
##################################################################

# Load data to compare with other composite scores and classifiers in state-of-the-art
load("processed_data/other_composites_features.RData")
B <- merge(B, BCOMP, all.x = TRUE)
rm(BCOMP)

# List of classification experiments with neuropsychiatric scores
# 1. Domain composite scores
# 2. Domain composite scores + CDRSB + FAQTOTAL
# 3. Adas.Tree
# 4. comp.Huang
# 5. CC1
# 6. CC2 +
# 7. CFC1 +
# 8. CFC2 + 
# 9. NP Features selected in Pereira2018

expelist <- data.frame(fmlas = c(paste("GR", paste(c(namesfac), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c(namesfac, "CDRSB", "FAQTOTAL"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("ADAS.tree"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("comp.Huang"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("CC1"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("CC2"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("CFC1"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c("CFC2"), collapse = " + "), sep = " ~ "),
                                 paste("GR", paste(c(selfeaturesPer), collapse = " + "), sep = " ~ ")),
                       expname = c("Domain scores", "Domain scores + CDR + FAQ",  "ADAS Tree\n(Llano, 2011)", 
                                   "Composite\n(Huang, 2015)", "CC1", "CC2", "CFC1", "CFC2", 
                                   "Selected features\n(Pereira, 2018)"),
                       stringsAsFactors = FALSE)

# Filter B to have all required info
selidx <- (complete.cases(B[, c(selfeaturesPer, "CDRSB")]) & B$Forget.index != Inf)
B2 <- B[selidx,]

# Cross-validation params
trsplit <- 0.7

# data.frame with sample sizes 
sampsize <- samplesizes_train_test(B2, seltw, MCIlabs, trsplit)

# RF params
ntr <- 100

# Cross-validation iterations
k_iter <- 200

# Run cross-validation
SOTAresults <- run_crossval_RF_expelist(B2, MCIlabs, seltw, expelist, k_iter, rf.ntrees = ntr)

# Plot results
SOTAresults <- SOTAresults %>%
    mutate(time = factor(time, levels = seltw, labels = paste(seltw, "months")),
           input = factor(expname,
                          levels = expelist$expname))
levels(SOTAresults$input)[2] <- "Domain scores,\nCDR and FAQ"

perf_metric <- "auc"
ylab <- "AUC value"
pl <- MCIpred.boxplot(SOTAresults, "auc", ylab) +
    theme(plot.margin = unit(c(0,0.5,0,0.25), "cm"))
    
ggsave("plots/RFauc_boxplot_comparesota_revjun2021.png", pl, width = 18, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_comparesota_revjun2021.eps", pl, width = 18, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Save results
save(SOTAresults, file = "results/results_sota_MCIprogression_prediction_vJun2021.RData")

