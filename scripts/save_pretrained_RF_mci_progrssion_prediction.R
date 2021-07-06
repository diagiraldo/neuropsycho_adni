# Train Random forest to predict progression in new data
### D. Giraldo, Jun 2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(randomForest)

source("scripts/MCIpred_functions.R")

# Load data
load("processed_data/domain_scores_srb_firstvisit.RData")
# Filter data for EVALUATION
B <- S %>%
    filter(DIAGNOSIS == "MCI")

# 
expelist <- data.frame(fmlas = c(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ ")),
                       expname = c("Domain composite scores"),
                       stringsAsFactors = FALSE)

# data.frame with sample sizes 
MCIlabs <- c("stableMCI", "converterMCI")
seltw <- seq(12,60,12)
sampsize <- samplesizes_train_test(B, seltw, MCIlabs, trsplit = 1)

# RF params
ntr <- 200

# Train classifiers 
set.seed(2021)
RFlist = list()

for (i in 1:length(seltw)){
    tw <- seltw[i]
    tmpB <- assign_MCIlabs_filter(B, tw, MCIlabs)
    sampinfo <- filter(sampsize, time == tw)
    setTRAIN <- rids4train(tmpB, sampinfo, MCIlabs)
    idxTRAIN <- which(tmpB$RID %in% setTRAIN)
    RFlist[[i]] <- randomForest(as.formula(c(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "))),
                          data = tmpB, subset = idxTRAIN, ntree = ntr)
}
names(RFlist) <- paste(seltw, "months")

save(RFlist, sampsize, file = "example/pretrained_RF_MCIprogression_prediction.RData")


