# Functions to use in prediction of MCI progression
# Jun 2021

# Assign labels for classification
assign_MCIlabs_filter <- function(B, tw, MCIlabs){
    tmp <- mutate(B, 
                  GR = ifelse((new.is.stable | new.time.change > tw) & new.months.fu >= tw , "stableMCI", NA),
                  GR = ifelse((new.n.prog > 0 & new.time.change <= tw), "converterMCI", GR)) %>%
        filter(!is.na(GR)) %>%
        mutate(GR = factor(GR, labels = MCIlabs, levels = MCIlabs))
    return(tmp)
}

# Calculate sample sizes for training and testing
samplesizes_train_test <- function(B, seltw, MCIlabs, trsplit){
    sampsize <- data.frame(time = seltw, 
                           tot.sMCI = numeric(length(seltw)), 
                           tot.cMCI = numeric(length(seltw)), 
                           under.gr = character(length(seltw)),
                           under.n = numeric(length(seltw)),
                           train.n = numeric(length(seltw)),
                           test.sMCI = numeric(length(seltw)), 
                           test.cMCI = numeric(length(seltw)),
                           stringsAsFactors = FALSE)
    for(i in 1:length(seltw)){
        tmpB <- assign_MCIlabs_filter(B, seltw[i], MCIlabs)
        tab <- table(tmpB$GR)
        sampsize[i, 2:3] <- tab
        sampsize$under.gr[i] <- names(tab)[which.min(tab)]
        sampsize$under.n[i] <- min(tab)
        sampsize[i, 6:7] <- round(sampsize$under.n[i]*trsplit)
    }
    sampsize <- mutate(sampsize,
                       test.sMCI = tot.sMCI - train.n,
                       test.cMCI = tot.cMCI - train.n)
    return(sampsize)
}

# Sample RIDS for training
rids4train <- function(tmpB, sampinfo, MCIlabs){
    tr.under <- sample(tmpB$RID[tmpB$GR == sampinfo$under.gr], size = sampinfo$train.n)
    tr.other <- sample(tmpB$RID[tmpB$GR %in% setdiff(MCIlabs, sampinfo$under.gr)], size = sampinfo$train.n)
    return(c(tr.under, tr.other))
}

# Run one iteration of the cross-validation for the experiments in a list
run_train_test_RF_expelist <- function(tmpB, MCIlabs, setTRAIN, setTEST, expelist, rf.ntrees){
    idxTRAIN <- which(tmpB$RID %in% setTRAIN)
    idxTEST <- which(tmpB$RID %in% setTEST)
    results <- data.frame(expname = expelist$expname,
                          auc = numeric(nrow(expelist)),
                          sens = numeric(nrow(expelist)),
                          spec = numeric(nrow(expelist)),
                          acc = numeric(nrow(expelist)),
                          stringsAsFactors = FALSE)
    for(ex in 1:nrow(expelist)){
        trainedRF <- randomForest(as.formula(expelist$fmlas[ex]), data = tmpB, subset = idxTRAIN, ntree = rf.ntrees)
        confmat <- table(tmpB$GR[idxTEST], predict(trainedRF, tmpB[idxTEST,]))
        results$sens[ex] <- confmat[MCIlabs[2], MCIlabs[2]]/rowSums(confmat)[MCIlabs[2]]
        results$spec[ex] <- confmat[MCIlabs[1], MCIlabs[1]]/rowSums(confmat)[MCIlabs[1]]
        results$acc[ex] <- sum(diag(confmat))/sum(confmat)
        pred_prob <- predict(trainedRF, tmpB[idxTEST,], type = "prob")[, "converterMCI"]
        pred <- prediction(pred_prob, tmpB$GR[idxTEST], label.ordering = MCIlabs)
        results$auc[ex] <- performance(pred, "auc")@y.values[[1]]
    }
    return(results)
}

# Run cross-validation with data in B for the experiments in a list
run_crossval_RF_expelist <- function(B, MCIlabs, seltw, expelist, k_iter, rf.ntrees){
    # Init data.frame to save results
    results <- data.frame(time = rep(seltw, each = k_iter*nrow(expelist)), 
                          iter = rep(1:k_iter, each = nrow(expelist), times = length(seltw)),
                          expname = rep(expelist$expname, times = k_iter*length(seltw)),
                          auc = numeric(length(seltw)*k_iter*nrow(expelist)),
                          sens = numeric(length(seltw)*k_iter*nrow(expelist)),
                          spec = numeric(length(seltw)*k_iter*nrow(expelist)),
                          acc = numeric(length(seltw)*k_iter*nrow(expelist)))
    # Run Classification experiments
    for (i in 1:length(seltw)){
        tw <- seltw[i]
        tmpB <- assign_MCIlabs_filter(B, tw, MCIlabs)
        sampinfo <- filter(sampsize, time == tw)
        cat("Starting ", k_iter, "iterations of cross-validation for prediction within ", tw, "months\n")
        # Run cross-validation iterations
        for(k in 1:k_iter){
            set.seed(1987 + tw + k)
            setTRAIN <- rids4train(tmpB, sampinfo, MCIlabs)
            setTEST <- setdiff(tmpB$RID, setTRAIN)
            # Run train and test each RF 
            tmpres <- run_train_test_RF_expelist(tmpB, MCIlabs, setTRAIN, setTEST, expelist, rf.ntrees = ntr)
            results[results$time == tw & results$iter == k, 3:7] <- tmpres
        }
    }
    return(results)
}
