####################################################
### MCI conversion prediction with domain scores ###
####################################################
### D. Giraldo, SEP2020

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(randomForest)
library(ROCR)
library(reshape2)

# Load data
load("processed_data/domain_scores_srb_firstvisit.RData")

# Load test totals for baseline 
load("processed_data/neuropsycho_seltests.RData")
sel_tests <- c("ADAS", "MMSE", "MOCA", "AVLT")
tot_tests <- unname(unlist(seltotals[sel_tests]))
S <- merge(S, dplyr::select(DF, c(1:5, tot_tests)), all.x = TRUE)
rm(DF, selitems, seltotals)

# Select totals and RAVLT-Immediate
tot <- tot_tests[c(1:4)]

B <- filter(S, set == "EVAL" & !is.na(CDRSB))

#############################
# Classification experiments
#############################

cv <- c("stableMCI", "converterMCI")
seltw <- seq(12,60,12)

trsplit <- 0.7
k_iter <- 1000

sampsize <- data.frame(tw = seltw, stable = numeric(length(seltw)), converter = numeric(length(seltw)))
results <- data.frame(time = rep(seltw, each = k_iter), iter = rep(1:k_iter, length(seltw)), 
                      auc_fac = numeric(length(seltw)*k_iter), auc_tot = numeric(length(seltw)*k_iter))
impo <- data.frame(time = rep(seltw, each = length(namesfac) + 3), 
                   namevar = character((length(namesfac) + 3)*length(seltw)), sumimpo = numeric((length(namesfac) + 3)*length(seltw)))

for (i in 1:length(seltw)){
    tw <- seltw[i]
    tmp <- mutate(B, GR = ifelse((new.is.stable | new.time.change > tw) & new.months.fu >= tw , "stableMCI", NA),
                  GR = ifelse((new.n.prog > 0 & new.time.change <= tw), "converterMCI", GR)) %>%
        filter(!is.na(GR)) %>%
        mutate(GR = factor(GR, labels = cv, levels = cv))
    sampsize[i, 2:3] <- table(tmp$GR)
    urn <- min(table(tmp$GR))
    urc <- names(table(tmp$GR))[which.min(table(tmp$GR))]
    tmp_impo <- vector("numeric", length = length(namesfac) + 3)
    for (k in 1:k_iter){
        set.seed(1987 + tw + k)
        setTRAIN <- c(sample(tmp$RID[tmp$GR == urc], size = round(urn*trsplit)),
                      sample(tmp$RID[tmp$GR %in% setdiff(cv, urc)], size = round(urn*trsplit)))
        setTEST <- setdiff(tmp$RID, setTRAIN)
        idxTRAIN <- which(tmp$RID %in% setTRAIN)
        idxTEST <- which(tmp$RID %in% setTEST)
        fmlafac <- as.formula(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "))
        rf <- randomForest(fmlafac, data = tmp, subset = idxTRAIN, importance = TRUE)
        tmp_impo <- tmp_impo + rf$importance[,3]
        pred_prob <- predict(rf, tmp[idxTEST,], type = "prob")[, "converterMCI"]
        pred <- prediction(pred_prob, tmp$GR[idxTEST], label.ordering = cv)
        aucfac = performance(pred, "auc")@y.values[[1]]
        fmlatot <- as.formula(paste("GR", paste(c(tot, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "))
        rf <- randomForest(fmlatot, data = tmp, subset = idxTRAIN)
        pred_prob <- predict(rf, tmp[idxTEST,], type = "prob")[, "converterMCI"]
        pred <- prediction(pred_prob, tmp$GR[idxTEST], label.ordering = cv)
        auctot = performance(pred, "auc")@y.values[[1]]
        results$auc_fac[results$time == tw & results$iter == k] <- aucfac
        results$auc_tot[results$time == tw & results$iter == k] <- auctot
    }
    impo$sumimpo[impo$time == tw] <- tmp_impo
}

impo$namevar <- rep(as.character(names(tmp_impo)), times = length(seltw))
impo <- dcast(impo, namevar ~ time, value.var = "sumimpo") %>%
    slice(match(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), namevar))
impo[, 2:6] <- sweep(impo[, 2:6], 2, colSums(impo[, 2:6]), '/')
impo[, 2:6] <- round(impo[, 2:6]*100, digits = 1)
apply(impo, 1, function(x) paste(x, collapse = "$ & $"))

save(results, sampsize, impo, file = "results/AUC_MCIprogression_predction.RData")

m <- melt(results, measure.vars = c("auc_fac", "auc_tot"))
m <- mutate(m, input = factor(variable, levels = c("auc_fac", "auc_tot"), labels = c("Domain Scores", "Other outcomes")),
            time = factor(time, levels = seltw, labels = paste(seltw, "months")))

pl <- ggplot(m, aes(time, value, colour = as.factor(input))) + geom_boxplot(outlier.size = 0.5) +
    theme_bw() + labs(x = "Time period",y = "AUC value", colour = "RF classifiers\ntrained with: ") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
          legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_hline(yintercept = 0.5, colour="gray40", linetype = "dotted")
pl
ggsave("plots/RFauc_boxplot.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")


tmp <- filter(results, time == 60)
tt1 <- with(tmp, t.test(auc_fac, auc_tot))
d <- (tt1$estimate[1] - tt1$estimate[2])/sqrt((var(tmp$auc_fac) + var(tmp$auc_tot))/2)
