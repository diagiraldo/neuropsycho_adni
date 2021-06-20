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

B <- filter(S, set == "EVAL")

# Load Neuropsychiatric data
load("processed_data/neuropsychiatric_MCIinEval.RData")
B <- merge(B, PVIS, all = TRUE)

#############################
# Classification experiments
#############################

cv <- c("stableMCI", "converterMCI")
seltw <- seq(12,60,12)

# Cross-validation params
trsplit <- 0.7
k_iter <- 1000

# RF params
ntr <- 200

# Init data.frames to save results
sampsize <- data.frame(tw = seltw, 
                       stable = numeric(length(seltw)), 
                       converter = numeric(length(seltw)))

results <- data.frame(time = rep(seltw, each = k_iter), 
                      iter = rep(1:k_iter, length(seltw)), 
                      auc_dom_cov = numeric(length(seltw)*k_iter), 
                      auc_raw_cov = numeric(length(seltw)*k_iter),
                      auc_adas_cov = numeric(length(seltw)*k_iter),
                      auc_dom_cov_psy = numeric(length(seltw)*k_iter), 
                      auc_raw_cov_psy = numeric(length(seltw)*k_iter),
                      auc_adas_cov_psy = numeric(length(seltw)*k_iter))

# List of formulas for classification
# 1. composite scores + edu + sex + age
# 2. tests scores + edu + sex + age
# 3. ADAS-Cog + edu + sex + age
# 4. composite scores + edu + sex + age + neuropsychiatric scores
# 5. tests scores + edu + sex + age + neuropsychiatric scores
# 6. ADAS-Cog + edu + sex + age + neuropsychiatric scores

fmlas <- c(paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "),
           paste("GR", paste(c(tot, "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "),
           paste("GR", paste(c("TOTAL13", "PTEDUCAT", "PTGENDER", "AGE"), collapse = " + "), sep = " ~ "),
           paste("GR", paste(c(namesfac, "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ "),
           paste("GR", paste(c(tot, "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ "),
           paste("GR", paste(c("TOTAL13", "PTEDUCAT", "PTGENDER", "AGE", "GDTOTAL", "NPIQSCORE"), collapse = " + "), sep = " ~ "))

# Run Classification experiments
for (i in 1:length(seltw)){
    tw <- seltw[i]
    # Assign labels
    tmp <- mutate(B, 
                  GR = ifelse((new.is.stable | new.time.change > tw) & new.months.fu >= tw , "stableMCI", NA),
                  GR = ifelse((new.n.prog > 0 & new.time.change <= tw), "converterMCI", GR)) %>%
        filter(!is.na(GR)) %>%
        mutate(GR = factor(GR, labels = cv, levels = cv))
    # Save sample size
    sampsize[i, 2:3] <- table(tmp$GR)
    urn <- min(table(tmp$GR))
    urc <- names(table(tmp$GR))[which.min(table(tmp$GR))]
    tmp_impo <- vector("numeric", length = length(namesfac) + 3)
    # Run cross-validation for each time
    for (k in 1:k_iter){
        set.seed(1987 + tw + k)
        setTRAIN <- c(sample(tmp$RID[tmp$GR == urc], size = round(urn*trsplit)),
                      sample(tmp$RID[tmp$GR %in% setdiff(cv, urc)], size = round(urn*trsplit)))
        setTEST <- setdiff(tmp$RID, setTRAIN)
        idxTRAIN <- which(tmp$RID %in% setTRAIN)
        idxTEST <- which(tmp$RID %in% setTEST)
        # Train and test each RF
        for (ex in 1:length(fmlas)){
            rf <- randomForest(as.formula(fmlas[ex]), data = tmp, subset = idxTRAIN, ntree = ntr)
            pred_prob <- predict(rf, tmp[idxTEST,], type = "prob")[, "converterMCI"]
            pred <- prediction(pred_prob, tmp$GR[idxTEST], label.ordering = cv)
            results[results$time == tw & results$iter == k, 2 + ex] <- performance(pred, "auc")@y.values[[1]]
        }
    }
}

save(results, sampsize, file = "results/AUC_MCIprogression_predction_vJun2021.RData")

#############################
# Classification results
#############################

load("results/AUC_MCIprogression_predction_vJun2021.RData")

# Plot AUC boxplots per classifier and time window
m <- melt(results, id.vars = c("time", "iter")) %>%
    mutate(time = factor(time, levels = seltw, labels = paste(seltw, "months")),
           input = factor(variable, 
                          levels = c("auc_dom_cov", "auc_raw_cov", "auc_adas_cov", "auc_dom_cov_psy", "auc_raw_cov_psy", "auc_adas_cov_psy"), 
                          labels = c("Domain scores", "Tests scores", "ADAS-Cog", "Domain scores + NP", "Tests scores + NP", "ADAS-Cog + NP")))

# Original comparison: Domain scores Vs. Raw tests scores
pl <- ggplot(filter(m, variable %in% c("auc_dom_cov", "auc_raw_cov")), aes(time, value, colour = as.factor(input))) + 
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() + 
    labs(x = "Time period",y = "AUC value", colour = "RF classifiers\ntrained with: ") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
          legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_hline(yintercept = 0.5, colour="gray40", linetype = "dotted")
pl
ggsave("plots/RFauc_boxplot_revjun2021.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_revjun2021.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Comparison with ADAS-Cog: Domain scores Vs. Raw tests scores
pl <- ggplot(filter(m, variable %in% c("auc_dom_cov", "auc_raw_cov", "auc_adas_cov")), aes(time, value, colour = as.factor(input))) + 
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() + 
    labs(x = "Time period",y = "AUC value", colour = "RF classifiers\ntrained with: ") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
          legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_hline(yintercept = 0.5, colour="gray40", linetype = "dotted")
pl
ggsave("plots/RFauc_boxplot_adas_revjun2021.png", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_adas_revjun2021.eps", pl, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Including Neuropsychiatric symptoms
tmp_m <- m %>%
    mutate(labA = ifelse(grepl("dom", variable), "Domain scores", ifelse(grepl("raw", variable), "Tests scores", "ADAS-Cog")),
           labB = ifelse(grepl("psy", variable), "with psychiatric symptoms", "without psychiatric symptoms"),
           labA = factor(labA, levels = c("Domain scores", "Tests scores", "ADAS-Cog")))
pl <- ggplot(tmp_m, aes(time, value, colour = as.factor(labB))) + 
    geom_boxplot(outlier.size = 0.5) +
    facet_grid(labA ~ .) +
    theme_bw() + 
    labs(x = "Time period",y = "AUC value", colour = "Information in\nclassifier: ") +
    theme(plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
          legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_hline(yintercept = 0.5, colour="gray40", linetype = "dotted")
pl
ggsave("plots/RFauc_boxplot_psychiatric_revjun2021.png", pl, width = 16, height = 20, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/RFauc_boxplot_psychiatric_revjun2021.eps", pl, width = 16, height = 20, units = "cm", dpi = 300, bg = "transparent")

tmp <- filter(results, time == 60)
tt1 <- with(tmp, t.test(auc_fac, auc_tot))
d <- (tt1$estimate[1] - tt1$estimate[2])/sqrt((var(tmp$auc_fac) + var(tmp$auc_tot))/2)
