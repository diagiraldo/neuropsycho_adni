###############################################
### Standardized regression based subscores ###
###############################################
### D. Giraldo, Sept 2020

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(reshape2)

# Load data
load("processed_data/selsubscores_firstvisit.RData")

# Learn parameters for score standarization
coef <- data.frame(it = it, intercept = NA, PTEDUCAT = NA, AGE = NA)
rsq <- data.frame(it = it, sigma = NA, r.squared = NA, adj.r.squared = NA)

CN <- filter(A, set == "STD") %>%
    dplyr::select(it, PTEDUCAT, AGE)

# Linear regression with controlling variables
for (s in seq_along(it)){
    fmla <- as.formula(paste(it[s], paste(c("PTEDUCAT", "AGE"), collapse = " + "), sep = " ~ "))
    tmp <- lm(fmla, data = CN)
    coef[s, 2:4] <- summary(tmp)$coefficients[,1]
    rsq[s, 2:4] <- c(summary(tmp)$sigma, summary(tmp)$r.squared, summary(tmp)$adj.r.squared)
    rm(tmp)
}

# Calculate SRB Z-scores
A_std <- data.frame(matrix(ncol = length(it), nrow = nrow(A)))
names(A_std) <- it

tmp <- A %>% 
    mutate(intercept = 1) %>%
    dplyr::select(c("intercept", "PTEDUCAT", "AGE")) %>%
    as.matrix()

for (s in seq_along(it)){
    res <- A[, it[s]] - tmp %*% t(as.matrix(coef[s, 2:4]))
    A_std[, it[s]] <- res/rsq$sigma[s]
    attr(A_std[, it[s]], "dimnames") <- NULL
}

A_srb <- A
A_srb[, it] <- A_std

# Save data 
save(coef, rsq, file = "results/srb_parameters_cn.RData")
save(A_srb, it, file = "processed_data/srb_subscores_firstvisit.RData")

rm(list = ls())

