###############################################
### Standardized regression based subscores ###
###############################################
### D. Giraldo, Sept 2020
### modified Jun 2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(reshape2)

# Load data
load("processed_data/selsubscores_firstvisit.RData")

# Filter normative data
CN <- filter(A, set == "STD") %>%
    dplyr::select(it, PTEDUCAT, AGE)

# Init parameters for score standarization
coef <- data.frame(it = it, intercept = NA, PTEDUCAT = NA, AGE = NA, sigma = NA, r.squared = NA, fstat = NA, pval = NA)

# Linear regression with controlling variables fo all sub-scores
for (s in seq_along(it)){
    fmla <- as.formula(paste(it[s], paste(c("PTEDUCAT", "AGE"), collapse = " + "), sep = " ~ "))
    it_lm <- lm(fmla, data = CN)
    coef[s, 2:4] <- it_lm$coefficients
    coef[s, 5:7] <- c(summary(it_lm)$sigma, summary(it_lm)$r.squared, summary(it_lm)$fstatistic[1])
    coef$pval[s] <- pf(summary(it_lm)$fstatistic[1], summary(it_lm)$fstatistic[2], summary(it_lm)$fstatistic[3], lower.tail=F)
    rm(it_lm)
}

# Print table with parameters
i <- 1
tabstr <- paste(unlist(c(as.character(coef[i,1]), round(coef[i,c(2:5)], digits = 3))), collapse = " & ")
for (i in 2:nrow(coef)) {
    tabstr <- paste(tabstr, 
                    paste(unlist(c(as.character(coef[i,1]), round(coef[i,c(2:5)], digits = 3))), collapse = " & "),
                    sep = "\\\\\n")
}
cat(tabstr)
rm(tabstr)

# Save table with parameters
write.table(coef, file = "results/srb_parameters.csv", sep = ",", row.names = FALSE)

##########################
# Calculate SRB Z-scores
##########################
A_std <- data.frame(matrix(ncol = length(it), nrow = nrow(A)))
names(A_std) <- it

data4reg <- A %>% 
    mutate(intercept = 1) %>%
    dplyr::select(c("intercept", "PTEDUCAT", "AGE")) %>%
    as.matrix()

for (s in seq_along(it)){
    res <- A[, it[s]] - data4reg %*% t(as.matrix(coef[s, 2:4]))
    A_std[, it[s]] <- res/coef$sigma[s]
    attr(A_std[, it[s]], "dimnames") <- NULL
}

A_srb <- A
A_srb[, it] <- A_std

# Save data 
save(A_srb, it, file = "processed_data/srb_subscores_firstvisit.RData")
rm(list = ls())
