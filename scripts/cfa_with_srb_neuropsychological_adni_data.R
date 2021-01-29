##########################################
### Confirmatory Factor Analysis (CFA) ###
### With SRB subscores ###################
##########################################
### D. Giraldo, Sept 2020

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(reshape2)
library(lavaan)

source("scripts/plot_functions.R")

# Load data
load("processed_data/srb_subscores_firstvisit.RData")

# Filter data for CFA
B <- A_srb %>%
    filter(set == "FA")

# Examine Correlations
COR <- B %>%
    dplyr::select(it, -MMNAM, - BMCUED) %>%
    cor(.)

corr.plot(COR, cut = FALSE, corrange = c(-0.1,1))

# Kaiser-Meyer-Olkin measure of sample adequacy and Bartlett test of sphericity
library(psych)
KMO(COR)
cortest.bartlett(COR, n = nrow(B))

# Confirmatory factor analysis (CFA)
model1 <- '
MEMORY =~ Q1SCORE + Q4SCORE + MOCADLREC + RAVLT.IMMED + AVTOT6 + AVTOTB + AVDEL30MIN + AVDELTOT + LIMMTOTAL + LDELTOTAL + MMRECALL
LANGUAGE =~ Q5SCORE + MOCANAM + BMNOCUE + CATANIMSC
EXECUTIVE =~ Q13SCORE + TRAASCOR + TRABSCOR + MOCASERIAL + TRAILS
VISUOSPATIAL =~ CLOCKSCOR + COPYSCOR + MOCACLOCK + Q3SCORE + CUBE + MMDRAW
ORIENTATION =~ Q7SCORE + MMORITIME + MMORISPACE + MOCAORI
ATTENTION =~ Q9SCORE + Q10SCORE + Q11SCORE + Q12SCORE + Q2SCORE
'
tmp <- dplyr::select(B, it)
fitcfa <- cfa(model1, data = tmp, sample.cov = cov(tmp), estimator = "ULS")
fitMeasures(fitcfa)

save(fitcfa, file = "results/cfa_parameters.RData")

# Examine structure and factor loadings
L <- inspect(fitcfa, what = "est")$lambda
pl <- weight.plot(L, nam = c("Sub-score", "Factor", "Factor Loading"), 
                  valrange = c(0,2.25), flip = TRUE, withvalues = TRUE)
pl

# Extract matrix of weights for domain score calculation
U2 <- lavInspect(fitcfa, what = "est")$theta
Uinv <- diag(1/diag(U2))
tmp <- t(L)%*%(Uinv)
W <- solve(tmp%*%L)%*%tmp
colnames(W) <- rownames(L)
W <- W[rev(1:nrow(W)),]
W <- W[,rev(1:ncol(W))]

pl <- weight.plot(W, nam = c("Domain score", "Sub-score", "Weight"), 
                  valrange = c(0,0.5), withvalues = TRUE)
pl
ggsave("plots/DS_weights_cfa.png", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/DS_weights_cfa.eps", pl, width = 10, height = 16, units = "cm", dpi = 300, bg = "transparent")

W1 <- as.data.frame(t(W)) %>%
    mutate(subscore = colnames(W)) %>%
    select(c(subscore, rownames(W)))
write.table(W1, file = "results/weights_for_domainscores.csv", sep = ",", row.names = FALSE)

# Calculate domain scores
S <- lavPredict(fitcfa, newdata = A_srb, method = "Bartlett") %>%
    as.data.frame() %>%
    cbind(dplyr::select(A_srb, c(1:6, DIAGNOSIS:set))) %>%
    mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("Dementia", "MCI", "CN")))

namesfac <- colnames(inspect(fitcfa)$lambda)
pl <- scores.boxplots(filter(S, set == "FA"), namesfac, nam = c("Domain", "Score", "Diagnostic Group: "), flip = TRUE)
pl

save(S, namesfac, file = "processed_data/domain_scores_srb_firstvisit.RData")

# Compare dds between CN and MCI
library(rcompanion)
n_tests <- 6

i <- 6
fmla <- as.formula(paste(namesfac[i], "DIAGNOSIS", sep = " ~ "))
tmp <- wilcox.test(fmla, data = S)
tmp$data.name
tmp$p.value*n_tests # Bonferroni correction

x <- S[, which(names(S) == namesfac[i])]
g <- as.character(S$DIAGNOSIS)
wilcoxonR(x, g, ci = TRUE)
