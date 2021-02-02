#####################################
### Clustering with domain scores ###
#####################################
### D. Giraldo, Sept 2020

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(reshape2)

# Packages for cluster analysis
library(NbClust)
library(cluster)

# Packages for Survival analysis
library(survival)
library(survminer)

source("scripts/plot_functions.R")

# Load data
load("processed_data/domain_scores_srb_firstvisit.RData")

# Filter MCI patients in Evaluation set for clustering 
B <- S %>%
    filter(set == "EVAL") %>%
    mutate(event = ifelse(new.n.prog > 0, 1, 0), 
           time = ifelse(new.n.prog > 0, new.time.change, new.months.fu)) 

######################################
# Calculate distance between subjects
######################################

# Load factor covariance matrix
library(lavaan)
load("results/cfa_parameters.RData")
factorcov <- inspect(fitcfa, what = "est")$psi
rm(fitcfa)

# Get weighted distance
tmp <- dplyr::select(B, all_of(namesfac))
n <- nrow(tmp) 
dwei <- matrix(data = NA, nrow = n, ncol = n)
for (i in 1:n){
    for (j in 1:i){
        x <- as.matrix(tmp[i,] - tmp[j,])
        dwei[i,j] <- x %*% factorcov %*% t(x)
    }
}
d <- as.dist(sqrt(dwei))

##############################
# Explore number of clusters
##############################

NK <- NbClust(data = tmp, diss = d, distance = NULL, min.nc = 2, max.nc = 10, method = "ward.D2", index = "all")
# --> 4 clusters

#Clusted stability index
n_it <- 1000
bfrac <- 0.8

cbi <- list()
for (k_clust in 2:9){
    pam.res <- pam(d, k = k_clust, diss = TRUE, pamonce = 5)
    B <- mutate(B, GR = factor(pam.res$clustering))
    cum_jac <- rep(0, k_clust)
    for (bit in 1:n_it){
        set.seed(1987 + bit)
        samp <- sample(B$RID, size = round(nrow(B)*bfrac))
        tmp_idx <- which(B$RID %in% samp)
        tmp_B <- B[tmp_idx,]
        tmp_d <- as.dist(sqrt(dwei[tmp_idx, tmp_idx]))
        tmp_gr <- pam(tmp_d, k = k_clust, diss = TRUE, pamonce = 5, cluster.only = TRUE)
        tmp_tab <- table(tmp_B$GR, tmp_gr)
        tmp_jac <- unname(tmp_tab/(matrix(rowSums(tmp_tab), nrow = k_clust, ncol = k_clust) +
                                       matrix(colSums(tmp_tab), nrow = k_clust, ncol = k_clust, byrow = TRUE) -
                                       tmp_tab))
        cum_jac <- cum_jac + apply(tmp_jac, 1, max)
    }
    rm(samp, tmp_B, tmp_gr, tmp_tab, tmp_jac)
    cum_jac <- cum_jac/n_it
    cbi[[(k_clust - 1)]] <- cum_jac
}

CSI <- as.data.frame(t(sapply(cbi, "[", i = 1:9)))
names(CSI) <- paste("cluster", 1:9, sep = "_")
CSI$n_clusters <- 2:9
CSI <- CSI[, c(10, 1:9)]
write.table(CSI, file = "results/cluster_stability_index_k2-9.csv", sep = ",", row.names = FALSE)

############################
# Partition around medoids
############################
k_clust = 4

pam.res <- pam(d, k = k_clust, diss = TRUE, pamonce = 5)
meds <- dplyr::select(B, namesfac)[pam.res$medoids,]

B <- mutate(B, GR0 = pam.res$clustering) %>%
    mutate(GR = ifelse(GR0 == 1, "MCI 1", ifelse(GR0 == 2, "MCI 4", ifelse(GR0 == 3, "MCI 2", "MCI 3")))) %>%
    mutate(GR = factor(GR))

pb <- scores.boxplots2(B, namesfac, nam = c("Domain", "Score", "MCI Subgroup"), flip = FALSE, by.diagnosis = FALSE) +
    ylim(-2.5, 8.5)
pb

meds <- mutate(meds, GR = paste("MCI", c(1, 4, 2, 3)))
meds <- meds[, c(5,1:4)]
write.table(meds, file = "results/medoids_domainscores_MCIsubgroups_k4.csv", sep = ",", row.names = FALSE)

####################################################
# Subgroups description with controls as reference
####################################################

CN <- filter(S, DIAGNOSIS == "CN") %>%
    mutate(event = ifelse(new.n.prog > 0, 1, 0), time = ifelse(new.n.prog > 0, new.time.change, new.months.fu),
           GR0 = 0, GR = "CN")

MCICN <- B %>%
    mutate(GR = as.character(GR)) %>%
    rbind(CN) %>%
    mutate(GR = factor(GR, levels = c("CN", paste("MCI", 1:4))))

pb <- scores.boxplots2(MCICN, namesfac, nam = c("Domain", "Score", "Group"), flip = FALSE, by.diagnosis = FALSE) +
    ylim(-2.5, 8.5)
pb
ggsave("plots/Scores_MCIsubgroups_CN.png", pb, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/Scores_MCIsubgroups_CN.eps", pb, width = 16, height = 12, units = "cm", dpi = 300, bg = "transparent")

# Subgroups description
x = c("X", levels(MCICN$GR))
sdes <- data.frame(matrix(ncol = length(x), nrow = 1))
colnames(sdes) <- x

sdes$X[1] <- "female percentage"
tab <- with(MCICN, table(PTGENDER,GR))
sdes[1, 2:ncol(sdes)] <- round(tab[1,]/colSums(tab)*100, digits = 2)

for(varn in c("AGE", namesfac)){
    fmla <- as.formula(paste(varn, "GR", sep = " ~ "))
    tmp <- MCICN %>%
        group_by(GR) %>% 
        summarise(mean = mean(!!sym(varn)), sd = sd(!!sym(varn)))
    tmpdf <- data.frame(X = paste(c("mean", "std"), tolower(varn))) 
    tmpdf <- cbind(tmpdf, t(round(tmp[,2:3], digits = 2)))
    rownames(tmpdf) <- NULL 
    colnames(tmpdf) <- c("X", levels(MCICN$GR))
    sdes <- rbind(sdes, tmpdf)    
}
sdes
write.table(sdes, file = "results/description_MCIsubgroups_k4.csv", sep = ",", row.names = FALSE)

# Comparison between pairs of groups
# 10 comparisons between groups x 6 domains

varn <- namesfac[1]
fmla <- as.formula(paste(varn, "GR", sep = " ~ "))
sgr <- c("MCI 1", "MCI 2")

tmp <- filter(MCICN, GR %in% sgr)
N <- nrow(tmp)

MW1 <- wilcox.test(fmla, tmp)
Z = qnorm(MW1$p.value/2)
r = abs(Z)/sqrt(N)
round(r, 2)
MW1$p.value*60

#####################
# Survival Analysis
#####################
surv_object <- Surv(time = B$time, event = B$event)
fitsv <- survfit(surv_object ~ GR, data = B)

# Kaplan-Meier curves
pl <- survival.plot(fitsv, B, k = k_clust)
pl
ggsave("plots/KMcurves_MCIsubgroups_k4.png", pl, width = 12, height = 10, units = "cm", dpi = 300, bg = "transparent")
ggsave("plots/KMcurves_MCIsubgroups_k4.eps", pl, width = 12, height = 10, units = "cm", dpi = 300, bg = "transparent")

# Cox model
fit.coxph <- coxph(surv_object ~ GR + PTGENDER + PTEDUCAT + AGE, data = B)
ggforest(fit.coxph, data = B)

survdiff(surv_object ~ GR, data = B)
pairwise_survdiff(Surv(time, event) ~ GR, data = B, p.adjust.method = "fdr")


