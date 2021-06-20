################################################
### Select items and split data for analysis ###
################################################
### D. Giraldo, Oct 2019

setwd("~/neuropsycho_adni")

library(dplyr)
library(reshape2)
library(tidyr)

source("scripts/auxiliar_functions_adni_data.R")

# Load pre-processed data
load("processed_data/neuropsycho_reversed_seltests.RData")
load("processed_data/diagnostic_data.RData")

sel_tests <- c("ADAS", "MMSE", "MOCA", "CLOCK", "TMT", "LM", "CATFL", "AVLT", "BNTMINT")
it <- unname(unlist(selitems[sel_tests]))

# Extract info and first visit in DF per subject
A <- dplyr::select(DF, c(1:5, it)) %>%
  merge(dplyr::select(DX, c(1:4, DIAGNOSIS, VISMONTH)), all.x = TRUE) %>%
  filter(complete.cases(.)) %>%
  merge(firstvisit(.)) %>%
  filter(is.firstvis) %>%
  merge(dplyr::select(D, c(RID, months.fu:vis.1stchange)), all.x = TRUE) %>%
  mutate(new.months.fu = months.fu - VISMONTH)

# Adjust data for selected visit
DX <- merge(DX, dplyr::select(A, c(RID, FIRSTVIS)), all.x = TRUE) %>%
  filter(VISMONTH >= FIRSTVIS) %>%
  arrange(VISMONTH)
A$new.n.dx = sapply(A$RID, function(x) length(unique(DX$DIAGNOSIS[DX$RID == x & !is.na(DX$DIAGNOSIS)])))
A$new.is.stable = ifelse(A$new.n.dx > 0, (A$new.n.dx == 1), NA)
A$new.str.dx = sapply(A$RID, function(x) paste(rle(as.character(DX$DIAGNOSIS[DX$RID == x & !is.na(DX$DIAGNOSIS)]))$values, collapse = "; "))
A$new.n.prog = sapply(A$RID, function(x) sum(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) > 0))
A$new.n.rev = sapply(A$RID, function(x) sum(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) < 0))
A$new.vis.change = sapply(A$RID, function(x) DX$VISMONTH[DX$RID == x & !is.na(DX$DIAGNOSIS)][which(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) != 0)[1] + 1])
A$new.time.change = A$new.vis.change - A$FIRSTVIS

# Include covariates and filter data
A <- A %>% 
    merge(dplyr::select(D, RID, PTEDUCAT, AGE.bl, PTGENDER), all.x = TRUE) %>%
    mutate(AGE = AGE.bl + VISMONTH/12, PTGENDER = as.factor(PTGENDER)) %>%
    filter(!is.na(AGE.bl)) %>%
    filter(DIAGNOSIS != "Dementia")
rm(D, DX)

# Include APOE info 
library(ADNIMERGE)
AP <- rbind(dplyr::select(apoeres, ORIGPROT, RID, APGEN1, APGEN2),
            dplyr::select(apoego2, ORIGPROT, RID, APGEN1, APGEN2),
            dplyr::select(apoe3, ORIGPROT, RID, APGEN1, APGEN2)) %>%
    filter(RID %in% A$RID) %>%
    mutate(apoee4 = (APGEN1 == "4" | APGEN2 == "4"))
A <- merge(A, dplyr::select(AP, RID, apoee4), all.x = TRUE)

# Data description

table(A$DIAGNOSIS)

x = c("X", levels(A$DIAGNOSIS))
sdes <- data.frame(matrix(ncol = length(x), nrow = 2))
colnames(sdes) <- x

sdes$X[1] <- "female percentage"
tab <- with(A, table(PTGENDER,DIAGNOSIS))
sdes[1, 2:ncol(sdes)] <- round(tab[1,]/colSums(tab)*100, digits = 2)

sdes$X[2] <- "APOE-e4 percentage"
tab <- with(A, table(apoee4,DIAGNOSIS))
sdes[2, 2:ncol(sdes)] <- round(tab[2,]/colSums(tab)*100, digits = 2)

sdes

#########################
# Include data partition
#########################

tmp <- read.csv("data_partition.csv", stringsAsFactors = FALSE)
names(tmp)[2] <- "set"
A <- merge(A, tmp) 

with(A, table(set, DIAGNOSIS))

# Save data first visit
save(A, it, file = "processed_data/selsubscores_firstvisit.RData")

# Data description per set

# Include test totals in data description
# Load test totals for baseline 
load("processed_data/neuropsycho_seltests.RData")
sel_tests <- c("ADAS", "MMSE", "CDR")
tot_tests <- unname(unlist(seltotals[sel_tests]))
A <- merge(A, dplyr::select(DF, c(1:5, tot_tests)), all.x = TRUE)
rm(DF, selitems, seltotals)

tmp <- filter(A, DIAGNOSIS == "MCI" & set == "EVAL")
tab <- with(tmp, table(PTGENDER))
round(tab/sum(tab)*100, digits = 1)
round(mean(tmp$AGE), digits = 1)
round(sd(tmp$AGE), digits = 1)
tab <- with(tmp, table(apoee4))
round(tab/sum(tab)*100, digits = 1)
round(mean(tmp$CDRSB, na.rm = TRUE), digits = 2)
round(sd(tmp$CDRSB, na.rm = TRUE), digits = 2)
round(mean(tmp$MMTOTAL, na.rm = TRUE), digits = 1)
round(sd(tmp$MMTOTAL, na.rm = TRUE), digits = 1)
round(mean(tmp$TOTAL13, na.rm = TRUE), digits = 1)
round(sd(tmp$TOTAL13, na.rm = TRUE), digits = 1)

rm(list = ls())


