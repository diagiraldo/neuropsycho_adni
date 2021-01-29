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

###############
# Split groups
###############

splitFA <- 0.5
set.seed(1995) # Release of Toy Story 1 ;)

setFA <- A$RID[A$new.months.fu < 12 | A$DIAGNOSIS == "Dementia"]
setFA <- c(setFA, sample(setdiff(A$RID, setFA), size = round(length(A$RID)*splitFA)-length(setFA)))
setEVAL <- setdiff(A$RID, setFA)
A <- mutate(A, set1 = ifelse(RID %in% setFA, "FA", "EVAL"))
rm(setFA, setEVAL)

# Include covariates and filter data
A <- A %>% 
    merge(dplyr::select(D, RID, PTEDUCAT, AGE.bl, PTGENDER), all.x = TRUE) %>%
    mutate(AGE = AGE.bl + VISMONTH/12, PTGENDER = as.factor(PTGENDER)) %>%
    filter(!is.na(AGE.bl)) %>%
    filter(DIAGNOSIS != "Dementia")
rm(D, DX)

# Reasign groups
setSTD <- A$RID[A$DIAGNOSIS == "CN" & A$set1 == "EVAL"]
setCFA <- setdiff(A$RID[A$DIAGNOSIS == "CN"], setSTD)
setEVAL <- A$RID[A$DIAGNOSIS == "MCI" & A$set1 == "EVAL"]

splitmci <- 0.6
nmci <- sum(A$DIAGNOSIS == "MCI") 
set.seed(123) 
tmp <- sample(setdiff(A$RID[A$DIAGNOSIS == "MCI" & A$new.months.fu > 12], setEVAL),
              size = round(nmci*splitmci - length(setEVAL)))
setEVAL <- c(setEVAL, tmp)
setCFA <- c(setCFA, setdiff(A$RID[A$DIAGNOSIS == "MCI"], setEVAL))

A <- mutate(A, set = ifelse(RID %in% setSTD, 
                            "STD", ifelse(RID %in% setCFA, "FA", "EVAL")))

with(A, table(DIAGNOSIS, set))
rm(tmp, setSTD, setCFA, setEVAL)

# Save data first visit
save(A, it, file = "processed_data/selsubscores_firstvisit.RData")

rm(list = ls())


