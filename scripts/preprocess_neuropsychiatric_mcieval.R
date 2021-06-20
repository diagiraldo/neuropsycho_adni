##########################################
### Extract Neuropsychiatric Inventory ###
##########################################
### D. Giraldo, Jun2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(ADNIMERGE)

# Load data
load("processed_data/domain_scores_srb_firstvisit.RData")
B <- filter(S, set == "EVAL")
PVIS <- select(B, RID, VISCODE, COLPROT) 

# Check GDS info
tmp_gds <- filter(gdscale, RID %in% B$RID) %>%
    select(COLPROT:VISCODE, GDTOTAL) %>%
    mutate(VISCODE = ifelse(VISCODE == "sc", "bl", VISCODE))

PVIS <- merge(PVIS, tmp_gds, all.x = TRUE)
rm(tmp_gds)

# Check if there is NPI-Q info for pt at taken visit
tmp_npiq <- filter(npiq, RID %in% B$RID) %>%
    select(COLPROT:VISCODE, NPISCORE) %>%
    rename(NPIQ = NPISCORE)

PVIS <- merge(PVIS, tmp_npiq, all.x = TRUE)
with(PVIS, table(is.na(NPIQ), COLPROT))
# Missing for visits during ADNI2 and ADNI 3
rm(tmp_npiq)

# Check if there is NPI info for pt at taken visit
tmp_npi <- filter(npi, RID %in% B$RID) %>%
    select(COLPROT:VISCODE, NPITOTAL)

PVIS <- merge(PVIS, tmp_npi, all.x = TRUE)
with(PVIS, table(is.na(NPITOTAL), COLPROT))

# Compute NPI-Q from NPI 
sel_npi <- paste0("NPI", LETTERS[1:12])

tmp_npi <- filter(npi, RID %in% B$RID) %>%
    select(COLPROT:VISCODE, all_of(sel_npi), matches("\\dB$"), NPITOTAL)
names(tmp_npi)[17:28] <- paste0(sel_npi, "SEV")

npisevnum <- function(sevec) {0 + 1 * grepl("1:", sevec) + 2 * grepl("2:", sevec) + 3 * grepl("3:", sevec)}

tmp_npi <- tmp_npi %>%
    mutate(across(NPIASEV:NPILSEV, npisevnum)) %>%
    mutate(NPIQ.calc = rowSums(select(., NPIASEV:NPILSEV))) %>%
    select(COLPROT:VISCODE, NPIQ.calc)

PVIS <- merge(PVIS, tmp_npi, all.x = TRUE) %>%
    mutate(NPIQSCORE = ifelse(is.na(NPIQ), NPIQ.calc, NPIQ)) %>%
    select(RID:SITEID, GDTOTAL, NPIQSCORE)

save(PVIS, file = "processed_data/neuropsychiatric_MCIinEval.RData")