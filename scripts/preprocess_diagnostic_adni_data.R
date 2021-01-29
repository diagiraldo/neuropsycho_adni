##################################################
### Script to pre-process ADNI diagnostic data ###
###### Using ADNIMERGE installed on Sep 2019 #####
##################################################
### D. Giraldo, Oct 2019

setwd("~/neuropsycho_adni/")
library(dplyr)
library(ADNIMERGE)

source("scripts/auxiliar_functions_adni_data.R")

# Info in diagnostic summary, not every visit 
DX <- dxsum %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("uns1", "y1", "y2"))) %>%
  dplyr::select(c(1:5, DIAGNOSIS, DXCHANGE, DXCONFID, DXMDUE, DXDDUE, DXDEP))

DXA <- adnimerge %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("m0", "y1"))) %>%
  dplyr::select(ORIGPROT, COLPROT, RID, VISCODE, DX, DX.bl, AGE, PTGENDER, PTEDUCAT, Years.bl, Month.bl)

DX <- merge(DX, DXA, all = TRUE) %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("CN", "MCI", "Dementia")),
         cod.DX = as.numeric(DIAGNOSIS),
         VISMONTH = as.numeric(ifelse(VISCODE == "sc", -1, ifelse(VISCODE == "bl", 0, gsub("m", "", VISCODE))))) %>%
  arrange(VISMONTH)
rm(DXA)

rids_bl <- unique(filter(DX, VISCODE == "bl")$RID)
rids_scnobl <- setdiff(unique(filter(DX, VISCODE == "sc")$RID), rids_bl)

D <- filter(DX, (VISCODE == "bl" & RID %in% rids_bl) | (VISCODE == "sc" & RID %in% rids_scnobl)) %>%
  dplyr::select(c(1:5, DIAGNOSIS, DX.bl, AGE, PTGENDER, PTEDUCAT)) %>%
  rename(DX3.bl = DIAGNOSIS, AGE.bl = AGE) %>%
  mutate(DX3.bl = factor(DX3.bl, levels = c("CN", "MCI", "Dementia")),
         cod.DX3.bl = as.numeric(DX3.bl))

D$n.visits = sapply(D$RID, function(x) length(unique(DX$VISCODE[DX$RID == x])))
D$months.fu = sapply(D$RID, function(x) max(DX$VISMONTH[DX$RID == x], na.rm = TRUE))
D$n.dx = sapply(D$RID, function(x) length(unique(DX$DIAGNOSIS[DX$RID == x & !is.na(DX$DIAGNOSIS)])))
D$is.stable = ifelse(D$n.dx > 0, (D$n.dx == 1), NA)
D$str.dx = sapply(D$RID, function(x) paste(rle(as.character(DX$DIAGNOSIS[DX$RID == x & !is.na(DX$DIAGNOSIS)]))$values, collapse = "; "))
D$n.prog = sapply(D$RID, function(x) sum(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) > 0))
D$n.rev = sapply(D$RID, function(x) sum(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) < 0))
D$vis.1stchange = sapply(D$RID, function(x) DX$VISMONTH[DX$RID == x & !is.na(DX$DIAGNOSIS)][which(diff(DX$cod.DX[DX$RID == x & !is.na(DX$DIAGNOSIS)]) != 0)[1] + 1])

save(D, DX, file = "processed_data/diagnostic_data.RData")
rm(DX, D, rids_bl, rids_scnobl)
