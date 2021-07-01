#######################################################
### Extract Scores to compare with other composites ###
#######################################################
### D. Giraldo, Jun2021

setwd("~/neuropsycho_adni")

library(dplyr)
library(ggplot2)
library(ADNIMERGE)

# Load data
load("processed_data/neuropsycho_seltests.RData")

# ADAS.Tree from Llano, 2011.
# Composite score from Huang, 2015.
# Composite scores in Raghvan,2013: ADAS3, CC1, CC2, 
DF <- DF %>%
    mutate(ADAS.tree = 1.05*Q1SCORE + 0.38*Q2SCORE + 0*Q3SCORE + 1.17*Q4SCORE +
               0.61*Q5SCORE + 0.13*Q6SCORE + 1.13*Q7SCORE + 0.41*Q8SCORE +
               0.54*Q9SCORE + 0.49*Q10SCORE + 0.69*Q11SCORE + 0.39*Q12SCORE + 0.68*Q13SCORE,
           comp.Huang = Q1SCORE + Q4SCORE + Q7SCORE + CDRSB + FAQTOTAL,
           ADAS3 = Q1SCORE + Q4SCORE + Q7SCORE,
           CC1 = ADAS3 + (75-RAVLT.IMMED) + (30 - MMTOTAL),
           CC2 = ADAS3 + CDMEMORY,
           CFC1 = CC1 + FAQTOTAL, CFC2 = CC2 + FAQTOTAL, CFC3 = CDMEMORY + FAQTOTAL,
           Forget.index = (LDELTOTAL - LIMMTOTAL)/LIMMTOTAL*100)
selcompsc <- c("ADAS.tree", "comp.Huang", "CC1", "CC2", "CFC1", "CFC2", "CFC3")

# Feature selection in Pereira, 2018: inclutes totals and parts
# RAVLT.IMMED = AVTOT15
selfeaturesPer <- c("TRABSCOR", "Forget.index", "RAVLT.IMMED", "TOTAL13", "TRAASCOR", "AVTOT6", "LIMMTOTAL",
                    "CATANIMSC", "AVDEL30MIN", "FAQTOTAL", "LDELTOTAL", "MOCADLREC", "AVDELTOT", "BNTTOTAL",
                    "Q4SCORE", "Q8SCORE", "MMTOTAL", "Q1SCORE", "MOCAFLUEN", "CDORIENT", "CDHOME", "AVTOTB")
# In Pereira,2018 np features + age

# Select variables of interest
DF <- select(DF, ORIGPROT:VISCODE, selcompsc, selfeaturesPer, CDGLOBAL, CDRSB)

# Keep data used
load("processed_data/selsubscores_firstvisit.RData")

BCOMP <- select(A, RID, VISCODE, COLPROT) 
BCOMP <- merge(BCOMP, DF, all.x = TRUE)
rm(A, DF)

save(BCOMP, selcompsc, selfeaturesPer, file = "processed_data/other_composites_features.RData")
