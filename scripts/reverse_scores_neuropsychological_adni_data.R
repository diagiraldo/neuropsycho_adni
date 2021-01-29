######################################################
### Script to reverse some neuropsychological data ###
############### measures of impairment ###############
######################################################
### D. Giraldo, Sept 2019

setwd("~/neuropsycho_adni/")
library(dplyr)

load("processed_data/neuropsycho_seltests.RData")

# Reverse MMSE items
DF <- DF %>%
  mutate(MMORITIME = 5-MMORITIME, MMORISPACE = 5-MMORISPACE, MMREGI = 3-MMREGI, MMRECALL = 3-MMRECALL,
         MMSPELLBKW = 5-MMSPELLBKW, MMNAM = 2-MMNAM, MMCOMMAND = 3-MMCOMMAND,
         MMREPEAT = 1-MMREPEAT, MMREAD = 1-MMREAD, MMWRITE = 1-MMWRITE, MMDRAW = 1-MMDRAW,
         MMTOTAL = 30-MMTOTAL)

# Reverse Logical Memory 25-x
DF <- DF %>%
  mutate(LIMMTOTAL = 25-LIMMTOTAL, LDELTOTAL = 25-LDELTOTAL)

# Reverse Rey Auditory Verbal Learning Test and Delayed version 15-x
DF <- DF %>%
  mutate_at(selitems$AVLT[2:5], function(x) 15-x) %>%
  mutate(RAVLT.IMMED = 75-RAVLT.IMMED, 
         RAVLT.LEARN = -RAVLT.LEARN)

# Reverse Clock Drawing - Copying
DF <- DF %>%
  mutate_at(selitems$CLOCK, function(x) 5-x)

# Reverse Category fluency 
DF <- DF %>%
  mutate(CATANIMSC = 60-CATANIMSC)

# Reverse Boston Naming Test 30-x
DF <- DF %>%
  mutate_at(selitems$BNT, function(x) 30-x) %>%
  mutate(BNTTOTAL = 30-BNTTOTAL)

# Reverse MINT 32-x
DF <- DF %>%
  mutate_at(selitems$MINT, function(x) 32-x) %>%
  mutate(MINTTOTAL = 32-MINTTOTAL)

# Reverse BNTMINT 30-x
DF <- DF %>%
  mutate_at(selitems$BNTMINT, function(x) 30-x) %>%
  mutate(BMTOTAL = 30-BMTOTAL)

# Reverse MoCA 
DF <- DF %>%
  mutate(TRAILS = 1-TRAILS, CUBE = 1-CUBE, MOCACLOCK = 3-MOCACLOCK, MOCANAM = 3-MOCANAM,
         MOCADIG = 2-MOCADIG, MOCALET = 1-MOCALET, MOCASERIAL = 3-MOCASERIAL, MOCAREP = 2-MOCAREP,
         MOCAFLUEN = 1-MOCAFLUEN, MOCAABS = 2-MOCAABS, MOCADLREC = 5-MOCADLREC, MOCAORI = 6-MOCAORI,
         MOCATOTAL = 30-MOCATOTAL)

# Reverse FCI
DF <- DF %>%
  mutate(MC_SCORE = 1-MC_SCORE, FC_SCORE = 8-FC_SCORE, SNGLCHKSCORE = 20-SNGLCHKSCORE,
         COMPLXCHKSCORE = 28-COMPLXCHKSCORE, BANKSTATSCORE = 14-BANKSTATSCORE, FCITOTAL= 74-FCITOTAL)

save(DF, selitems, seltotals, file = "processed_data/neuropsycho_reversed_seltests.RData")

rm(DF, selitems, seltotals)