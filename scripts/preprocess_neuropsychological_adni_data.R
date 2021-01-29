#####################################################
### Script to pre-process Neuropsychological data ###
####### Using ADNIMERGE installed on Sep 2019 #######
#####################################################
### D. Giraldo, Sept 2019

setwd("~/neuropsycho_adni/")

library(dplyr)
library(ADNIMERGE)

# Function to detect double registers RID-VISCODE
reg2del <- function(DF){
  tmp <- with(DF, table(RID, VISCODE))
  idx <- which(tmp > 1, arr.ind = TRUE)
  nam <- dimnames(tmp[idx[,1],idx[,2]])
  l <- paste(DF$RID, DF$VISCODE, sep=":") %in% paste(nam$RID, nam$VISCODE, sep=":")
  return(l)
}

# Init Data frames and lists
DF <- data.frame()
selitems <- list()
seltotals <- list()

# ADAS-Cog 
ADAS <- adas %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("uns1", "y1", "y2"))) %>%
  dplyr::select(c(1:5, ends_with("SCORE"), TOTAL13)) %>%
  filter(!is.na(TOTAL13)) %>%
  unique()

DF <- ADAS
selitems$ADAS <- paste("Q", 1:13, "SCORE", sep = "")
seltotals$ADAS <- c("TOTAL13")
#rm(ADAS)

# CDR
CDR <- cdr %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  mutate(VISCODE = ifelse(VISCODE == "sc", "bl", VISCODE)) %>%
  dplyr::select(c(1:5, CDCARE, CDCOMMUN, CDHOME, CDJUDGE, CDMEMORY, CDORIENT, CDGLOBAL, CDRSB)) %>%
  filter(!is.na(CDRSB)) %>%
  unique()

DF <- merge(DF, CDR, all = TRUE)
selitems$CDR <- c("CDCARE", "CDCOMMUN", "CDHOME", "CDJUDGE", "CDMEMORY", "CDORIENT")
seltotals$CDR <- c("CDGLOBAL", "CDRSB")
#rm(CDR)

# FAQ
faq2number <- function(faqchr){
  n <- as.numeric(gsub("[^0-9]", "",  faqchr))
  return(n)
}

FAQ <- faq %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  dplyr::select(c(1:5, FAQBEVG:FAQTV)) %>%
  mutate_at(vars(FAQBEVG:FAQTV), faq2number) %>%
  mutate(FAQTOTAL = rowSums(select(., FAQBEVG:FAQTV))) %>%
  filter(!is.na(FAQTOTAL)) %>%
  unique()

DF <- merge(DF, FAQ, all = TRUE)
selitems$FAQ <- c("FAQBEVG", "FAQEVENT", "FAQFINAN", "FAQFORM", "FAQGAME",
                  "FAQMEAL", "FAQREM", "FAQSHOP", "FAQTRAVL", "FAQTV")
seltotals$FAQ <- c("FAQTOTAL")
#rm(FAQ)

# MMSE (to be reversed)
mmse_itnames <- c("MMYEAR", "MMSEASON", "MMDATE", "MMMONTH", "MMDAY",
                  "MMSTATE", "MMAREA", "MMCITY", "MMHOSPIT", "MMFLOOR",
                  "WORD1", "WORD2", "WORD3",
                  "WORD1DL", "WORD2DL", "WORD3DL",
                  "MMWATCH", "MMPENCIL",
                  "MMHAND", "MMFOLD", "MMONFLR",
                  "MMREPEAT", "MMREAD", "MMWRITE", "MMDRAW",
                  "MMLTR1", "MMLTR2", "MMLTR3", "MMLTR4", "MMLTR5")
mmse_itwo <- c("MMBALL", "MMFLAG", "MMTREE",
               "MMW", "MMO", "MMR", "MML", "MMD",
               "MMBALLDL", "MMFLAGDL", "MMTREEDL")

MMSE <- mmse %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  mutate(VISCODE = ifelse(VISCODE == "sc", "bl", VISCODE)) %>%
  dplyr::select(1:5, mmse_itnames, mmse_itwo, MMSCORE) %>%
  mutate(WORD1 = ifelse(COLPROT == "ADNI3", WORD1, MMBALL),
         WORD2 = ifelse(COLPROT == "ADNI3", WORD2, MMFLAG),
         WORD3 = ifelse(COLPROT == "ADNI3", WORD3, MMTREE),
         WORD1DL = ifelse(COLPROT == "ADNI3", WORD1DL, MMBALLDL),
         WORD2DL = ifelse(COLPROT == "ADNI3", WORD2DL, MMFLAGDL),
         WORD3DL = ifelse(COLPROT == "ADNI3", WORD3DL, MMTREEDL),
         MMLTR1 = ifelse(COLPROT == "ADNI3", ifelse(MMLTR1 %in% c("d", "D"), "Correct", "Incorrect"), as.character(MMD)),
         MMLTR2 = ifelse(COLPROT == "ADNI3", ifelse(MMLTR2 %in% c("l", "L"), "Correct", "Incorrect"), as.character(MML)),
         MMLTR3 = ifelse(COLPROT == "ADNI3", ifelse(MMLTR3 %in% c("r", "R"), "Correct", "Incorrect"), as.character(MMR)),
         MMLTR4 = ifelse(COLPROT == "ADNI3", ifelse(MMLTR4 %in% c("o", "O"), "Correct", "Incorrect"), as.character(MMO)),
         MMLTR5 = ifelse(COLPROT == "ADNI3", ifelse(MMLTR5 %in% c("w", "W"), "Correct", "Incorrect"), as.character(MMW))) %>%
  dplyr::select(-mmse_itwo) %>%
  mutate_at(mmse_itnames, function(x) as.numeric(ifelse(grepl("Incorrect", x), 0, ifelse(grepl("[C|c]orrect", x), 1, x)))) %>%
  mutate(MMORITIME = MMYEAR + MMSEASON + MMMONTH + MMDATE + MMDAY,
         MMORISPACE = MMSTATE + MMAREA + MMCITY + MMHOSPIT + MMFLOOR,
         MMREGI = WORD1 + WORD2 + WORD3,
         MMRECALL = WORD1DL + WORD2DL + WORD3DL,
         MMSPELLBKW = MMLTR1 + MMLTR2 + MMLTR3 + MMLTR4 + MMLTR5,
         MMNAM = MMWATCH + MMPENCIL,
         MMCOMMAND = MMHAND + MMFOLD + MMONFLR) %>%
  mutate(MMTOTAL = rowSums(select(., MMORITIME:MMCOMMAND, MMREPEAT, MMREAD, MMWRITE, MMDRAW))) %>%
  filter(!is.na(MMTOTAL)) %>%
  unique()

DF <- merge(DF, MMSE, all = TRUE)
selitems$MMSE <- c("MMORITIME", "MMORISPACE", "MMREGI", "MMRECALL", "MMSPELLBKW", 
                  "MMNAM", "MMCOMMAND", "MMREPEAT", "MMREAD", "MMWRITE", "MMDRAW")
seltotals$MMSE <- c("MMTOTAL")
# rm(MMSE)

# Logical Memory (to be reversed) 25-x
LM <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  dplyr::select(c(1:5, LIMMTOTAL, LDELTOTAL)) %>%
  mutate(VISCODE = ifelse(VISCODE == "sc", "bl", VISCODE)) %>%
  filter(!is.na(LIMMTOTAL) & !is.na(LDELTOTAL)) %>%
  unique()

DF <- merge(DF, LM, all = TRUE)
selitems$LM <- c("LIMMTOTAL", "LDELTOTAL") # No items, just totals
seltotals$LM <- c("LIMMTOTAL", "LDELTOTAL") # No items, just totals
# rm(LM)

# Rey Auditory Verbal Learning Test and Delayed version (to be reversed) 15-x
AVLT <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, starts_with("AVTOT"), AVDEL30MIN, AVDELTOT)) %>%
  mutate(RAVLT.IMMED = AVTOT1 + AVTOT2 + AVTOT3 + AVTOT4 + AVTOT5,
         RAVLT.LEARN = AVTOT5 - AVTOT1,
         RAVLT.FORGET = AVTOT5 - AVDEL30MIN) %>%
  unique()

DF <- merge(DF, AVLT, all = TRUE)
selitems$AVLT <- c("RAVLT.IMMED", paste("AVTOT", c(6, "B"), sep = ""), "AVDEL30MIN", "AVDELTOT")
# Some subjects has only a subset of trials, why?
seltotals$AVLT <- c("RAVLT.IMMED", "RAVLT.LEARN", "RAVLT.FORGET")
# rm(AVLT)

# Clock Drawing - Copying (to be reversed)
clock_itnames <- c("CLOCKCIRC", "CLOCKSYM", "CLOCKNUM", "CLOCKHAND", "CLOCKTIME",
                   "COPYCIRC", "COPYSYM", "COPYNUM", "COPYHAND", "COPYTIME")

CLOCK <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, clock_itnames, CLOCKSCOR, COPYSCOR)) %>%
  mutate_at(clock_itnames, function(x) as.numeric(ifelse(x %in% c("Correct", "Yes"), 1, ifelse(x %in% c("Incorrect", "No"), 0, x)))) %>%
  filter(!is.na(rowSums(select(., clock_itnames[1:5]))) | !is.na(rowSums(select(., clock_itnames[6:10])))) %>%
  unique()

DF <- merge(DF, CLOCK, all = TRUE)
selitems$CLOCK <- c("CLOCKSCOR", "COPYSCOR")
seltotals$CLOCK <- c("CLOCKSCOR", "COPYSCOR")
# rm(CLOCK)

# Category Fluency (to be reversed)
CATFL <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, CATANIMSC, CATVEGESC)) %>%
  filter(!is.na(CATANIMSC) | !is.na(CATVEGESC)) %>%
  unique()

DF <- merge(DF, CATFL, all = TRUE)
selitems$CATFL <- c("CATANIMSC") # No items, just totals
seltotals$CATFL <- c("CATANIMSC") # No items, just totals
# rm(CATFL)

# Trail Making Test, Scores are in seconds A: max 150s, B: max 300s 
TMT <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, starts_with("TRAA"), starts_with("TRAB"))) %>%
  mutate(TRAASCOR = ifelse(TRAASCOR > 150, 150, TRAASCOR),
         TRABSCOR = ifelse(TRABSCOR > 300, 300, TRABSCOR)) %>%
  filter(!is.na(TRAASCOR) | !is.na(TRABSCOR)) %>%
  unique()

DF <- merge(DF, TMT, all = TRUE)
selitems$TMT <- c("TRAASCOR", "TRABSCOR") # No items, just totals
seltotals$TMT <- c("TRAASCOR", "TRABSCOR") # No items, just totals
# rm(TMT)

# American National Adult Reading Test
ANART <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, ANARTERR)) %>%
  filter(!is.na(ANARTERR)) %>%
  unique()

DF <- merge(DF, ANART, all = TRUE)
selitems$ANART <- c("ANARTERR") # No items, just totals
seltotals$ANART <- c("ANARTERR") # No items, just totals
# rm(ANART)

# Boston Naming Test (to be reversed) 30-x
BNT <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, BNTSPONT, BNTCSTIM, BNTTOTAL)) %>%
  filter(!is.na(BNTSPONT)) %>%
  mutate(BNTCSTIM = ifelse(is.na(BNTCSTIM), 0, BNTCSTIM)) %>%
  mutate(BNTTOTAL = ifelse(BNTCSTIM == 0, BNTSPONT, BNTTOTAL)) %>%
  filter((BNTSPONT + BNTCSTIM) == BNTTOTAL) %>%
  unique()

DF <- merge(DF, BNT, all = TRUE)
selitems$BNT <- c("BNTSPONT", "BNTCSTIM")
seltotals$BNT <- c("BNTTOTAL") 
# rm(BNT)

# Mulitlingual	Naming	test	(MINT) (to be reversed) only in ADNI 3
# It replaces the instrument BNT from ADNI2 Neuropsych battery
MINT <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, MINTUNCUED, MINTSEMCUE, MINTTOTAL)) %>%
  filter(!is.na(MINTUNCUED)) %>%
  mutate(MINTSEMCUE = ifelse(is.na(MINTSEMCUE), 0, MINTSEMCUE)) %>%
  mutate(MINTTOTAL = ifelse(MINTSEMCUE == 0, MINTUNCUED, MINTTOTAL)) %>%
  filter((MINTUNCUED + MINTSEMCUE) == MINTTOTAL) %>%
  unique()

DF <- merge(DF, MINT, all = TRUE)
selitems$MINT <- c("MINTUNCUED", "MINTSEMCUE")
seltotals$MINT <- c("MINTTOTAL")
# rm(MINT)

# Naming	test	(BNTMINT) (to be reversed) over 30
BNTMINT <- neurobat %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2", "sc"))) %>%
  dplyr::select(c(1:5, BNTSPONT, BNTCSTIM, BNTTOTAL, MINTUNCUED, MINTSEMCUE, MINTTOTAL)) %>%
  filter(!is.na(BNTSPONT) | !is.na(MINTUNCUED)) %>%
  mutate(BNTCSTIM = ifelse(is.na(BNTCSTIM), 0, BNTCSTIM),
         MINTSEMCUE = ifelse(is.na(MINTSEMCUE), 0, MINTSEMCUE), 
         BNTTOTAL = ifelse(BNTCSTIM == 0, BNTSPONT, BNTSPONT + BNTCSTIM),
         MINTTOTAL = ifelse(MINTSEMCUE == 0, MINTUNCUED, MINTUNCUED + MINTSEMCUE),
         BMNOCUE = ifelse(COLPROT == "ADNI3", MINTUNCUED/32*30, BNTSPONT),
         BMCUED = ifelse(COLPROT == "ADNI3", MINTSEMCUE/32*30, BNTCSTIM),
         BMTOTAL = ifelse(COLPROT == "ADNI3", MINTTOTAL/32*30, BNTTOTAL)) %>%
  filter(BMTOTAL <= 30) %>%
  dplyr::select(c(1:5, BMNOCUE, BMCUED, BMTOTAL)) %>%
  unique()

DF <- merge(DF, BNTMINT, all = TRUE)
selitems$BNTMINT <- c("BMNOCUE", "BMCUED")
seltotals$BNTMINT <- c("BMTOTAL")
# rm(BNTMINT)

# MoCA (to be reversed)
moca_itnames <- c("TRAILS", "CUBE",
                  "CLOCKCON", "CLOCKNO", "CLOCKHAN", # Up to 3 pts
                  "LION", "RHINO", "CAMEL", # Up to 3 pts
                  "IMMT1W1", "IMMT1W2", "IMMT1W3", "IMMT1W4", "IMMT1W5", #No points
                  "IMMT2W1", "IMMT2W2", "IMMT2W3", "IMMT2W4", "IMMT2W5", #No points
                  "DIGFOR", "DIGBACK", # Up to 2 pts
                  "LETTERS", #ERRORS, INCORRECT WHEN >=2 (1pt)
                  "SERIAL1", "SERIAL2", "SERIAL3", "SERIAL4", "SERIAL5", # 4-5 correct: 3pts, 2-3 correct: 2pts, 1 correct: 1pt, 0:0
                  "REPEAT1", "REPEAT2", # Up to 2 pts
                  "FFLUENCY", # 1pt if >= 11
                  "ABSTRAN", "ABSMEAS", # Up to 2 pts
                  "DELW1", "DELW2", "DELW3", "DELW4", "DELW5", # Up to 5 pts
                  "DATE", "MONTH", "YEAR", "DAY", "PLACE", "CITY") # Up to 6 pts

MOCA <- moca %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  dplyr::select(1:5, moca_itnames) %>%
  mutate_at(moca_itnames[c(1:20, 22:28, 30:31, 37:42)], function(x) as.numeric(ifelse(x == "Correct", 1, ifelse(x == "Incorrect", 0, x)))) %>%
  mutate_at(vars(starts_with("DELW")), function(x) as.numeric(ifelse(is.na(x), NA, ifelse(x == "Correct with No Cue", 1, 0)))) %>%
  mutate(IMMT1 = IMMT1W1 + IMMT1W2 + IMMT1W3 + IMMT1W4 + IMMT1W5,
         IMMT2 = IMMT2W1 + IMMT2W2 + IMMT2W3 + IMMT2W4 + IMMT2W5,
         SERIAL = SERIAL1 + SERIAL2 + SERIAL3 + SERIAL4 + SERIAL5,
         MOCACLOCK = CLOCKCON + CLOCKNO + CLOCKHAN,
         MOCANAM = LION + RHINO + CAMEL,
         MOCADIG = DIGFOR + DIGBACK,
         MOCALET = as.numeric(LETTERS < 2),
         MOCASERIAL = ifelse(SERIAL >= 4, 3, ifelse(SERIAL >= 2, 2, ifelse(SERIAL == 1, 1, 0))),
         MOCAREP = REPEAT1 + REPEAT2,
         MOCAFLUEN = as.numeric(FFLUENCY >= 11),
         MOCAABS = ABSTRAN + ABSMEAS,
         MOCADLREC = DELW1 + DELW2 + DELW3 + DELW4 + DELW5,
         MOCAORI = DATE + MONTH + YEAR + DAY + PLACE + CITY) %>%
  mutate(MOCATOTAL = rowSums(select(., TRAILS, CUBE, MOCACLOCK:MOCAORI)))%>%
  filter(!is.na(MOCATOTAL)) %>%
  unique()

DF <- merge(DF, MOCA, all = TRUE)
selitems$MOCA <- c("TRAILS", "CUBE", "MOCACLOCK", "MOCANAM", "MOCADIG", "MOCALET", "MOCASERIAL",
                   "MOCAREP", "MOCAFLUEN", "MOCAABS", "MOCADLREC", "MOCAORI")
seltotals$MOCA <- c("MOCATOTAL")
# rm(MOCA)

# Financial	Capacity	Instrument	Short-Form	(FCI-SF) for ADNI3 
fci_itnames <- c("NICKELS", "QUARTERS", # MC_SCORE : Mental Calculation
               "BUDGET", "INSURANCE", "TAXCREDA", "TAXCREDB", "PAYEE", # FC_SCORE : Financial Conceptual Knowledge
               paste("SNGLCHK", 7:16, sep = ""), # SNGLCHKSCORE : Single Checkbook Register Task
               paste("COMPLXCHK", 17:30, sep = ""), # COMPLXCHKSCORE : Complex Checkbook/Register Task
               paste("BANKSTAT", c("A", "B", 32:37), sep = ""), # BANKSTATSCORE : Bank Statement Management
               "FCISCORE")

FCI <- fci %>%
  filter(!reg2del(.)) %>%
  filter(!(VISCODE %in% c("f", "uns1", "y1", "y2"))) %>%
  filter(DONE == "Yes") %>%
  dplyr::select(1:5, fci_itnames) %>%
  mutate_at(fci_itnames[1:39], function(x) ifelse(x == "OT", 0, as.numeric(x))) %>%
  mutate(TAXCRED = rowSums(select(., TAXCREDA, TAXCREDB), na.rm = TRUE),
         BANKSTAT = rowSums(select(., BANKSTATA, BANKSTATB), na.rm = TRUE)) %>%
  mutate(MC_SCORE = (NICKELS + QUARTERS),
         FC_SCORE = (BUDGET + INSURANCE + TAXCRED + PAYEE),
         SNGLCHKSCORE = rowSums(select(., paste("SNGLCHK", 7:16, sep = ""))),
         COMPLXCHKSCORE = rowSums(select(., paste("COMPLXCHK", 17:30, sep = ""))),
         BANKSTATSCORE = rowSums(select(., paste("BANKSTAT", 32:37, sep = ""), BANKSTAT)),
         FCITOTAL = MC_SCORE + FC_SCORE + SNGLCHKSCORE + COMPLXCHKSCORE + BANKSTATSCORE) %>%
  filter(!is.na(FCITOTAL)) %>%
  unique()

DF <- merge(DF, FCI, all = TRUE)
selitems$FCI <- c("MC_SCORE", "FC_SCORE", "SNGLCHKSCORE", "COMPLXCHKSCORE", "BANKSTATSCORE")
seltotals$FCI <- c("FCITOTAL")
# rm(FCI)

rm(ADAS, ANART, AVLT, BNT, CATFL,CDR, CLOCK, FAQ, FCI, LM, MINT, MMSE, MOCA, TMT, BNTMINT, 
   clock_itnames, fci_itnames, mmse_itnames, mmse_itwo, moca_itnames,
   reg2del, faq2number)

save(DF, selitems, seltotals, file = "~/ADNI_meta_oct2018/processed_data/neuropsycho_seltests.RData")

rm(DF, selitems, seltotals)

# Cogstate	Brief	Battery	(CBB)	Computerized	Testing for ADNI3 -> No viscode
# Everyday	Cognition (Self report and partner) from ADNI 3 -> compared to 10 years :/

