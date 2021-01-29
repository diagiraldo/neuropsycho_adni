###############################################
### Auxiliar functions to process ADNI data ###
###############################################
### D. Giraldo, Oct 2019

# Function to detect double registers RID-VISCODE
reg2del <- function(DF){
    tmp <- with(DF, table(RID, VISCODE))
    idx <- which(tmp > 1, arr.ind = TRUE)
    nam <- dimnames(tmp[idx[,1],idx[,2]])
    l <- paste(DF$RID, DF$VISCODE, sep=":") %in% paste(nam$RID, nam$VISCODE, sep=":")
    return(l)
}

# Function to detect first visit in data.frame, requires dplyr
firstvisit <- function(DF){
    tmp1 <- DF %>%
        dplyr::select(RID, VISCODE) %>%
        mutate(VISMONTH = as.numeric(ifelse(VISCODE == "bl", 0, gsub("m", "", VISCODE)))) 
    tmp2 <- data.frame(RID = unique(DF$RID), FIRSTVIS = NA) 
    tmp2$FIRSTVIS <- sapply(tmp2$RID, function(x) min(tmp1$VISMONTH[tmp1$RID == x]))
    tmp1 <- merge(tmp1, tmp2, all.x = TRUE) %>%
        mutate(is.firstvis = (VISMONTH == FIRSTVIS))
    return(tmp1)
}
