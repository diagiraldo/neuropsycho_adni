# Pre-process functions

library(dplyr)

it <- c("Q1SCORE" , "Q2SCORE" , "Q3SCORE" , "Q4SCORE" , "Q5SCORE" , "Q6SCORE" , 
        "Q7SCORE" , "Q8SCORE" , "Q9SCORE" , "Q10SCORE" , "Q11SCORE" , "Q12SCORE" , 
        "Q13SCORE" , "MMORITIME" , "MMORISPACE" , "MMREGI" , "MMRECALL" , 
        "MMSPELLBKW" , "MMNAM" , "MMCOMMAND" , "MMREPEAT" , "MMREAD" , "MMWRITE" , 
        "MMDRAW" , "TRAILS" , "CUBE" , "MOCACLOCK" , "MOCANAM" , "MOCADIG" , 
        "MOCALET" , "MOCASERIAL" , "MOCAREP" , "MOCAFLUEN" , "MOCAABS" , "MOCADLREC" , 
        "MOCAORI" , "CLOCKSCOR" , "COPYSCOR" , "TRAASCOR" , "TRABSCOR" , 
        "LIMMTOTAL" , "LDELTOTAL" , "CATANIMSC" , "RAVLT.IMMED" , 
        "AVTOT6" , "AVTOTB" , "AVDEL30MIN" , "AVDELTOT" , "BMNOCUE" , "BMCUED")

# Reverse scores
reverse_scores <- function(DF, selitems){
    # Reverse MMSE items
    DF <- DF %>%
        mutate(MMORITIME = 5-MMORITIME, MMORISPACE = 5-MMORISPACE, MMREGI = 3-MMREGI, MMRECALL = 3-MMRECALL,
               MMSPELLBKW = 5-MMSPELLBKW, MMNAM = 2-MMNAM, MMCOMMAND = 3-MMCOMMAND,
               MMREPEAT = 1-MMREPEAT, MMREAD = 1-MMREAD, MMWRITE = 1-MMWRITE, MMDRAW = 1-MMDRAW)
    
    # Reverse Logical Memory 25-x
    DF <- DF %>%
        mutate(LIMMTOTAL = 25-LIMMTOTAL, LDELTOTAL = 25-LDELTOTAL)
    
    # Reverse Rey Auditory Verbal Learning Test and Delayed version 15-x
    DF <- DF %>%
        mutate_at(selitems$AVLT[2:5], function(x) 15-x) %>%
        mutate(RAVLT.IMMED = 75-RAVLT.IMMED)
    
    # Reverse Clock Drawing - Copying
    DF <- DF %>%
        mutate_at(selitems$CLOCK, function(x) 5-x)
    
    # Reverse Category fluency 
    DF <- DF %>%
        mutate(CATANIMSC = 60-CATANIMSC)
    
    # Reverse BNTMINT 30-x
    DF <- DF %>%
        mutate_at(selitems$BNTMINT, function(x) 30-x) 
    
    # Reverse MoCA 
    DF <- DF %>%
        mutate(TRAILS = 1-TRAILS, CUBE = 1-CUBE, MOCACLOCK = 3-MOCACLOCK, MOCANAM = 3-MOCANAM,
               MOCADIG = 2-MOCADIG, MOCALET = 1-MOCALET, MOCASERIAL = 3-MOCASERIAL, MOCAREP = 2-MOCAREP,
               MOCAFLUEN = 1-MOCAFLUEN, MOCAABS = 2-MOCAABS, MOCADLREC = 5-MOCADLREC, MOCAORI = 6-MOCAORI)
    
    return(DF)
}

# Calculate SRB Z-scores
calculate_srbz <- function(A, it, coef){
    A_std <- data.frame(matrix(ncol = length(it), nrow = nrow(A)))
    names(A_std) <- it
    tmpregre <- A %>% 
        mutate(intercept = 1) %>%
        dplyr::select(c("intercept", "PTEDUCAT", "AGE")) %>%
        as.matrix()
    for (s in seq_along(it)){
        res <- A[, it[s]] - tmpregre %*% t(as.matrix(coef[s, 2:4]))
        A_std[, it[s]] <- res/coef$sigma[s]
        attr(A_std[, it[s]], "dimnames") <- NULL
    }
    A_std <- cbind(select(A, -all_of(it)), A_std)
    return(A_std)
}

# Calculate composite scores 
composite_scores <- function(newdata, weights, centerit = 0){
    itorder <- colnames(weights)
    facorder <- rownames(weights)
    dmat <- as.matrix(scale(newdata[ ,itorder], center = centerit[itorder], scale = FALSE))
    S <- dmat %*% t(weights)
    S <- as.data.frame(S[, rev(facorder)])
    return(S)
}

# Calculate distances to MCI subgroups medoids
distance2meds <- function(S, factorcov, grmeds){
    retdf <- select(S, -all_of(namesfac)) %>%
        mutate(GR = NA)
    tmpdist <- as.data.frame(matrix(NA, nrow = nrow(S), ncol = nrow(grmeds)))
    names(tmpdist) <- gsub(" ", "_", paste("dist", grmeds$GR))
    namesfac <- colnames(factorcov)
    for(i in 1:nrow(S)){
        for(j in 1:nrow(grmeds)){
            x <- as.matrix(S[i, namesfac] - grmeds[j, namesfac])
            tmpdist[i, j] <- sqrt(x %*% factorcov %*% t(x))
        }
        retdf$GR[i] <- grmeds$GR[which.min(tmpdist[i, ])]
    }
    return(cbind(retdf, tmpdist))
}

# Predict progression to dementia with a list of pre-trained random forest
predict_MCIprogression <- function(S, RFlist){
    retdf <- select(S, -all_of(namesfac))
    for (i in 1:length(RFlist)){
        tmppred <- as.data.frame(predict(RFlist[[i]], S, type = "prob")) %>%
            mutate(prediction = predict(RFlist[[i]], S))
        names(tmppred) <- paste(gsub(" ", "", names(RFlist)[i]), names(tmppred), sep = "_")
        tmppred <- select(tmppred, ends_with("prediction"), ends_with("converterMCI"))
        retdf <- cbind(retdf, tmppred)
    }
    return(retdf)
}

