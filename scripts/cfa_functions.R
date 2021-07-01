# Functions to use in CFA
# Jun 2021

# it requires lavaan, dplyr

# Extract sub-scores in model
its_in_model <- function(model){
    x <- lavaanify(model)
    return(unique(x$rhs[x$user == 1]))
}

# Extract factor loadings, factor covariance, and calculate weights from CFA result
loadings_weights_fcov <- function(fitcfa){
    L <- lavInspect(fitcfa, what = "est")$lambda
    Cz <- lavInspect(fitcfa, what = "est")$psi
    Ce <- lavInspect(fitcfa, what = "est")$theta
    Ceinv <- diag(1/diag(Ce))
    tmp <- t(L)%*%(Ceinv)
    W <- solve(tmp%*%L)%*%tmp
    colnames(W) <- rownames(L)
    W <- W[rev(1:nrow(W)),]
    W <- W[,rev(1:ncol(W))]
    return(list(L = L, W = W, Cz = Cz))
}

# Save weights in csv
save_weights_csv <- function(W, filedir){
    W1 <- as.data.frame(round(t(W), digits = 3)) %>%
        mutate(subscore = colnames(W)) %>%
        select(c(subscore, rev(rownames(W))))
    W1 <- W1[nrow(W1):1,]
    write.table(W1, file = filedir, sep = ",", row.names = FALSE)
    cat("Table with weights saved in", filedir, "\n")
}

# Calculate composite scores for new data
composite_scores <- function(newdata, weights, centerit = 0){
    itorder <- colnames(weights)
    facorder <- rownames(weights)
    dmat <- as.matrix(scale(newdata[ ,itorder], center = centerit[itorder], scale = FALSE))
    S <- dmat %*% t(weights)
    S <- as.data.frame(S[, rev(facorder)])
    return(S)
}

# Pairwise comparisons od domain composite scores
compare_domain_scores <- function(S, vars, grvar){
    nvars <- length(vars)
    comps <- combn(as.character(unique(S[,grvar])), m = 2)
    n_tests <- nvars*ncol(comps)
    grcomp <- data.frame(variable = character(n_tests), 
                         groups = character(n_tests),
                         pval = numeric(n_tests), 
                         estimatediff = numeric(n_tests),
                         stringsAsFactors = FALSE)
    for(i in 1:ncol(comps)){
        grcomp$groups[(i-1)*nvars + (1:nvars)] <- paste(comps[, 1], collapse = " vs. ")
        tmpS <- filter(S, get(grvar) %in% comps[, i])
        for (v in 1:nvars){
            fmla <- as.formula(paste(vars[v], grvar, sep = " ~ "))
            tmptest <- wilcox.test(fmla, data = tmpS, conf.int = TRUE)
            grcomp$variable[(i-1)*nvars + v] <- vars[v]
            grcomp$pval[(i-1)*nvars + v] <- tmptest$p.value
            grcomp$estimatediff[(i-1)*nvars + v] <- tmptest$estimate
        }
    }
    grcomp <- mutate(grcomp,
                     pval.bonferroni = pval*n_tests)
    return(grcomp)
}

