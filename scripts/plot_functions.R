######################################
### Plot functions - using ggplot2 ###
######################################
### D. Giraldo, Oct 2019

# Function to plot correlation matrix, requires reshape2 and ggplot2
corr.plot <- function(COR, cut = TRUE, corrange = c(-1,1)){
    if (cut) {
        COR[upper.tri(COR, diag = FALSE)] <- NA
    }
    m <- melt(COR)
    p <- ggplot(m, aes(Var1, Var2)) + geom_tile(aes(fill = value), colour = "white") + 
        theme_bw() + labs(x = "",y = "") +
        scale_x_discrete(expand = c(0, 0), position = "top") + 
        scale_y_discrete(expand = c(0, 0), limits = rev(levels(m$Var2))) +
        scale_fill_distiller(limits = corrange, na.value = "white", palette = "Spectral", name = "Correlation value ") +
        theme(axis.text.y = element_text(size = 6, vjust = 0.5, colour = "grey30", hjust = 1), 
              axis.text.x = element_text(angle = 90, size = 6, colour = "grey30", hjust = 0),
              axis.ticks = element_blank(), axis.text.x.top = element_text(vjust = 0.5),
              plot.background = element_rect(fill = "transparent", colour = NA), plot.margin = unit(c(-0.45,0,0,-0.45), "cm"), 
              legend.background = element_rect(fill = "transparent",colour = NA),
              legend.justification = c(0, 0), legend.direction = "horizontal",
              legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = ggplot2::margin(0,0,0,0)) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.75, title.position = "top", title.hjust = 0.5))
    if (cut) {
        p <- p + theme(legend.position = c(0.1, 0.1))
    } else {
        p <- p + theme(legend.position = "bottom", legend.justification = c(0.5, 0))
    }
    return(p)
}

# Function to plot weight matrix, requires reshape2 and ggplot2
weight.plot <- function(W, nam = c("Var1", "Var2", "Weight"), valrange = c(-1,1), flip = FALSE, withvalues = FALSE){
    m <- melt(W)
    p <- ggplot(m, aes(Var1, Var2)) + geom_tile(aes(fill = value), colour = "white") +
        theme_bw() + labs(x = nam[1], y = nam[2]) +
        scale_x_discrete(expand = c(0, 0), limits = rev(levels(m$Var1))) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradient2(limits = valrange, na.value = "white", name = nam[3], 
                             guide = guide_colorbar(title.hjust = 0, title.vjust = 1, label.theme = element_text(size = 7), 
                                                    barwidth = 10, barheight = 0.5, title.theme = element_text(size = 8))) +
        theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, colour = "grey30", vjust = 0.5),
              axis.text.y = element_text(size = 6, vjust = 0.5, colour = "grey30", hjust = 1),
              legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
              legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
              plot.background = element_rect(fill = "transparent",colour = NA), plot.margin = unit(c(0,0,0,0), "cm")) +
        guides(fill = guide_colorbar(barwidth = 8, barheight = 0.75, title.position = "top", title.hjust = 0.5))
    if (flip) {
        p <- p + coord_flip()
    }
    if (withvalues) {
        p <- p + geom_text(aes(Var1, Var2, label = ifelse(value > 0, round(value, digits = 2), "")), color = "grey30", size = 3)
    }
    return(p)
}

# Function to plot Scores boxplots
scores.boxplots <- function(S, namesfac, nam = c("variable", "value", "group"), flip = FALSE, by.diagnosis = TRUE){
    m <- melt(S, measure.vars = namesfac)
    if (by.diagnosis){
        m <- mutate(m, GR = DIAGNOSIS)
    }
    p <- ggplot(m, aes(variable, value, colour = GR)) + geom_boxplot(outlier.size = 0.5) +
        theme_bw() + labs(x = nam[1], y = nam[2], colour = nam[3]) +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
              legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
              plot.margin = unit(c(0,0,0,0), "cm"))
    if (flip) {
        p <- p + coord_flip() + scale_x_discrete(limits = rev(levels(m$variable)))
    }
    return(p)
}

# Function to plot Scores boxplots
scores.boxplots2 <- function(S, namesfac, nam = c("variable", "value", "group"), flip = FALSE, by.diagnosis = TRUE){
    m <- melt(S, measure.vars = namesfac)
    if (by.diagnosis){
        m <- mutate(m, GR = DIAGNOSIS)
    }
    p <- ggplot(m, aes(GR, value, colour = variable)) + geom_boxplot(outlier.size = 0.5) +
        theme_bw() + labs(x = nam[3], y = nam[2], colour = nam[1]) +
        theme(plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "top", legend.background = element_rect(fill = "transparent",colour = NA),
              legend.margin = ggplot2::margin(0,0,0,0), legend.box.margin = unit(c(0.2,0,0,0), "cm"),
              plot.margin = unit(c(0,0,0,0), "cm"))
    if (flip) {
        p <- p + coord_flip() + scale_x_discrete(limits = rev(levels(m$variable)))
    }
    return(p)
}
