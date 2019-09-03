#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
library(GGally)
library(e1071) 
library(sinaplot)
require(cowplot)

#from https://stackoverflow.com/questions/47496364/rounding-digits-in-ggpairs
mycor <- function(data, mapping, alignPercent = 0.6, method = "pearson", 
    use = "complete.obs", corAlignPercent = NULL, corMethod = NULL, 
    corUse = NULL, sgnf=2, ...) {
    if (!is.null(corAlignPercent)) {
        stop("'corAlignPercent' is deprecated.  Please use argument 'alignPercent'")
    }
    if (!is.null(corMethod)) {
        stop("'corMethod' is deprecated.  Please use argument 'method'")
    }
    if (!is.null(corUse)) {
        stop("'corUse' is deprecated.  Please use argument 'use'")
    }
    useOptions <- c("all.obs", "complete.obs", "pairwise.complete.obs", 
        "everything", "na.or.complete")
    use <- pmatch(use, useOptions)
    if (is.na(use)) {
        warning("correlation 'use' not found.  Using default value of 'all.obs'")
        use <- useOptions[1]
    } else {
        use <- useOptions[use]
    }
    cor_fn <- function(x, y) {
        cor(x, y, method = method, use = use)
    }
    xCol <- deparse(mapping$x)
    yCol <- deparse(mapping$y)
    if (GGally:::is_date(data[[xCol]]) || GGally:::is_date(data[[yCol]])) {
        if (!identical(class(data), "data.frame")) {
            data <- fix_data(data)
        }
        for (col in c(xCol, yCol)) {
            if (GGally:::is_date(data[[col]])) {
                data[[col]] <- as.numeric(data[[col]])
            }
        }
    }
    if (is.numeric(GGally:::eval_data_col(data, mapping$colour))) {
        stop("ggally_cor: mapping color column must be categorical, not numeric")
    }
    colorCol <- deparse(mapping$colour)
    singleColorCol <- ifelse(is.null(colorCol), NULL, paste(colorCol, 
        collapse = ""))
    if (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete")) {
        if (length(colorCol) > 0) {
            if (singleColorCol %in% colnames(data)) {
                rows <- complete.cases(data[c(xCol, yCol, colorCol)])
            } else {
                rows <- complete.cases(data[c(xCol, yCol)])
            }
        } else {
            rows <- complete.cases(data[c(xCol, yCol)])
        }
        if (any(!rows)) {
            total <- sum(!rows)
            if (total > 1) {
                warning("Removed ", total, " rows containing missing values")
            } else if (total == 1) {
                warning("Removing 1 row that contained a missing value")
            }
        }
        data <- data[rows, ]
    }
    xVal <- data[[xCol]]
    yVal <- data[[yCol]]
    if (length(names(mapping)) > 0) {
        for (i in length(names(mapping)):1) {
            tmp_map_val <- deparse(mapping[names(mapping)[i]][[1]])
            if (tmp_map_val[length(tmp_map_val)] %in% colnames(data)) 
                mapping[[names(mapping)[i]]] <- NULL
            if (length(names(mapping)) < 1) {
                mapping <- NULL
                break
            }
        }
    }
    if (length(colorCol) < 1) {
        colorCol <- "ggally_NO_EXIST"
    }
    if ((singleColorCol != "ggally_NO_EXIST") && (singleColorCol %in% 
        colnames(data))) {
        cord <- plyr::ddply(data, c(colorCol), function(x) {
            cor_fn(x[[xCol]], x[[yCol]])
        })
        colnames(cord)[2] <- "ggally_cor"
        cord$ggally_cor <- signif(as.numeric(cord$ggally_cor), 
            sgnf)
        lev <- levels(data[[colorCol]])
        ord <- rep(-1, nrow(cord))
        for (i in 1:nrow(cord)) {
            for (j in seq_along(lev)) {
                if (identical(as.character(cord[i, colorCol]), 
                  as.character(lev[j]))) {
                  ord[i] <- j
                }
            }
        }
        cord <- cord[order(ord[ord >= 0]), ]
        cord$label <- GGally:::str_c(cord[[colorCol]], ": ", cord$ggally_cor)
        xmin <- min(xVal, na.rm = TRUE)
        xmax <- max(xVal, na.rm = TRUE)
        xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
            (xmax - xmin))
        ymin <- min(yVal, na.rm = TRUE)
        ymax <- max(yVal, na.rm = TRUE)
        yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
            (ymax - ymin))
        p <- ggally_text(label = GGally:::str_c("Corr: ", signif(cor_fn(xVal, 
            yVal), sgnf)), mapping = mapping, xP = 0.5, yP = 0.9, 
            xrange = xrange, yrange = yrange, color = "black", 
            ...) + theme(legend.position = "none")
        xPos <- rep(alignPercent, nrow(cord)) * diff(xrange) + 
            min(xrange, na.rm = TRUE)
        yPos <- seq(from = 0.9, to = 0.2, length.out = nrow(cord) + 
            1)
        yPos <- yPos * diff(yrange) + min(yrange, na.rm = TRUE)
        yPos <- yPos[-1]
        cordf <- data.frame(xPos = xPos, yPos = yPos, labelp = cord$label)
        cordf$labelp <- factor(cordf$labelp, levels = cordf$labelp)
        p <- p + geom_text(data = cordf, aes(x = xPos, y = yPos, 
            label = labelp, color = labelp), hjust = 1, ...)
        p
    }  else {
        xmin <- min(xVal, na.rm = TRUE)
        xmax <- max(xVal, na.rm = TRUE)
        xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * 
            (xmax - xmin))
        ymin <- min(yVal, na.rm = TRUE)
        ymax <- max(yVal, na.rm = TRUE)
        yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * 
            (ymax - ymin))
        p <- ggally_text(label = paste("Corr:\n", signif(cor_fn(xVal, 
            yVal), sgnf), sep = "", collapse = ""), mapping, xP = 0.5, 
            yP = 0.5, xrange = xrange, yrange = yrange, ...) + 
            theme(legend.position = "none")
        p
    }
}

plot_pref_sina <- function(table_pref){
  png('graph/Lib_pref_sina.png', res=600, height=2400, width=3000)
  sinaplot(pref~strain, data=table_pref, pch=20, cex=0.3)
  dev.off()
  }

plot_pairs <- function(table_pref_cast){
  textsize <- 7
  p  <- ggpairs(select(table_pref_cast,-ID),
		lower=list(continuous=wrap(ggally_points,size=0.2,color="black",alpha=0.3)),
		upper = list(continuous=wrap(mycor,sgnf=2,size=4,color="black")),
		diag = list(continuous = wrap("densityDiag",fill="grey"))) +
	  theme(plot.title=element_text(size=textsize,face="bold"),
		panel.grid.major = element_blank(),
		axis.title=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"))
  ggsave('graph/LibCorPairs.png',p,width=6.5,height=6.5)
  }

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Bris07','NDako16')
table_pref <- read_tsv('result/data_pref.tsv') %>%
                mutate(strain=factor(strain,levels=StrainLevels)) %>%
                filter(strain!='NA')
table_pref_cast <- cast(table_pref,ID~strain,value='pref')
plot_pref_sina(table_pref)
plot_pairs(table_pref_cast)
