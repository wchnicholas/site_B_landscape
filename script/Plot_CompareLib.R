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
library(ggforce)
require(cowplot)

plot_pref_sina <- function(table_pref, StrainLevels){ 
  textsize <- 7
  table_pref <- table_pref %>%
                  mutate(strain=factor(strain, levels=rev(StrainLevels)))
  colorscale <- c(brewer.pal(9,"Set1"))
  p <- ggplot(table_pref, aes(x=strain, y=pref, color=mut_type, group=strain)) +
         geom_sina(size=0.2) +
         theme_cowplot(12) +
         scale_color_manual(values=c('gray', colorscale)) +
         theme(plot.title=element_text(size=textsize,face="bold"),
               panel.grid.major = element_blank(),
               legend.position = "none", 
               legend.text=element_text(size=textsize,face="bold"),
               legend.key.size = unit(0.5,"line"), 
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold")) +
         coord_flip() +
         ylab(bquote(bold("Preference")))
  ggsave('graph/Lib_pref_sina.png',p,width=2,height=3.5)
  }

plot_pairs <- function(table_pref_cast){
  textsize <- 7
  colorscale <- c(brewer.pal(9,"Set1"))
  table_pref_cast <- table_pref_cast 
  p  <- ggpairs(data=table_pref_cast, columns=3:8,
                mapping=ggplot2::aes(colour=mut_type),
		lower =list(continuous=wrap(ggally_points,size=0.1,alpha=1)),
		upper = list(continuous="blank"),
		diag  = list(continuous ="blank")) +
          theme_cowplot(12) +
	  theme(plot.title=element_text(size=textsize,face="bold"),
		panel.grid.major = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
		axis.title=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"))
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c('gray', colorscale)) +
        scale_color_manual(values=c('gray', colorscale))
      }
    }
  ggsave('graph/LibCorPairs.png',p,width=4,height=4)
  }

coloring <- function(ID){
  if (ID=='KGSESV'){return ('HK68')}
  else if (ID=='EESENV'){return ('Bk79')}
  else if (ID=='EEYENV'){return ('Bei89')}
  else if (ID=='QKYDST'){return ('Mos99')}
  else if (ID=='HKFDFA'){return ('Bris07')}
  else if (ID=='HNSDFA'){return ('NDako16')}
  else {return ('mut')}
  }

sizing <- function(mut_type){
  if (mut_type=='mut'){return (1)}
  else {return (3)}
  }

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Bris07','NDako16')
table_pref <- read_tsv('result/data_pref.tsv') %>%
                mutate(strain=factor(strain,levels=StrainLevels)) %>%
                filter(strain!='NA') %>%
                mutate(mut_type=mapply(coloring, ID)) %>%
                mutate(size=mapply(sizing, mut_type)) %>%
                mutate(mut_type=factor(mut_type, levels=c('mut', StrainLevels))) %>%
                arrange(mut_type)
table_pref_cast <- cast(table_pref,ID+mut_type~strain,value='pref') %>%
                     arrange(mut_type)
plot_pref_sina(table_pref, StrainLevels)
plot_pairs(table_pref_cast)
