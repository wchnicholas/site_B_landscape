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
require(cowplot)
require(reshape2)

plot_heatmap <- function(epi_table,model,class){
  textsize <- 7
  if (class=='pos'){
    epi_table_class <- select(epi_table,mut1,mut2,CI_pos) %>% rename(epi_count=CI_pos)
    max_color <- "blue"
    }
  if (class=='neg'){
    epi_table_class <- select(epi_table,mut1,mut2,CI_neg) %>% rename(epi_count=CI_neg)
    max_color <- "red"
    }
  p  <- ggplot(epi_table_class, aes(x=mut1,y=mut2))+
          geom_point(aes(fill=epi_count,size=epi_count),color='black',pch=21) + 
          scale_fill_gradientn(colours=c("white", max_color),
                limits=c(0,6),
                values=rescale(c(0,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,2.4)) +
          theme_classic() +
          theme(#panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                text = element_text(size=textsize,face="bold"),
                legend.key.size = unit(0.5, 'lines'),
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.title=element_blank(),
                axis.text.x.top=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          scale_x_discrete(position = "top") +
	  xlab("") + 
	  ylab("")
  ggsave(paste('graph/Inf_heatmap_overall_',model,'_',class,'.png',sep=''),p,height=2.0,width=2.8)
  }

wrapper <- function(model){
  epi_table  <- read_tsv(paste('result/Inf_heatmap_overall_',model,'.tsv', sep=''))
  mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
  epi_table  <- epi_table %>%
		  mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		  mutate(mut2=factor(mut2,levels=mut_levels))
  plot_heatmap(epi_table,model,'pos')
  plot_heatmap(epi_table,model,'neg')
  }

wrapper('linearity')
wrapper('nonlinearity')
