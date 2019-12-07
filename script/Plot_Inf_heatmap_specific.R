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

convert_label_to_numeric <- function(resi, mut_factor){
  i <- unique(as.numeric(mut_factor[mut_factor==resi]))
  return (i)
  }

plot_heatmap <- function(epi_table,WTresibox,strainname){
  textsize <- 9
  epi_table_class  <- epi_table %>%
                        filter(strain==strainname)
  WTresibox_strain <- WTresibox %>%
                        filter(strain==strainname)
  p  <- ggplot()+
          geom_point(data=epi_table_class,aes(x=mut1,y=mut2,fill=class), color='black', pch=21) + 
          scale_fill_gradientn(colours=c('red',"white",'blue'),
                limits=c(-1,1),
                values=rescale(c(0,0.5,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,2.4)) +
          theme_classic() +
          theme(plot.title=element_text(size=textsize,face="bold"),
                text = element_text(size=textsize,face="bold"),
                legend.key.size = unit(0.5, 'lines'),
                legend.position='none',
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.title=element_blank(),
                axis.text.x.top=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          scale_x_discrete(position = "top") +
	  xlab("") + 
	  ylab("") +
          ggtitle(strainname) +
          geom_rect(data=WTresibox_strain, size=0.3, fill=NA, colour="black",
                    aes(xmin=x-0.5, xmax=x+0.5, ymin=y-0.5, ymax=y+0.5))
  return (p)
  }

wrapper <- function(model){
  epi_table  <- read_tsv(paste('result/Inf_heatmap_specific_',model,'.tsv',sep=''))
  WTresibox  <- read_tsv('data/WTheatmap.tsv')
  mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
  epi_table  <- epi_table %>%
		  mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		  mutate(mut2=factor(mut2,levels=mut_levels))
  WTresibox  <- WTresibox %>%
		  mutate(x=as.numeric(sapply(resi1, convert_label_to_numeric, mut_factor=epi_table$mut1))) %>%
		  mutate(y=as.numeric(sapply(resi2, convert_label_to_numeric, mut_factor=epi_table$mut2)))
  p1 <- plot_heatmap(epi_table,WTresibox,'HK68')
  p2 <- plot_heatmap(epi_table,WTresibox,'Bk79')
  p3 <- plot_heatmap(epi_table,WTresibox,'Bei89')
  p4 <- plot_heatmap(epi_table,WTresibox,'Mos99')
  p5 <- plot_heatmap(epi_table,WTresibox,'Bris07')
  p6 <- plot_heatmap(epi_table,WTresibox,'NDako16')
  p <- grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
  ggsave(paste('graph/Inf_heatmap_specific_',model,'.png',sep=''),p,height=8,width=5)
  }

wrapper('linearity')
wrapper('nonlinearity')
