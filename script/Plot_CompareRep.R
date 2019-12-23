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


PlotCompareFit_Rep <- function(table_data,strainname){
  textsize=7
  colorscale <- c(brewer.pal(9,"Set1"))
  table_data <- filter(table_data,strain==strainname)
  print (strainname)
  R1fit = log(table_data$Rep1Fitness)
  R2fit = log(table_data$Rep2Fitness)
  R1fit[!is.finite(R1fit)] <- NA
  R2fit[!is.finite(R2fit)] <- NA
  print (paste("Pearson Cor:", cor(R1fit,R2fit,use="complete.obs"),sep=' '))
  p <- ggplot(table_data,aes(x=log(Rep1Fitness),y=log(Rep2Fitness),color=mut_type)) +
	 geom_point(size=0.3) +
         theme_cowplot(12) +
         scale_color_manual(values=c('gray', colorscale)) +
	 theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
	       axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
	 xlab(bquote(bold('Replicate 1'))) +
	 ylab(bquote(bold('Replicate 2'))) +
         ggtitle(strainname)
  return (p)
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

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Bris07','NDako16')
table_data <- read_csv('result/data_all.csv') %>%
                filter(strain!='Bris07P194') %>%
                mutate(strain= replace(strain, which(strain=="Bris07L194"), 'Bris07')) %>%
                mutate(mut_type=mapply(coloring, ID)) %>%
                mutate(mut_type=factor(mut_type, levels=c('mut', StrainLevels))) %>%
                arrange(mut_type)
p1 <- PlotCompareFit_Rep(table_data,'HK68')
p2 <- PlotCompareFit_Rep(table_data,'Bk79')
p3 <- PlotCompareFit_Rep(table_data,'Bei89')
p4 <- PlotCompareFit_Rep(table_data,'Mos99')
p5 <- PlotCompareFit_Rep(table_data,'Bris07')
p6 <- PlotCompareFit_Rep(table_data,'NDako16')
p <- grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
ggsave('graph/Compare_Rep.png',p,height=6,width=4)
