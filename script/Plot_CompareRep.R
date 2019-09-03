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
  table_data <- filter(table_data,strain==strainname)
  print (strainname)
  R1fit = log10(table_data$Rep1Fitness)
  R2fit = log10(table_data$Rep2Fitness)
  R1fit[!is.finite(R1fit)] <- NA
  R2fit[!is.finite(R2fit)] <- NA
  print (paste("Pearson Cor:", cor(R1fit,R2fit,use="complete.obs"),sep=' '))
  p <- ggplot(table_data,aes(x=log10(Rep1Fitness),y=log10(Rep2Fitness))) +
	 geom_point(size=0.1) +
	 theme(plot.title=element_text(size=textsize,face="bold"),
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

table_data <- read_csv('result/data_all.csv') %>%
                filter(strain!='Bris07P194') %>%
                mutate(strain= replace(strain, which(strain=="Bris07L194"), 'Bris07'))
p1 <- PlotCompareFit_Rep(table_data,'HK68')
p2 <- PlotCompareFit_Rep(table_data,'Bk79')
p3 <- PlotCompareFit_Rep(table_data,'Bei89')
p4 <- PlotCompareFit_Rep(table_data,'Mos99')
p5 <- PlotCompareFit_Rep(table_data,'Bris07')
p6 <- PlotCompareFit_Rep(table_data,'NDako16')
p <- grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
ggsave('graph/Compare_Rep.png',p,height=6,width=4)
