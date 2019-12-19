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

plot_KLdist_vs_pref <- function(data_table, strain, color){
  graphname <- paste('graph/KLdist_vs_pref',strain,'.png',sep='')
  textsize <- 7
  scaleFUN <- function(x) sprintf("%.1f", x)
  data_table <- data_table %>%
                  filter(Background==strain)
  print (paste("correlation in ", strain, ": ", cor(data_table$KL, data_table$MeanPref, method = "spearman"), sep=''))
  p <-  ggplot(data_table, aes(x=KL,y=MeanPref)) +
	  geom_point(size=0.5, color=color) +
          geom_errorbar(aes(ymin=MeanPref-StdPref, ymax=MeanPref+StdPref), color=color, width=0) +
          geom_errorbarh(aes(xmin=Lowerbound, xmax=Upperbound), color=color, height=0) +
          theme_cowplot(12) + 
	  theme(plot.title=element_text(size=textsize,face="bold"),
		axis.title=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		legend.title=element_blank(),
		legend.text=element_text(size=textsize,face="bold"),
		legend.position='right') +
	  ylab(bquote(bold(Mean~preference))) +
	  xlab(bquote(bold(KL~distance))) +
          scale_y_continuous(labels=scaleFUN) +
          scale_x_continuous(breaks=c(0.0,0.2,0.4), labels=c('0.0','0.2','0.4'))
  ggsave(graphname, p, height=1.4, width=1.4)
  }

strain_name_converter <- function(strain){
  if (strain=='Bris07L194'){return ('Bris07')}
  else {return (strain)}
  }

strain_levels <- c('HK68','Bk79','Bei89','Mos99','Bris07','NDako16')
YearPrefTable <- read_tsv('result/Prefs_ByYear.tsv') %>%
                   filter(Background!='Bris07P194') %>%
                   mutate(Background=factor(Background,levels=strain_levels))
KLdist_table <- read_csv('result/AllDKLInfo_2018.csv')
KLdist_allbutBstar <- KLdist_table %>%
                        filter(Region=='OtherEpibutB*') %>%
                        rename(Year=year) %>%
                        rename(KL=Dkl) %>%
                        rename(Background=Strain) %>%
                        mutate(Background=mapply(strain_name_converter, Background)) %>%
                        mutate(Background=factor(Background,levels=strain_levels)) %>%
                        select(-Region)
data_table <- inner_join(KLdist_allbutBstar, YearPrefTable)
palette  <- c(brewer.pal(9,"Set3"))
palette  <- c(palette[3:12])
for (n in 1:length(strain_levels)){
  strain <- strain_levels[n]
  color <- palette[n]
  plot_KLdist_vs_pref(data_table, strain, color)
  }
