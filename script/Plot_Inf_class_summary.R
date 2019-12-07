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
require(cowplot)
require(reshape2)

plot_pie <- function(table_summary, main, model){
  textsize <- 7
  colorscale  <- c(brewer.pal(8,"Set3"))[1:4]
  colorscale  <- c('purple','blue','red','white')
  p <-  ggplot(table_summary,aes(main, freq, fill=epi_class)) +
          geom_bar(stat="identity",width=1,color='black') +
          theme_cowplot(12) +
          theme(axis.title=element_text(size=textsize,face="bold"),
                axis.text=element_text(size=textsize,face="bold"),
                axis.text.x=element_blank(),
                legend.key.size=unit(0.2,'in'),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='right') +
          labs(y=expression(bold('# of residues')),x=expression()) +
          scale_fill_manual(values=colorscale,drop=FALSE) +
          coord_polar("y", start=0)
  ggsave(paste('graph/epi_class_summary_',model,'_',main,'.png',sep=''),p,height=2.2,width=4)
  }

plot_epi_summary <- function(epi_table, model){
  epi_class_levels <- c('mix','pos','neg','none')
  single_epi_table <- epi_table %>%
                        filter(mut_class=='single') %>%
                        select(epi_class) %>%
                        table() %>%
                        data.table() %>%
                        mutate(freq=N/sum(N))
  single_epi_table <- rename(single_epi_table,epi_class=.)
  double_epi_table <- epi_table %>%
                        filter(mut_class=='double') %>%
                        select(epi_class) %>%
                        table() %>%
                        data.table() %>%
                        mutate(freq=N/sum(N))
  double_epi_table <- rename(double_epi_table,epi_class=.)
  print (single_epi_table)
  print (double_epi_table)
  single_epi_table <- single_epi_table %>%
                        select(epi_class,freq) %>%
                        mutate(epi_class=factor(epi_class,levels=epi_class_levels))
  double_epi_table <- double_epi_table %>%
                        select(epi_class,freq) %>%
                        mutate(epi_class=factor(epi_class,levels=epi_class_levels))
  plot_pie(single_epi_table, 'single', model) 
  plot_pie(double_epi_table, 'double', model) 
  }

epi_classification <- function(CI_pos, CI_neg){
  if (CI_pos > 0 & CI_neg > 0){return ('mix')}
  else if (CI_pos > 0 & CI_neg==0){return ('pos')}
  else if (CI_neg > 0 & CI_pos==0){return ('neg')}
  else if (CI_pos==0 & CI_neg==0){return ('none')}
  else {print ('something is wrong with epi_classification scheme')}
  }

mut_classification <- function(mut1,mut2){
  if (mut1==mut2){return ('single')}
  else {return ('double')}
  }

wrapper <- function(model){
  epi_table  <- read_tsv(paste('result/Inf_heatmap_overall_',model,'.tsv',sep='')) %>%
		  mutate(epi_class=mapply(epi_classification,CI_pos,CI_neg)) %>%
		  mutate(mut_class=mapply(mut_classification, mut1, mut2))
  plot_epi_summary(epi_table, model)
  }

wrapper('linearity')
wrapper('nonlinearity')
