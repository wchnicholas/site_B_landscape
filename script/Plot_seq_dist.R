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

plot_dist_vs_cor <- function(data_table, y_label, graphname){
  textsize <- 7
  print (paste("Rank Cor:", cor(data_table$landscape_cor, data_table$identity, use="complete.obs", method='spearman'),sep=' '))
  print (cor.test(data_table$landscape_cor, data_table$identity, use="complete.obs", method='spearman'))
  p <- ggplot(data_table,aes(x=identity,y=landscape_cor)) +
         geom_point(size=1, color='black') +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold"),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_y_continuous(breaks=c(-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8),
                            labels=c('-0.4','-0.2', '0.0', '0.2', '0.4', '0.6', '0.8')) +
         ylab(bquote(bold("Correlation of\nsite B landscape"))) +
         xlab(y_label)
  ggsave(graphname,p,height=2,width=2)
  }

data_table <- read_tsv('result/EpiB_seq_dist.tsv')
HA_data_table  <- mutate(data_table, identity=HA_identity)
RBS_data_table <- mutate(data_table, identity=RBS_identity)
plot_dist_vs_cor(HA_data_table,  "Sequence identity of\nHA ectodomain (# of aa)", 'graph/Identity_vs_cor_HA.png')
plot_dist_vs_cor(RBS_data_table, "Sequence identity of\nRBS base (# of aa)", 'graph/Identity_vs_cor_RBS.png')
