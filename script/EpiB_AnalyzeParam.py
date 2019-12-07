#!/usr/bin/python
import os
import sys
import numpy as np
from math import exp, log
from operator import itemgetter
from collections import defaultdict

def reading_epi_file(epi_file):
  epi_dict = defaultdict(dict)
  ID_list  = []
  strains  = []
  infile = open(epi_file,'r')
  for line in infile.xreadlines():
    if 'strain' in line: continue
    line = line.replace('"','').rstrip().rsplit(',')
    mut1   = line[0]
    mut2   = line[1]
    strain = line[2]
    strain = 'Bris07' if strain == 'Bris07L194' else strain
    CI_low = -float(line[3])
    CI_up  = -float(line[4])
    e      = -float(line[5])
    ID = mut1+'-'+mut2
    epi_dict[strain][ID] = {'e':e, 'CI_up':CI_up, 'CI_low':CI_low}
    if ID not in ID_list: ID_list.append(ID)
    strains.append(strain)
  infile.close()
  return epi_dict, list(set(strains)), ID_list

def formating_epi_file(epi_dict, strains, ID_list, out_file1, out_file2):
  print "writing: %s" % out_file1
  print "writing: %s" % out_file2
  out_file1 = open(out_file1, 'w')
  out_file2 = open(out_file2, 'w')
  out_file1.write("\t".join(['ID','mut1','mut2'])+"\t"+
                 "\t".join(strains)+"\t"+
                 "\t".join(['epi_mean','epi_std','CI_pos','CI_neg','CI_zero','class'])+"\n")
  out_file2.write("\t".join(['ID','mut1','mut2','strain','class'])+"\n")
  for ID in ID_list:
    ID_reform = ID.rsplit('-')
    mut1 = ID_reform[0][-1]+ID_reform[0][0:-1]
    mut2 = ID_reform[1][-1]+ID_reform[1][0:-1]
    ID_reform = mut1+'-'+mut2
    out  = [ID_reform,mut1,mut2]
    epis = []
    count_CI_pos  = 0
    count_CI_neg  = 0
    count_CI_zero = 0
    for strain in strains:
      e = epi_dict[strain][ID]['e']
      CI_up  = epi_dict[strain][ID]['CI_up']
      CI_low = epi_dict[strain][ID]['CI_low']
      epi_class=0
      if CI_low > 0:  count_CI_pos+=1; epi_class=1
      elif CI_up < 0: count_CI_neg+=1; epi_class=-1
      elif CI_low < 0 and CI_up > 0: count_CI_zero+=1; epi_class=0
      else: print 'something is wrong with the CI classification settings'
      out.append(e)
      epis.append(e)
      out_file2.write("\t".join(map(str,[ID_reform,mut1,mut2,strain,epi_class]))+"\n")
    out+=[np.mean(epis), np.std(epis), count_CI_pos, count_CI_neg, count_CI_zero]
    if min(epis) > 0: out.append('pos')
    elif max(epis) < 0: out.append('neg')
    else: out.append('mix')
    out_file1.write("\t".join(map(str,out))+"\n")
  out_file1.close()
  out_file2.close()
      
def main():
  for model in ['nonlinearity','linearity']:
    epi_file  = 'result/pairwiseparameters_'+model+'.csv'
    out_file1  = 'result/Inf_heatmap_overall_'+model+'.tsv'
    out_file2  = 'result/Inf_heatmap_specific_'+model+'.tsv'
    epi_dict, strains, ID_list = reading_epi_file(epi_file)
    strains   = ['HK68', 'Bk79', 'Bei89', 'Mos99', 'Bris07', 'NDako16']
    formating_epi_file(epi_dict, strains, ID_list, out_file1, out_file2)

if __name__ == '__main__':
  main()
