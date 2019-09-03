#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from Bio import SeqIO
from math import log10, exp
from collections import defaultdict, Counter

def reading_data(filename):
  fit_dict = defaultdict(dict)
  infile = open(filename, 'r')
  for line in infile.xreadlines():
    if 'InputCount' in line: continue
    line = line.rstrip().rsplit(',')
    ID     = line[0]
    strain = line[10]
    strain = 'Bris07'if strain=='Bris07L194' else strain
    fit    = line[13]
    fit_dict[strain][ID] = fit
  infile.close()
  return fit_dict

def fit_to_pref(fit_dict):
  pref_dict = defaultdict(dict)
  for strain in fit_dict.keys():
    totalID  = len(fit_dict[strain].keys())
    avgfit   = np.mean(map(lambda x:float(x),fit_dict[strain].values()))
    sdfit    = np.std(map(lambda x:float(x),fit_dict[strain].values()))
    print "strain:", strain, "; # of variant:", totalID, "; avg fit:", avgfit, "; std fit", sdfit
    for ID in fit_dict[strain].keys():
      pref_dict[strain][ID] = (float(fit_dict[strain][ID])-avgfit)/sdfit
  return pref_dict

def write_file(fit_dict,pref_dict,outfile):
  print "writing: %s" % outfile 
  outfile = open(outfile,'w')
  outfile.write("\t".join(['ID','strain','fit','pref'])+"\n")
  for strain in fit_dict.keys():
    for ID in fit_dict[strain].keys():
      fit  = fit_dict[strain][ID]
      pref = pref_dict[strain][ID]
      outfile.write("\t".join(map(str,[ID,strain,fit,pref]))+"\n")
  outfile.close()

def main():
  filename  = "result/data_all.csv"
  outfile   = "result/data_pref.tsv"
  fit_dict  = reading_data(filename)
  pref_dict = fit_to_pref(fit_dict)
  write_file(fit_dict,pref_dict,outfile)

if __name__ == "__main__":
  main()
