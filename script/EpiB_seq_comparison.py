#!/usr/bin/python
import os
import sys
import glob
import string
import operator
import numpy as np
from Bio import SeqIO
from scipy import stats
from itertools import imap
from collections import defaultdict, Counter

def hamming(str1, str2):
  assert len(str1) == len(str2)
  return sum(imap(operator.ne, str1, str2))

def extract_RBS_seq(seq_dict, resi_offset, RBS_residues):
  RBS_dict = {}
  for ID in seq_dict.keys():
    seq = seq_dict[ID]
    RBS = ''.join([seq[RBS_residue-resi_offset] for RBS_residue in RBS_residues])
    RBS_dict[ID] = RBS
    print ID+"\t"+RBS.replace('',"\t")
  return RBS_dict

def reading_fasta(infile):
  records  = SeqIO.parse(infile,"fasta")
  seq_dict = {}
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    seq_dict[ID] = seq
  return seq_dict

def reading_pref(filename):
  pref_dict = defaultdict(dict)
  infile = open(filename,'r')
  for line in infile.xreadlines():
    if 'pref' in line: continue
    line = line.strip().rsplit("\t")
    ID     = line[0]
    strain = line[1]
    fit    = float(line[2])
    pref   = float(line[3])
    pref_dict[strain][ID] = pref
  infile.close()
  return pref_dict

def compute_pref_cor(pref_dict):
  cor_dict = {}
  strains = pref_dict.keys()
  muts    = pref_dict['HK68'].keys()
  for strain1 in strains:
    for strain2 in strains:
      if strain1 != strain2:
        strain1_prefs = [pref_dict[strain1][mut] for mut in muts]
        strain2_prefs = [pref_dict[strain2][mut] for mut in muts]
        cor = stats.pearsonr(strain1_prefs, strain2_prefs)[0]
        cor_dict[strain1+'-'+strain2] = cor
  return cor_dict

def strain_to_ID(strain):
  if strain=='A/HongKong/1/1968': return 'HK68'
  elif strain=='A/Bangkok/1/1979': return 'Bk79'
  elif strain=='A/Beijing/353/1989': return 'Bei89'
  elif strain=='A/Moscow/10/1999': return 'Mos99'
  elif strain=='A/Brisbane/10/2007': return 'Bris07'
  elif strain=='A/North_Dakota/26/2016': return 'NDako16'
  else: 
    print "strain cannot be converted to ID"
    sys.exit()

def out_pairwise_dist(seq_dict, RBS_dict, cor_dict, outfile):
  print "writing: %s" % outfile
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['strain1', 'strain2', 'HA_identity', 'RBS_identity', 'landscape_cor'])+"\n")
  strains = [strain for strain in seq_dict.keys()]
  for n in range(len(strains)):
    for m in range(len(strains)):
      if n < m:
        strain1 = strains[n]
        strain2 = strains[m]
        ID1 = strain_to_ID(strain1)
        ID2 = strain_to_ID(strain2)
        HA_identity  = len(seq_dict[strain1])-hamming(seq_dict[strain1], seq_dict[strain2])
        RBS_identity = len(RBS_dict[strain1])-hamming(RBS_dict[strain1], RBS_dict[strain2])
        cor      = cor_dict[ID1+'-'+ID2]
        outfile.write("\t".join(map(str,[ID1, ID2, HA_identity, RBS_identity, cor]))+"\n")
  outfile.close()

def main():
  fasta_infile = "Fasta/HA_ecto.fa"
  outfile      = "result/EpiB_seq_dist.tsv"
  preffile     = "result/data_pref.tsv"
  resi_offset  = 11
  RBS_residues = [135, 137, 145, 222, 225, 226]
  pref_dict = reading_pref(preffile)
  seq_dict  = reading_fasta(fasta_infile)
  RBS_dict  = extract_RBS_seq(seq_dict, resi_offset, RBS_residues)
  cor_dict  = compute_pref_cor(pref_dict)
  out_pairwise_dist(seq_dict, RBS_dict, cor_dict, outfile)

if __name__ == "__main__":
  main()
