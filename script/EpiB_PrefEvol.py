#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from Bio import SeqIO
from scipy import stats
from math import log10, log, exp
from collections import defaultdict, Counter

def reading_pref_file(filename):
  pref_dict = defaultdict(dict)
  motif_576 = []
  infile = open(filename,'r')
  for line in infile.xreadlines():
    if 'pref' in line: continue
    line = line.strip().rsplit("\t")
    ID     = line[0]
    strain = line[1]
    fit    = float(line[2])
    pref   = float(line[3])
    pref_dict[strain][ID] = pref
    motif_576.append(ID)
  infile.close()
  return pref_dict, list(set(motif_576))

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def seqs2consensus(seqlist):
  consensus = ''
  for n in range(len(seqlist[0])):
    resi = []
    for seq in seqlist:
      resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

def motifextracting(alnfile, residues):
  records   = [record for record in SeqIO.parse(alnfile,"fasta")]
  motifdict = defaultdict(list)
  for record in records:
    header = str(record.id)
    if header.count('|')!=4: continue
    ID    = header.rsplit('|')[0]
    PSG   = header.rsplit('|')[1]
    year  = header.rsplit('|')[-1][1:5]
    seq   = str(record.seq)
    assert(isInt(year))
    motif = ''.join([seq[pos] for pos in residues])
    motifdict[year].append(motif)
  return motifdict

def writing_motiffile(motifdict, motiffile, motif_576):
  natural_motifs = list(set([motif for year in motifdict.keys() for motif in motifdict[year]]))
  natural_motifs = set(natural_motifs).intersection(set(motif_576))
  print "A total of %i out of 576 motifs observed in naturally circulating strains " % len(natural_motifs)
  print "writing: %s" % motiffile
  outfile   = open(motiffile, 'w')
  outfile.write("\t".join(['year', 'motif', 'freq'])+"\n")
  for year in sorted(motifdict.keys(), key=lambda x:int(x)):
    motifs = Counter(motifdict[year])
    total_count = sum([int(motifs[motif]) for motif in motifs.keys()])
    for motif in natural_motifs:
      count = motifs[motif] if motif in motifs.keys() else 0
      freq = float(count)/float(total_count)
      outfile.write("\t".join(map(str, [year, motif, freq]))+"\n")
  outfile.close()

def ConsensusbyYear(motifdict, yearfile):
  yeardict = {}
  print "writing: %s" % yearfile
  outfile = open(yearfile, 'w')
  outfile.write('year'+"\t"+'num_seq'+"\n")
  for year in sorted(map(int,motifdict.keys())):
    seqs = motifdict[str(year)]
    yeardict[year] = seqs2consensus(seqs)
    num_seq = len(seqs)
    outfile.write(str(year)+"\t"+str(num_seq)+"\n")
  outfile.close()
  return yeardict

def analyze_pref_evol(pref_dict, outfile, motifdict, yeardict):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['Background','Year','MeanPref','StdPref'])+"\n")
  for strain in pref_dict.keys():
    print 'Working on %s' % strain
    for year in sorted(map(int,motifdict.keys())):
      motifs = motifdict[str(year)]
      prefs  = []
      if year==1968: print Counter(motifs)
      for motif in motifs:
	if motif in pref_dict[strain].keys():
	  prefs.append(pref_dict[strain][motif])
      outfile.write("\t".join(map(str,[strain, year, np.mean(prefs), stats.sem(prefs)]))+"\n")
  outfile.close()

def main():
  filename  = "result/data_pref.tsv"
  alnfile   = 'Fasta/HumanH3N2_All_2018.fa'
  yearfile  = 'result/HumanH3N2_HAecto_year.tsv'
  motiffile = 'result/Motif_ByYear.tsv'
  outfile   = 'result/Prefs_ByYear.tsv'
  pos_offset = 15
  residues   = [156,158,159,190,193,196]
  residues   = [resi+pos_offset for resi in residues]
  pref_dict, motif_576 = reading_pref_file(filename)
  motifdict = motifextracting(alnfile, residues)
  yeardict  = ConsensusbyYear(motifdict, yearfile)
  writing_motiffile(motifdict, motiffile, motif_576)
  analyze_pref_evol(pref_dict, outfile, motifdict, yeardict)
    
if __name__ == "__main__":
  main()
