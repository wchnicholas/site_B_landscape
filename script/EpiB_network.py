#!/usr/bin/python
import os
import sys
import operator
import colorsys
import networkx as nx
import numpy as np
from math import exp
from itertools import imap
from collections import defaultdict, Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def codondistancemap():
  codondists = {}
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  aas = list(set(dnamap.values()))
  for aa1 in aas:
    for aa2 in aas:
      codons1 = [codon for codon in dnamap.keys() if dnamap[codon] == aa1]
      codons2 = [codon for codon in dnamap.keys() if dnamap[codon] == aa2]
      dist    = min([hamming(codon1, codon2) for codon1 in codons1 for codon2 in codons2])
      codondists[aa1+'-'+aa2] = dist
      codondists[aa2+'-'+aa1] = dist
  return codondists

def MinNucDist(pep1, pep2, codondists):
  return sum([codondists[aa1+'-'+aa2] for aa1, aa2 in zip(pep1, pep2)])

def reading_pref_file(filename):
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

def SelectMuts(pref_dict, strains):
  muts = []
  for mut in pref_dict['HK68'].keys():
    prefs = [float(pref_dict[strain][mut]) for strain in strains]
    if max(map(float,prefs)) > -99:
      muts.append(mut)
  return muts

def buildgraph(nodes):
  G = nx.Graph()
  [G.add_node(node) for node in nodes]
  for n1 in G.nodes():
    for n2 in G.nodes():
      if hamming(n1,n2) == 1: G.add_edge(n1,n2)
  return G

def labelnode(fit):
  high = float(2.5)
  low  = float(0)
  mid  = float(high-low)/2+low
  if fit > high:  return colorsys.rgb_to_hsv(1, 0, 0)
  elif fit > mid: return colorsys.rgb_to_hsv(1,(high-fit)/(high-mid),0)
  elif fit > low: return colorsys.rgb_to_hsv(1,1,(mid-fit)/(mid-low))
  elif fit <= low: return colorsys.rgb_to_hsv(1,1,1)
  else: print 'Something is wrong with the coloring function'; print fit; sys.exit()

def drawgraph(G,outfilename,pref_dict,strain,codondists):
  outfile = open(outfilename,'w')
  outfile.write('strict graph{'+"\n"+"\t"+'rankdir=LR'+
                                "\n"+"\t"+'node [shape=circle; label=""]'+
                                "\n"+"\t"+'overlap = false'+
                                "\n"+"\t"+'splines = true'+"\n")
  for var in G.nodes():
    fit = float(pref_dict[strain][var])
    col = labelnode(fit)
    col = ','.join(map(str,list(col)))
    outfile.write("\t"+var+' [fillcolor="'+col+'", color=black, penwidth=2, style="filled,rounded"];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    fit1 = float(pref_dict[strain][var1])
    fit2 = float(pref_dict[strain][var2])
    if hamming(var1,var2) == 1 and MinNucDist(var1, var2, codondists) == 1:
      outfile.write("\t"+var1+' -- '+var2+' [style=solid, color=grey, width=0.1];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def main():
  preffile = 'result/data_pref.tsv'
  outfile_prefix = 'dot/EvoNetwork_'
  strains   = ['HK68','Bk79','Bei89','Mos99','Bris07','NDako16']
  codondists = codondistancemap()
  pref_dict  = reading_pref_file(preffile)
  muts       = SelectMuts(pref_dict, strains)
  G = buildgraph(muts)
  for strain in strains:
    outfilename = outfile_prefix+strain+'.dot' 
    drawgraph(G,outfilename,pref_dict,strain,codondists)
    os.system('fdp -Tpng %s -o %s' % (outfilename,outfilename.replace('.dot','.png')))
    print 'Finish drawing for: %s' % strain

if __name__ == '__main__':
  main()
