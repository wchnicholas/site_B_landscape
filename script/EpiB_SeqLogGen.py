#!/usr/bin/python
import os
import sys
import glob
from math import exp
from collections import defaultdict

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
    pref_dict[strain][ID] = fit
  infile.close()
  return pref_dict

def main():
  filename  = "result/data_pref.tsv"
  outfile   = 'result/seqlogo.fa'
  pref_dict = reading_pref_file(filename)
  strains   = ['HK68','Bk79','Bei89','Mos99','Bris07','NDako16']
  multiply  = 100
  for strain in strains:
    outfile_tmp   = outfile.replace('.fa','_'+strain+'.fa')
    print "writing: %s" % outfile_tmp
    outfile_tmp   = open(outfile_tmp,'w')
    for mut in pref_dict[strain].keys():
      fit  = pref_dict[strain][mut]
      N    = int(exp(fit)*multiply)
      for n in range(0,N):
	outfile_tmp.write('>'+strain+"\n"+mut+"\n")
    outfile_tmp.close()
    os.system('/Users/wchnicholas/Bioinformatics/weblogo-3.6.0/weblogo -A protein -f result/seqlogo_'+strain+'.fa  -U probability --resolution 600 -F png -s large --aspect-ratio 2 -c chemistry -X NO -Y NO -o graph/seqlogo_'+strain+'.png')

if __name__ == "__main__":
  main()
