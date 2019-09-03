#!/usr/bin/python
import os
import sys
import glob
from collections import defaultdict

def Mut2ID(Mut, WTseq, residues):
  ID = ''
  for residue, aa in zip(residues, WTseq):
    if str(residue) in Mut: ID += Mut.rsplit(str(residue))[1][0]
    else: ID += aa
  return ID

def ProcessTable(filename, Sample, inputcutoff, WTseq, residues, outfilename):
  print "Reading file %s" % filename
  infile  = open(filename,'r')
  Count_InputDict = defaultdict(int)
  Count_Rep1Dict  = defaultdict(int)
  Count_Rep2Dict  = defaultdict(int)
  CountMut = 0
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    if 'Mut' in line: continue
    Mut   = line[0]
    InputCount = float(line[2])
    Rep1Count  = float(line[3])
    Rep2Count  = float(line[4])
    if InputCount >= inputcutoff:
      CountMut += 1
      ID = Mut2ID(Mut, WTseq, residues)
      Count_InputDict[ID] = InputCount
      Count_Rep1Dict[ID]  = Rep1Count
      Count_Rep2Dict[ID]  = Rep2Count
  infile.close()
  print "Total Mut in Sample %s: %i (Input cutoff = %i)" % (Sample, CountMut, inputcutoff)
  return Count_InputDict, Count_Rep1Dict, Count_Rep2Dict

def writeIDfreqfile(Count_InputDict, Count_Rep1Dict, Count_Rep2Dict, outfilename): 
  print "Writing file %s" % outfilename
  outfile = open(outfilename,'w')
  outfile.write("\t".join(['ID','InputCount','Rep1Count','Rep2Count',
                           'Rep1Fitness', 'Rep2Fitness', 'Fitness',
                           'Rep1_Preference', 'Rep2_Preference', 'Preference'])+"\n")
  IDs = Count_InputDict.keys()
  InputTotal = sum(Count_InputDict.values())
  Rep1Total = sum(Count_Rep1Dict.values())
  Rep2Total = sum(Count_Rep2Dict.values())
  TotalFit  = 0
  Rep1_TotalFit = 0
  Rep2_TotalFit = 0
  TotalID   = len(IDs)
  outdata   = []
  for ID in IDs:
    InputCount = Count_InputDict[ID]
    Rep1Count  = Count_Rep1Dict[ID]
    Rep2Count  = Count_Rep2Dict[ID]
    InputFreq  = InputCount/InputTotal
    Rep1Freq   = Rep1Count/Rep1Total
    Rep2Freq   = Rep2Count/Rep2Total
    Rep1Fitness = Rep1Freq/InputFreq
    Rep2Fitness = Rep2Freq/InputFreq
    Fitness     = (Rep1Fitness+Rep2Fitness)/2
    Rep1_TotalFit+=Rep1Fitness
    Rep2_TotalFit+=Rep2Fitness
    TotalFit += Fitness
    outdata.append(map(str,[ID, InputCount, Rep1Count, Rep2Count, Rep1Fitness, Rep2Fitness, Fitness]))
  for out in outdata:
    Rep1_Fitness = float(out[-3])
    Rep2_Fitness = float(out[-2])
    Fitness    = float(out[-1])
    Preference = Fitness/TotalFit*TotalID
    Rep1_Preference = Rep1_Fitness/Rep1_TotalFit*TotalID
    Rep2_Preference = Rep2_Fitness/Rep2_TotalFit*TotalID
    outfile.write("\t".join(out)+"\t"+"\t".join(map(str,[Rep1_Preference, Rep2_Preference, Preference]))+"\n")
  outfile.close()

def ReadingInfo(info_file):
  info_dict = {}
  infile = open(info_file, 'r')
  for line in infile.readlines():
    if 'ID' in line: continue
    line = line.rstrip().rsplit("\t")
    ID     = line[0]
    Sample = line[1]
    info_dict[ID] = Sample
  infile.close()
  return info_dict

def main():
  inputcutoff = 50
  filenames = glob.glob('result/*MultiMutLib*.tsv')
  WTseqfile = 'data/WTseq.tsv'
  residues  = [156,158,159,190,193,196]
  WTseqs    = ReadingInfo(WTseqfile)
  for filename in filenames:
    outfilename = filename.replace('MultiMutLib','Index')
    Sample      = filename.rsplit('.')[0].rsplit('_')[2]
    WTseq       = WTseqs[Sample]
    Count_InputDict, Count_Rep1Dict, Count_Rep2Dict = ProcessTable(filename, Sample, inputcutoff, WTseq, residues, outfilename)
    writeIDfreqfile(Count_InputDict, Count_Rep1Dict, Count_Rep2Dict, outfilename)

if __name__ == "__main__":
  main()
