## ANALYSIS FOR H3N2 ANTIGENIC SITE B DEEP MUTATIONAL SCANNING
This study aims to understand how the local fitness landscape of antigenic site B in H3N2 HA evolves in the past 50 years. The repository here describes the analysis for the deep mutational scanning experiment that focuses on HA1 residues 156, 158, 159, 190, 193, 196.

### REQUIREMENTS
* [Python](https://www.python.org/) version 2.7
* [Weblogo](https://weblogo.berkeley.edu) version 3.6

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA563320](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA563320), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2"
* [./data/SampleID.tsv](./data/SampleID.tsv): Describes the sample identity for each fastq file
* [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa): Reference (wild type) nucleotide sequences for the sequencing data
* [./data/WTseq.tsv](./data/WTseq.tsv): Amino acids for the wild type sequences at residues 156, 158, 159, 190, 193, 196
* [./Fasta/HumanH3N2\_All\_2018.aln](./Fasta/HumanH3N2\_All\_2018.aln): Full-length HA protein sequences from human H3N2 downloaded from [GISAID](https://www.gisaid.org/)

### ANALYSIS PIPELINE
1. [./script/EpiB\_fastq\_to\_fitness.py](./script/EpiB_fastq_to_fitness.py): Converts raw reads to variant counts and fitness measures.
    - Input files:
      - Raw sequencing reads in fastq/ folder
      - [./data/SampleInfo.tsv](./data/SampleInfo.tsv)
      - [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa)
    - Output files:
      - result/EpiB\_MultiMutLib\_\*.tsv
2. [./script/EpiB\_clean\_mut.py](./script/EpiB\_clean\_mut.py): Filter mutants of interest
    - Input files:
      - result/EpiB\_MultiMutLib\_\*.tsv
    - Output files:
      - result/EpiB\_Index\_\*.tsv
3. [./script/combine\_data.jl](./script/combine\_data.jl): Re-calculate mutant fitness. Written by [Jakub Otwinowski](https://github.com/jotwin)
    - Input files:
      - result/EpiB\_Index\_\*.tsv
    - Output files: 
      - [./result/data\_all.csv](./result/data\_all.csv):
4. [./script/EpiB\_fit\_to\_pref.py](./script/EpiB\_fit\_to\_pref.py): Preference normalization
    - Input files:
      - [./result/data\_all.csv](./result/data\_all.csv):
    - Output files:
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
5. [./script/EpiB\_PrefEvol.py](./script/EpiB\_PrefEvol.py): Amino acid sequences of HA residues 156, 158, 159, 190, 193, and 196 in naturally occurring strains were extracted
    - Input files:
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
      - [./Fasta/HumanH3N2\_All\_2018.fa](./Fasta/HumanH3N2\_All\_2018.fa)
    - Output files:
      - [./result/HumanH3N2_HAecto_year.tsv](./result/HumanH3N2_HAecto_year.tsv)
      - [./result/Prefs\_ByYear.tsv](./result/Prefs\_ByYear.tsv)
      - [./result/Motif\_ByYear.tsv](./result/Motif\_ByYear.tsv)
6. [./script/EpiB\_seq\_comparison.py](./script/EpiB\_seq\_comparison.py):
    - Input files:
      - [./Fasta/HA_ecto.fa](./Fasta/HA_ecto.fa)
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
    - Output files:
      - [./result/EpiB_seq_dist.tsv](./result/EpiB_seq_dist.tsv)

### PLOTTING
1. [./script/Plot\_CompareRep.R](./script/Plot\_CompareRep.R): Compare mutant fitness (i.e. enrichment ratio) from replicates
    - Input files:
      - [./result/data\_all.csv](./result/data\_all.csv)
    - Output files:
      - [./graph/Compare\_Rep.png](./graph/Compare\_Rep.png)
2. [./script/Plot\_CompareLib.R](./script/Plot\_CompareLib.R): Compare mutant fitness from different genetic backgrounds
    - Input files:
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
    - Output files:
      - [./graph/Lib\_pref\_sina.png](./graph/Lib\_pref\_sina.png)
      - [./graph/LibCorPairs.png](./graph/LibCorPairs.png)
3. [./script/EpiB\_SeqLogGen.py](./script/EpiB\_SeqLogGen.py): Generate sequence logo based on mutant preference
    - Input files:
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
    - Output files:
      - result/seqlogo\_\*.fa
      - graph/seqlogo\_\*.png
4. [./script/EpiB\_network.py](./script/EpiB\_network.py): Plot fitness landscape (network graph)
    - Input files:
      - [./result/data\_pref.tsv](./result/data\_pref.tsv)
    - Output files:
      - dot/EvoNetwork\_*.dot
      - dot/EvoNetwork\_*.png
5. [./script/Plot\_TrackPref.R](./script/Plot\_TrackPref.R): Plot the normalized preference of naturally occurring sequences in different genetic backgrounds
    - Input files:
      - [./result/Prefs\_ByYear.tsv](./result/Prefs\_ByYear.tsv)
    - Output files:
      - [./graph/PrefsByYear.png](./graph/PrefsByYear.png)
6. [./script/Plot_seq_dist.R](./script/Plot_seq_dist.R): Plot the relationship between pairwise correlation of fitness landscape and pairwise sequence identity
    - Input files:
      - [./result/EpiB_seq_dist.tsv](./result/EpiB_seq_dist.tsv)
    - Output files:
      - [./graph/Identity_vs_cor_HA.png](./graph/Identity_vs_cor_HA.png)
      - [./graph/Identity_vs_cor_RBS.png](./graph/Identity_vs_cor_RBS.png)
7. [./script/Plot_TrackFreq.R](./script/Plot_TrackFreq.R): Plot the frequency of different haplotypes over time
    - Input files:
      - [./result/Motif_ByYear.tsv](./result/Motif_ByYear.tsv)
    - Output files:
      - [./graph/FreqByYear.png](./graph/FreqByYear.png)
