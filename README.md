## ANALYSIS FOR H3N2 ANTIGENIC SITE B DEEP MUTATIONAL SCANNING
This study aims to understand how the local fitness landscape of antigenic site B in H3N2 HA evolves in the past 50 years. The repository here describes the analysis for the deep mutational scanning experiment that focuses on HA1 residues 156, 158, 159, 190, 193, 196.

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA563320](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA563320), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [./data/SampleID.tsv](./data/SampleID.tsv): Describes the sample identity for each fastq file
* [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa): Reference (wild type) sequences

### ANALYSIS PIPELINE FOR HK68 (H3N2), PERTH09 (H3N2), and WSN (H1N1)
1. [./script/EpiB\_fastq\_to\_fitness.py](./script/EpiB_fastq_to_fitness.py): Converts raw reads to variant counts and fitness measures.
    - Input files:
      - Raw sequencing reads in fastq/ folder
      - [./data/SampleInfo.tsv](./data/SampleInfo.tsv)
      - [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa)
    - Output files: result/EpiB\_MultiMutLib\_\*.tsv
