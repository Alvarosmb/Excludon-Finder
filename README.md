# Operon-Finder
Easy to use pipeline for overlapping UTRs and Non-contiguous operons annnotation using Oxford Nanopore direct RNA sequencing.

The code used in this project is adapted from (https://github.com/HPCBio/lasa-ONT-2022June) created by Lindsay Clark (https://github.com/lvclark). The original code has been modified to create a less time-consuming and more general pipeline that requires only input and uses a conda environment to simplify usage.

## Getting started
### Installing conda enviroment
```
git clone https://github.com/username/repo.git  
cd repo  
conda env create -f environment.yml  
conda activate env_name  
```
## Alignment and depth covergare annotation 
The script ... is used to:  
Align reads to references using minimap2  
Sort and index BAMs and determine the total coverage at each genomic position using Samtools

## Overview of the steps performed for operon identification:
Identifying overlapping UTRs and non-contiguous operon with the script ...    
  1. Determine start and stop points for the alignment for each read
  2. Count of mapped reads
  3. Filter reads with more than 80% alignment identity to the reference  and overlapping CDS, tRNA, or tmRNA genes
  4. Overlapping reads (on the same strand) were combined into ranges representing putative operons
  5. Operon identification: 
    - Iteratioin through each gene of the putative operon taking the median coverage of 20 nt upstream of the TTS of the gene. 
    - Divide this median coverage by a split ratio of 20. Check if, for each nt of the intergenic region to the next gene, the coverage is less than this ratio. 
    - If 3 consecutive nt have less coverage thant thse ratio, the operon is splitted.
  6. Operons were labeled as non-contiguous if an operon on the opposite strand was flanked on both sides by genes in the operon in question.
  7. Count of reads that overlapp genes flinking the reverse gene(s) of the operon.
  ## UTRs annotation:
  Count of 3'UTRs and 5' UTRs in the genome regardless of whether they form an operon or not




