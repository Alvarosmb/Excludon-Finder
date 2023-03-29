# Excludon-Finder
Easy to use pipeline for annnotation of overlapping UTRs and Non-contiguous operons  using Oxford Nanopore direct RNA sequencing.

This project's code is adapted from Lindsay Clark's (https://github.com/lvclark) repository, (https://github.com/HPCBio/lasa-ONT-2022June), with modifications made to create a more efficient and versatile pipeline that only requires input and utilizes a conda environment for ease of use.

## Getting started
### Installing conda enviroment
```
git clone https://github.com/Alvarosmb/Excludon-Finder.git
cd Excludon-Finder 
conda env create -f excludon_finder_environment.yml 
conda activate excludon_finder
```
## Alignment and depth covergare annotation 
The script ```Scripts/ExcludonFinder ``` is used to:  
Align reads to references using minimap2  
Sort and index BAMs and determine the total coverage at each genomic position using Samtools

 ## Excludon annotation:
Identification of  overlapping 3' UTRs and 5' UTRs in the genome regardless of whether they form an operon or not with the script ```Scripts/Overlapping_reads.R ```.
 * Identify convergent and divergent pairs of genes in the genome
 * Calculate median read end and start for these genes
 * For divergent genes, median read start of the pair of genes must overlap
 * For convergent genes, median read end of the pair of genes must overlap  
 
Overlapping length must be at least of 20 nt
 

## Overview of the steps performed for Non-contiguous operon identification:
  * Determine start and stop points for the alignment for each read
  *  Count of mapped reads
  * Filter reads with more than 80% alignment identity to the reference  and overlapping CDS, tRNA, or tmRNA genes
  * Overlapping reads (on the same strand) were combined into ranges representing putative operons
  * Non-contiguous Operon identification:   
    - Iteratioin through each gene of the putative operon taking the median coverage of 20 nt upstream of the TTS of the gene  
    - Divide this median coverage by a split ratio (default value set at 10).  
    - Check if, for each nt of the intergenic region to the next gene, the coverage is less than this ratio.   
    - If 3 consecutive nt have less coverage thant this ratio, the operon is splitted.  
  * Operons were labeled as non-contiguous if al least one gene on the opposite strand was flanked on both sides by genes in the operon in question.
  * Count of reads that overlapp genes flinking the reverse gene(s) of the operon.
 

## Usage
```
Excludon_finder -f <path_to_fasta_input> -q <path_to_fastq_input> -g <path_to_gff_input> [-m mode] [-r split ratio]

```
 * ```-f --fasta FASTA FILE PATH```
    * Reference genome 
 * ```-q --fastq FASTQ FILE PATH```
    * Sequencing data
 * ```-g --gff GFF FILE PATH```
    * Annotation file 
 * ```-m --mode MODE OPTION ```
    * Specify if it is only desired to identify either Operons or Excludons (OPTIONAL). Options "Operon" or "Excludon". Both will be identifyied by default.
  * ```-r --split_ratio SPLIT RATIO ```
      * Asses the value of split ratio (OPTIONAL). Default value: 10

## Examples
In the folder ExampleFiles, input files from _Autographa Californica nucleopolyhedrovirus_ (ACMNPV), from Baculoviridae family, are provided. This organism was chosen as an exemplar due to its relatively small genome size and therefore, the low computational complexity required for annotating non-contiguous operons. SPLIT RATIO [-r option ] must be set at 20 to detect a non-contiguous operon
#### Expected output files:  
***{Sample}NonContiguousOperons.csv***      -Non-contiguous operon predictions file  
***{Sample}Overlapping_3UTRs.csv***         -3'UTRs predictions file  
***{Sample}Overlapping_5UTRs.csv***         -5' UTRs predictions file  



