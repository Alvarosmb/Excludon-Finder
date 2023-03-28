# Operon-Finder
Easy to use pipeline for overlapping UTRs and Non-contiguous operons annnotation using Oxford Nanopore direct RNA sequencing.

This project's code is adapted from Lindsay Clark's (https://github.com/lvclark) repository, (https://github.com/HPCBio/lasa-ONT-2022June), with modifications made to create a more efficient and versatile pipeline that only requires input and utilizes a conda environment for ease of use.

## Getting started
### Installing conda enviroment
```
git clone https://github.com/Alvarosmb/Operon-Finder.git
cd Operon-Finder 
conda env create -f environment.yml  
conda activate env_name  
```
## Alignment and depth covergare annotation 
The script ```Scripts/ExcludonFinder ``` is used to:  
Align reads to references using minimap2  
Sort and index BAMs and determine the total coverage at each genomic position using Samtools

## Overview of the steps performed for Non-contiguous operon identification:
Identifying overlapping UTRs and non-contiguous operon with the script ```Scripts/Overlapping_reads.R ```
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
  ## Excludon annotation:
  Count of overlapping 3'UTRs and 5' UTRs in the genome regardless of whether they form an operon or not. Median read ends must overlapp at least 20 nt of the neighbor gene.


## Usage
```
Excludon_finder -f <fasta_input> -q <fastq_input> -g <gff_input> [-m mode] [-r split ratio]

```
 * ```-f --fasta FASTA FILE ```
    * Reference genome 
 * ```-q --fastq FASTQ FILE ```
    * Sequence data
 * ```-g --gff GFF FILE ```
    * Annotation file 
 * ```-m --mode MODE OPTION ```
    * Specify if it is only desired to identify either Operons or Excludons (OPTIONAL). Options "Operon" or "Excludon". Both will be identifyied by default.
  * ```-r --split_ratio SPLIT RATIO ```
      * Asses the value of split ratio (OPTIONAL). Default value: 10

## Examples
In the folder ExampleFiles, the following Use the following to obtain output files from the example:input files are provided as an example:  

### Expected output files:  
***{Sample}NonContiguousOperons.csv***      -Non-contiguous operon predictions file  
***{Sample}Overlapping_3UTRs.csv***         -3'UTRs predictions file  
***{Sample}Overlapping_5UTRs.csv***         -5' UTRs predictions file  
***{Sample}Operons.bed***                   -Genomic regions file for operons visualization


