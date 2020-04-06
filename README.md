# Complete and validated genomes from a metagenome

Daniel J Giguere, Alexander T Bacheli

This repository is intended to accompany the pre-print (TODO insert link) by providing the code used, as well as a few explanations. 

### Counting long read coverage at beginning and end of fasta reference

When filtering Nanopore reads by > 90% query coverage, apparent read coverage at the beginning and end of the reference genome will appear to decrease because `minimap2` does not consider circular genomes when mapping. Therefore we applied the following strategy:

### Orienting genomes on dnaA gene

We applied the following strategy to orient the genomes: 

### Generating genome plots using the circlize R package 

The R package `circlize` was used to generate each of the genome figure. A general script is available for each of the data. 
