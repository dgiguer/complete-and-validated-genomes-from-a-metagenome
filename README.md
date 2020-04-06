# Complete and validated genomes from a metagenome

Daniel J Giguere, Alexander T Bacheli

This repository is intended to accompany the pre-print (TODO insert link) by providing the code used, as well as a few explanations. 

### Orienting genomes on dnaA gene

We applied the following strategy to orient the genomes:

1. Annotate the genome using prokka.

2. Identify the start of the dnaA genes using prokka (using a grep "dnaA" command). If multiple dnaA genes exist, extract the first dnaA identified by prokka. 

3. Using a python script (orient.py), re-orient the genome so that the start of the dnaA gene correspond with the start of the genome. 

### Counting long read coverage at beginning and end of fasta reference

After filtering Nanopore reads by > 90% query coverage, apparent read coverage at the beginning and end of the reference genome will appear to decrease because `minimap2` does not consider circular genomes when mapping. Therefore we applied the following strategy:

1. Using the oriented genome, identify the (approximately) halfway point of the genome based on its length.

2. Using python (modify-genome.py), make a re-oriented copy of the genome called "reverse-start" that begins at the halfway point of the oriented genome.

3. Map long reads to both the forward and reverse-oriented genomes, separately. 

4. To determine the long read coverage for the first and last quarters of the complete genome calculate the coverage of every 1000bp of the second and third quarters of the reverse-oriented genome. Because the second and third quarters of the reverse-oriented genome is relatively far from the origin of the reverse-oriented genome, we can assume the filtered reads were accruately mapped to this region by `minimap2`.

5. To determine the long read coverage for the second and third quarters of the complete genome calculate the coverage of every 1000bp of the second and third quarters of the forward-oriented genome. 

6. Add the coverage from each quarter of the genome together to get the complete genome coverage.

7. The final 1000bp segment of the genome is assumed to be 1000bp in length, and it likely includes coverage of part of the start of the genome. To accruately determine the coverage of this end segment of the genome, multiply the calculated coverage by the fraction of the size of the end segment of the genome divided by 1000bp. 

### Generating genome plots using the circlize R package 

The R package `circlize` was used to generate each of the genome figure. A general script is available for each of the data. 
