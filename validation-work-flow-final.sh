#!/bin/bash
# Dan Giguere, Alec Bahcheli, Ben Joris
# Before beginning, make sure you have the following:
  # all the forward and reverse short reads are in accessible directories
  # all the long (nanopore) reads are in a directory
  # the path to the original metaFlye assembly with all files maintained there

# set a size cutoff for the mags to meet in order to be analyzed: 300kb eliminates most chloroplast and mitochondria genomes
SIZE=300000
# dot plot parameters: setting to 'Yes' will generate smaller dot plots for every 50kb of the genome in addition to the whole genome dot plot
SMALLPLOTS='Yes'
# add an argument for the number of threads
THREADS=25
# define the short forward and reverse genes from illumina and the long reads from nanopore
FWD='/Volumes/data/algae/data/reads/illumina/algae_R1_strict.fastq.gz'
REV='/Volumes/data/algae/data/reads/illumina/algae_R2_strict.fastq.gz'
LONGREADS='/Volumes/data/algae/data/reads/guppy_3_3_0_hac_output/algae_hac_long_reads.fastq'
# define the directory containing the metaFlye assembly
RAWFASTA='/Volumes/data/cmags/data/flye_assembly_meta_hac'
# directory to work for the relevant project
LOCAL=''
# directory to work for reassemblies, annotations and dot plots
LOCAL1=`echo $LOCAL/polishing`
LOCAL2=`echo $LOCAL/re-assemblies`
LOCAL3=`echo $LOCAL/contig-annotations`
LOCAL4=`echo $LOCAL/dot-plots`
TABLE1=`echo $LOCAL/output/table1.txt`
TABLE2=`echo $LOCAL/output/table2.txt`
TABLE3=`echo $LOCAL/output/table3.txt`
# path to necessary tools: bowtie2, bowtie2-build, minimap2, samtools, bcftools, rebaler, reapr, pilon, unicycler, flye, tRNAscan-SE-2.0
BOWTIE='/Volumes/data/bin/bowtie2.3.5/bowtie2'
BUILD='/Volumes/data/bin/bowtie2.3.5/bowtie2-build'
MAP='/Volumes/data/bin/minimap2/minimap2'
SAMTOOLS='/Volumes/data/bin/samtools-1.10/samtools-1.10/samtools'
BCFTOOLS='/Volumes/data/bin/bcftools-1.10.2/bcftools'
REBALER='/usr/bin/python3.6 /Volumes/data/bin/Rebaler/rebaler-runner.py'
REAPR='/Volumes/data/bin/Reapr_1.0.18/reapr'
PILON='/Volumes/data/bin/pilon-1.23.jar'
UNICYCLER='/Volumes/data/bin/Unicycler/unicycler-runner.py'
FLYE='/Volumes/data/bin/Flye/bin/flye'
SPADES='/Volumes/data/bin/SPAdes-3.11.1-Linux/bin/spades.py'
TRNASCAN='/Volumes/data/bin/tRNAscan-SE-2.0/tRNAscan-SE'
# path to python
PYTHON='/usr/bin/python3.6'




# make a directory for the project (as defined in the variable LOCAL)
# output: $LOCAL
# output: $LOCAL/raw-contigs
# output: $LOCAL/pilon-logs
mkdir -p $LOCAL $LOCAL/raw-contigs $LOCAL1/pilon-logs $LOCAL/output $LOCAL/output/circos $LOCAL/output/dotplots





###################################

# snake make pipeline already determine for this step

###################################

# extract circular contigs that are over 300kb in size; input requires the metaFlye directory
# input: $RAWFASTA/assembly_info.txt
# input: $RAWFASTA/assembly.fasta 
# output: $LOCAL/raw-contigs/$name (series of different files, one for every circular contig larger than 200kb)
# requires: extract-circularized-fasta.py
./extract-circularized-fasta.py $RAWFASTA/assembly_info.txt $RAWFASTA/assembly.fasta $LOCAL/raw-contigs/mags-to-analyze.txt $LOCAL/raw-contigs $SIZE

# define the mags to analyze from the output of the last script
# input: $LOCAL/raw-contigs/mags-to-analyze.txt
# output: defined MAGS variable for future use
MAGS=`awk '{print $1}' $LOCAL/raw-contigs/mags-to-analyze.txt`




###################################

# validation and polishing processes

###################################

# repeat the following validation and polishing steps for each contig
for name in $MAGS
do

# create a directory for the contig validation, polishing, reassembly, etc
# output: $LOCAL/$name
# output: $LOCAL/$name/output
# output: $name/output/fcd_errors
# output: $name/output/total_coverage
# output: $LOCAL2
# output: $LOCAL2/$name
# output: $LOCAL2/$name/flye
# output: $LOCAL2/$name/unicycler
# output: $LOCAL2/$name/unicycler-long
# output: $LOCAL2/$name/spades
# output: $LOCAL3
# output: $LOCAL3/$name
# output: $LOCAL3/$name/output
# output: $LOCAL4
# output: $LOCAL4/$name
# output: $LOCAL/$name/circos
mkdir -p $LOCAL/$name $LOCAL/$name/output $LOCAL/$name/output/fcd_errors $LOCAL/$name/output/total_coverage $LOCAL2 $LOCAL2/$name $LOCAL2/$name/flye $LOCAL2/$name/unicycler $LOCAL2/$name/unicycler-long $LOCAL2/$name/spades $LOCAL3 $LOCAL3/$name $LOCAL3/$name/output $LOCAL4 $LOCAL4/$name $LOCAL/$name/circos

# link the raw fasta file to that directory for analysis
# input: $LOCAL/raw-contigs/$name
# output: $LOCAL/$name/$name-contigs.fa
cp $LOCAL/raw-contigs/$name $LOCAL/$name/$name-contigs.fa

# the header line needs to be changed to circular=true before rebaler
# input: $LOCAL/$name/$name-contigs.fa
# output: $LOCAL/$name/$name-contigs-edited.fa
sed '1 s/\(^>.*$\)/\1circular=true/' $LOCAL/$name/$name-contigs.fa > $LOCAL/$name/$name-contigs-edited.fa

# conncatenate the genomes into a single fasta file for minimap2 mapping
# input: $LOCAL/$name/$name-contigs-edited.fa
# output: $LOCAL1/raw-genomes.fa
cat $LOCAL/$name/$name-contigs-edited.fa >> $LOCAL1/raw-genomes.fa

done

# map the nanopore long reads to the mag
# input: $LOCAL1/raw-genomes.fa
# input: $LONGREADS
# output: $LOCAL1/raw-genome-nanopore.sam
$MAP -t $THREADS -aQLx map-ont --secondary=no --sam-hit-only $LOCAL1/raw-genomes.fa $LONGREADS > $LOCAL1/raw-genome-nanopore.sam

# extract relevant information on nanopore sequence alignments scores
# input: $LOCAL1/raw-genome-nanopore.sam
# output: $LOCAL1/sam.txt
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $LOCAL1/raw-genome-nanopore.sam > $LOCAL1/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
# input: $LOCAL1/cigars.txt
# output: $LOCAL1/cigar-results.txt
# requires: cigar-parse.py
./cigar-parse.py $LOCAL1/sam.txt $LOCAL1/final-read-names.txt

# filter the nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# this gets the sam headers in the file
# input: $LOCAL1/raw-genome-nanopore.sam
# output: $LOCAL1/nanopore-filtered-to-raw.sam
awk '$0 ~ "@"{print $0}' $LOCAL1/raw-genome-nanopore.sam > $LOCAL1/nanopore-filtered-to-raw.sam
# input: $LOCAL1/final-read-names.txt
# input: $LOCAL1/raw-genome-nanopore.sam
# output: $LOCAL1/nanopore-filtered-to-raw.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $LOCAL1/final-read-names.txt $LOCAL1/raw-genome-nanopore.sam >> $LOCAL1/nanopore-filtered-to-raw.sam

echo Finished filtering nanopore reads against raw assemblies

# activate conda virtual environment to run rebaler
source /Volumes/data/bin/miniconda3/bin/activate

for name in $MAGS
do

# get nanopore reads for a given mag
# input: $LOCAL1/nanopore-filtered-to-raw.sam
# output: $LOCAL/$name/nanopore-filtered-to-raw.sam
grep "^@" $LOCAL1/nanopore-filtered-to-raw.sam > $LOCAL/$name/nanopore-filtered-to-raw.sam
# input: $LOCAL1/nanopore-filtered-to-raw.sam
# output: $LOCAL/$name/nanopore-filtered-to-raw.sam
grep $name $LOCAL1/nanopore-filtered-to-raw.sam | grep -v "^@" >> $LOCAL/$name/nanopore-filtered-to-raw.sam

# intialize (sort) the sam file using samtools
# input: $LOCAL/$name/nanopore-filtered-to-raw.sam
# output: $LOCAL/$name/nanopore-filtered-to-raw-sorted.bam
$SAMTOOLS sort $LOCAL/$name/nanopore-filtered-to-raw.sam -o $LOCAL/$name/nanopore-filtered-to-raw-sorted.bam -@ $THREADS

# convert the sam file into fastq for use in polishing
# input: $LOCAL/$name/nanopore-filtered-to-raw-sorted.bam
# output: $LOCAL/$name/nanopore-filtered-to-raw.fastq
$SAMTOOLS fastq $LOCAL/$name/nanopore-filtered-to-raw-sorted.bam > $LOCAL/$name/nanopore-filtered-to-raw.fastq -@ $THREADS

# use rebaler to correct the assembly using nanopore
# changed to the edited file
# input:$LOCAL/$name/$name-contigs-edited.fa
# input: $LOCAL/$name/nanopore-filtered-to-raw.fastq
# output: $LOCAL/$name/rebaler-genome.fasta
$REBALER -t $THREADS $LOCAL/$name/$name-contigs-edited.fa $LOCAL/$name/nanopore-filtered-to-raw.fastq > $LOCAL/$name/rebaler-genome.fasta

# export the fastq and genomes
# input: $LOCAL/$name/rebaler-genome.fasta
# output: $LOCAL1/all-rebaler-genomes.fasta
cat $LOCAL/$name/rebaler-genome.fasta >> $LOCAL1/all-rebaler-genomes.fasta
# input: $LOCAL/$name/nanopore-filtered-to-raw.fastq
# output: $LOCAL/nanopore-filtered.fastq
cat $LOCAL/$name/nanopore-filtered-to-raw.fastq >> $LOCAL1/nanopore-filtered.fastq

done

echo Finished rebaler corrections 

# remap the nanopore against the rebaler assemblies
# input: $LOCAL1/all-rebaler-genomes.fasta
# input: $LOCAL1/nanopore-filtered.fastq
# output: $LOCAL1/nanopore-to-rebaler.sam
$MAP -t $THREADS -aLQx map-ont --sam-hit-only --secondary=no $LOCAL1/all-rebaler-genomes.fasta $LOCAL1/nanopore-filtered.fastq > $LOCAL1/nanopore-to-rebaler.sam

# make a bowtie2 database to map the illumina reads against the rebaler assemblies
# input: $LOCAL1/all-rebaler-genomes.fasta
# output: $LOCAL1/all-rebaler-polished has 6 files
# output files: all-rebaler-polished.1.bt2 all-rebaler-polished.2.bt2 all-rebaler-polished.3.bt2 all-rebaler-polished.4.bt2 all-rebaler-polished.rev.1.bt2 all-rebaler-polished.rev.2.bt2
$BUILD $LOCAL1/all-rebaler-genomes.fasta $LOCAL1/all-rebaler-polished

# map the illumina reads against the rebaler assemblies
# input: $LOCAL1/all-rebaler-polished database
# input: $FWD
# input: $REV
# output: $LOCAL1/all-illumina-to-rebaler.sam
$BOWTIE -x $LOCAL1/all-rebaler-polished -1 $FWD -2 $REV -p $THREADS -S $LOCAL1/all-illumina-to-rebaler.sam --no-unal --no-discordant --end-to-end --very-sensitive

# progress update
echo Mapped to rebaler

# sort the sam file
# input: $LOCAL1/nanopore-to-rebaler.sam
# output: $LOCAL1/nanopore-to-rebaler-sorted.bam
$SAMTOOLS sort -@ $THREADS $LOCAL1/nanopore-to-rebaler.sam -o $LOCAL1/nanopore-to-rebaler-sorted.bam

# index the bam file 
# input: $LOCAL1/nanopore-to-rebaler-sorted.bam
# output: $LOCAL1/nanopore-to-rebaler-sorted.bam.bai
$SAMTOOLS index -@ $THREADS $LOCAL1/nanopore-to-rebaler-sorted.bam

# remove non-primary alignments and convert to bam
# input: $LOCAL1/all-illumina-to-rebaler.sam
# output: $LOCAL1/all-illumina-to-rebaler.bam
$SAMTOOLS view -@ $THREADS -F 256 -bS $LOCAL1/all-illumina-to-rebaler.sam > $LOCAL1/all-illumina-to-rebaler.bam

# use samtools to sort the mapped reads
# input: $LOCAL1/all-illumina-to-rebaler.bam
# output: $LOCAL1/all-illumina-to-rebaler-sorted.bam
$SAMTOOLS sort $LOCAL1/all-illumina-to-rebaler.bam -o $LOCAL1/all-illumina-to-rebaler-sorted.bam -@ $THREADS

# index the sorted mapped reads creating and index
# input: $LOCAL1/all-illumina-to-rebaler-sorted.bam
# output: $LOCAL1/all-illumina-to-rebaler-sorted.bam.bai
$SAMTOOLS index -b $LOCAL1/all-illumina-to-rebaler-sorted.bam -@ $THREADS

# pilon does not requires subsetting the reads to those specific for a particular contig
# correct the rebaler assembly using both illumina and nanopore reads via pilon. Send logs to log directory. OUTPUTS MULTIPLE FILES
# input:$LOCAL1/all-rebaler-genomes.fasta
# input: $LOCAL1/all-illumina-to-rebaler-sorted.bam
# input: $LOCAL1/nanopore-to-rebaler-sorted.bam
# input: $LOCAL1/all-illumina-to-rebaler-sorted.bam.bai
# input: $LOCAL1/nanopore-to-rebaler-sorted.bam.bai
# output: $LOCAL1/pilon-logs/pilon.log
# output dir: $LOCAL1/pilon
# output includes 3 files: pilon.changes pilon.fasta pilon.fasta.fai
/usr/bin/java -Xmx32G -jar $PILON --genome $LOCAL1/all-rebaler-genomes.fasta --frags $LOCAL1/all-illumina-to-rebaler-sorted.bam --nanopore $LOCAL1/nanopore-to-rebaler-sorted.bam --fix snps,indels --threads $THREADS --verbose --outdir $LOCAL1/pilon --changes &> $LOCAL1/pilon-logs/pilon.log

sed 's/circular=true_pilon//g' $LOCAL1/pilon/pilon.fasta > $LOCAL1/all-pilon-genomes.fasta

# map the filtered nanopore reads to the pilon reassembly
# input: $LOCAL1/all-pilon-genomes.fasta
# input: $LOCAL1/nanopore-filtered.fastq
# output: $LOCAL1/nanopore-to-pilon.sam
$MAP -t $THREADS -aLQx map-ont --secondary=no $LOCAL1/all-pilon-genomes.fasta $LOCAL1/nanopore-filtered.fastq > $LOCAL1/nanopore-to-pilon.sam

# build a bowtie2 database to map the illumina reads against the pilon polished assemblies
# input: $LOCAL1/all-pilon-genomes.fasta
# output: $LOCAL1/all-pilon-genomes has 6 files
# output files: all-pilon-polished.1.bt2 all-pilon-polished.4.bt2 all-pilon-polished.2.bt2 all-pilon-polished.rev.1.bt2 all-pilon-polished.3.bt2 all-pilon-polished.rev.2.bt2
$BUILD $LOCAL1/all-pilon-genomes.fasta $LOCAL1/all-pilon-genomes

# re-map the short reads against the pilon polished references
# input: $LOCAL1/all-pilon-genomes database
# input: $FWD
# input: $REV
# output: $LOCAL1/illumina-to-all-pilon-genomes.sam
$BOWTIE -x $LOCAL1/all-pilon-genomes -1 $FWD -2 $REV -p $THREADS -S $LOCAL1/illumina-to-all-pilon-genomes.sam --no-unal --no-discordant --end-to-end --very-sensitive

# remove unmapped reads and convert to bam
# input: $LOCAL1/nanopore-to-pilon.sam
# output: $LOCAL1/nanopore-to-pilon.bam
$SAMTOOLS view -F 4 -bS $LOCAL1/nanopore-to-pilon.sam -o $LOCAL1/nanopore-to-pilon.bam -@ $THREADS

# use samtools to sort nanpore mapping to the pilon reassembly
# input: $LOCAL1/nanopore-to-pilon.bam
# output: $LOCAL1/nanopore-to-pilon-sorted.bam
$SAMTOOLS sort $LOCAL1/nanopore-to-pilon.bam -o $LOCAL1/nanopore-to-pilon-sorted.bam -@ $THREADS

# convert the illumina reads to bam and sort them
# input: $LOCAL1/illumina-to-all-pilon-genomes.sam
# output: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam
$SAMTOOLS sort $LOCAL1/illumina-to-all-pilon-genomes.sam > $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam -@ $THREADS

# progress update
echo Finished mapping to pilon

# index the reads in preparation for reapr
# input: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam
# output: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam.bai
$SAMTOOLS index $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam

# use reapr pipeline to check for errors in assembly according to illumina reads
# input: $LOCAL1/all-pilon-genomes.fasta
# input: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam
# input: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam.bai
# output dir: $LOCAL1/reapr has multiple files
$REAPR pipeline $LOCAL1/all-pilon-genomes.fasta $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam $LOCAL1/reapr

# progress update
echo Finished reapr

# generate mpileup files for variant calling
# input: $LOCAL1/all-pilon-genomes.fasta
# input: $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam
# output: $LOCAL1/mpileup.txt
$BCFTOOLS mpileup --threads $THREADS -o $LOCAL1/mpileup.txt -f $LOCAL1/all-pilon-genomes.fasta $LOCAL1/illumina-to-all-pilon-genomes-sorted.bam

# collect variants from mpileup
# input: $LOCAL1/mpileup.txt
# output: $LOCAL1/variants.txt
$BCFTOOLS call -cv -o $LOCAL1/variants.txt --threads $THREADS $LOCAL1/mpileup.txt

# filter for variants that have > 50% of reads supporting the call (indels)
# input: $LOCAL1/variants.txt
# output: $LOCAL1/variants-filtered.bcf
$BCFTOOLS filter -i '%MAX(IMF)>0.5' $LOCAL1/variants.txt -Ob > $LOCAL1/variants-filtered.bcf

# normalization is needed to left align each of the variants
# input: $LOCAL1/all-pilon-genomes.fasta
# input: $LOCAL1/variants-filtered.bcf
# output: $LOCAL1/variants-norm.bcf
$BCFTOOLS norm -d all -f $LOCAL1/all-pilon-genomes.fasta $LOCAL1/variants-filtered.bcf -Ob > $LOCAL1/variants-norm.bcf

# index these variants
# input: $LOCAL1/variants-norm.bcf
# output: $LOCAL1/variants-norm.bcf.csi
$BCFTOOLS index $LOCAL1/variants-norm.bcf -o $LOCAL1/variants-norm.bcf.csi

# apply all the variants to a final polished assembly
# input: $LOCAL1/all-pilon-genomes.fasta
# input: $LOCAL1/variants-norm.bcf
# output: $LOCAL1/final-assemblies.fasta
$BCFTOOLS consensus -f $LOCAL1/all-pilon-genomes.fasta $LOCAL1/variants-norm.bcf > $LOCAL1/final-assemblies.fasta

# extract the final reads from each polsihed assembly
# input: $LOCAL1/final-assemblies.fasta
# output: $LOCAL/$name/$name-final-assembly.fasta (same number of files as contigs defined in $MAGS)
./extract-circularized-fasta.py $LOCAL1/final-assemblies.fasta $LOCAL

echo Finished all validation and polishing


















###################################

# orient genome assemblies

###################################

for name in $MAGS
do

#################
# these steps are done to generate improved nanopore mapping to the final assembly

# calculating GC content and GC skew

# from the fasta file itself, use seqinr to calculate the GC content in 1000 base windows that can be added to the files in R

# activate conda virtual environment to run prokka
source /Volumes/data/bin/miniconda3/bin/activate

# run prokka on all genomes, make sure you have the most recent version of tbl2asn, or this will not work
# input: $name $LOCAL/$name/$name-final-assembly.fasta
# output (multiple files; relevant ones provided) $LOCAL/$name/circos/prokka-initial/$name.gff
prokka --outdir $LOCAL/$name/circos/prokka --prefix $name $LOCAL/$name/$name-final-assembly.fasta --cpu $THREADS

# re-orient the genome around the dnaA gene where there is a switch in GC-skew
# input: $LOCAL/$name/circos/prokka-initial/$name.gff
# output: $LOCAL/$name/circos/dnaA-loci.txt
grep 'dnaA\|DnaA\|dnaa' $LOCAL/$name/circos/prokka/$name.gff > $LOCAL/$name/circos/dnaA-loci.txt

# input: $LOCAL/$name/circos/dnaA-loci.txt
# input: $LOCAL/$name/$name-final-assembly.fasta
# output: $LOCAL/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL/$name/$name-orientation-shift.txt
# requires: orient-v2.py
./orient-v2.py $LOCAL/$name/circos/dnaA-loci.txt $LOCAL/$name/$name-final-assembly.fasta $LOCAL/$name/$name-final-assembly-oriented.fasta $name $LOCAL/$name/$name-orientation-shift.txt

# get the gc content, gc skew and culmulative gc skew in 1000 base windows
# input: $LOCAL/$name/$name-final-assembly.fasta
# output: $LOCAL/$name/$name-gc-info.txt
# requires: circlize_gc_information.R
./circlize_gc_information.R -i $LOCAL/$name/$name-final-assembly-oriented.fasta -o $LOCAL/$name/circos/$name-gc-info.txt

# rebuild the final assembly to have 80kb of the beginning of the genome copied to the end of the genome to improve remapping of long reads
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL/$name/$name-final-assembly-reverse-start.fasta
# requires: modify-genome-v4.py
./modify-genome-v4.py $LOCAL/$name/$name-final-assembly-oriented.fasta $LOCAL/$name/$name-final-assembly-reverse-start.fasta

# group all genomes oriented and their reverse pair
# input: $LOCAL/$name/$name-final-assembly-reverse-start.fasta
# output: $LOCAL1/all-reverse-start.fasta
cat $LOCAL/$name/$name-final-assembly-reverse-start.fasta >> $LOCAL1/all-reverse-start.fasta
# group all oriented genomes
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL1/all-oriented.fasta
cat $LOCAL/$name/$name-final-assembly-oriented.fasta >> $LOCAL1/all-oriented.fasta

done

# progress update
echo Finished genome orientation















###################################

# Remap the reads to the oriented and reverse-oriented genomes
# Filter according to the previous filtering parameters
# Concatenate the filtered oriented and reverse-oriented mapped reads and convert to fastq for use in reassembly

###################################


# remap unfiltered nanopore reads against the final forward and reverse assemblies
# input: $LOCAL1/all-reverse-start.fasta
# input: $LONGREADS
# output: $LOCAL1/nanopore-unfilteredr-to-final.sam
$MAP -t $THREADS -aLQx map-ont --secondary=no --sam-hit-only $LOCAL1/all-reverse-start.fasta $LONGREADS > $LOCAL1/nanopore-unfilteredr-to-final.sam

# remap unfiltered nanopore reads against the final forward and reverse assemblies
# input: $LOCAL1/all-oriented.fasta
# input: $LONGREADS
# output: $LOCAL1/nanopore-unfilteredf-to-final.sam
$MAP -t $THREADS -aLQx map-ont --secondary=no --sam-hit-only $LOCAL1/all-oriented.fasta $LONGREADS > $LOCAL1/nanopore-unfilteredf-to-final.sam

# extract relevant information on nanopore reverse sequence alignments scores
# input: $LOCAL1/nanopore-unfilteredr-to-final.sam
# output: $LOCAL1/sam.txt
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $LOCAL1/nanopore-unfilteredr-to-final.sam > $LOCAL1/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
# input: $LOCAL1/cigars.txt
# output: $LOCAL1/cigar-results.txt
# requires: cigar-parse.py
./cigar-parse.py $LOCAL1/sam.txt $LOCAL1/final-read-names.txt

# filter the reverse nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# this gets the sam headers in the file
# input: $LOCAL1/nanopore-unfilteredr-to-final.sam
# output: $LOCAL1/nanopore-filteredr-to-final.sam
awk '$0 ~ "@"{print $0}' $LOCAL1/nanopore-unfilteredr-to-final.sam > $LOCAL1/nanopore-filteredr-to-final.sam
# input: $LOCAL1/final-read-names.txt
# input: $LOCAL1/nanopore-unfilteredr-to-final.sam
# output: $LOCAL1/nanopore-filteredr-to-final.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $LOCAL1/final-read-names.txt $LOCAL1/nanopore-unfilteredr-to-final.sam >> $LOCAL1/nanopore-filteredr-to-final.sam

# extract relevant information on nanopore forward sequence alignments scores
# input: $LOCAL1/nanopore-unfilteredf-to-final.sam
# output: $LOCAL1/sam.txt
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $LOCAL1/nanopore-unfilteredf-to-final.sam > $LOCAL1/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
# input: $LOCAL1/cigars.txt
# output: $LOCAL1/cigar-results.txt
# requires: cigar-parse.py
./cigar-parse.py $LOCAL1/sam.txt $LOCAL1/final-read-names.txt

# filter the reverse nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# this gets the sam headers in the file
# input: $LOCAL1/nanopore-unfilteredf-to-final.sam
# output: $LOCAL1/nanopore-filteredf-to-final.sam
awk '$0 ~ "@"{print $0}' $LOCAL1/nanopore-unfilteredf-to-final.sam > $LOCAL1/nanopore-filteredf-to-final.sam
# input: $LOCAL1/final-read-names.txt
# input: $LOCAL1/nanopore-unfilteredf-to-final.sam
# output: $LOCAL1/nanopore-filteredf-to-final.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $LOCAL1/final-read-names.txt $LOCAL1/nanopore-unfilteredf-to-final.sam >> $LOCAL1/nanopore-filteredf-to-final.sam

# make a bowtie2 database for mapping illumina reads to the finally assembly
# input: $LOCAL1/all-oriented.fasta
# output: $LOCAL1/genomes-final-assembly
$BUILD $LOCAL1/all-oriented.fasta $LOCAL1/genomes-final-assembly

# map the illumina reads against the final assembly
# input: $LOCAL1/genomes-final-assembly
# input: $FWD
# input: $REV
# output: $LOCAL1/illumina-to-final.sam
$BOWTIE -x $LOCAL1/genomes-final-assembly -1 $FWD -2 $REV -p $THREADS -S $LOCAL1/illumina-to-final.sam --no-unal --no-discordant --end-to-end --very-sensitive



# sort the unfiltered-forward nanopore reads
# input: $LOCAL1/nanopore-unfilteredf-to-final.sam
# output: $LOCAL1/nanopore-unfilteredf-to-final-sorted.bam
$SAMTOOLS sort $LOCAL1/nanopore-unfilteredf-to-final.sam -o $LOCAL1/nanopore-unfilteredf-to-final-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-unfilteredf-to-final-sorted.bam
# output: $LOCAL1/nanopore-unfilteredf-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/nanopore-unfilteredf-to-final-sorted.bam -@ $THREADS

# sort the unfiltered-reverse nanopore reads
# input: $LOCAL1/nanopore-unfilteredr-to-final.sam
# output: $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam
$SAMTOOLS sort $LOCAL1/nanopore-unfilteredr-to-final.sam -o $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam
# output: $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam -@ $THREADS



# sort the filtered-forward nanopore reads
# input: $LOCAL1/nanopore-filteredf-to-final.sam
# output: $LOCAL1/nanopore-filtered-to-final-sorted.bam
$SAMTOOLS sort $LOCAL1/nanopore-filteredf-to-final.sam -o $LOCAL1/nanopore-filteredf-to-final-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-filteredf-to-final-sorted.bam
# output: $LOCAL1/nanopore-filteredf-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/nanopore-filteredf-to-final-sorted.bam -@ $THREADS

# sort the filtered-reverse nanopore reads
# input: $LOCAL1/nanopore-filteredr-to-final.sam
# output: $LOCAL1/nanopore-filteredr-to-final-sorted.bam
$SAMTOOLS sort $LOCAL1/nanopore-filteredr-to-final.sam -o $LOCAL1/nanopore-filteredr-to-final-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-filteredr-to-final-sorted.bam
# output: $LOCAL1/nanopore-filteredr-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/nanopore-filteredr-to-final-sorted.bam -@ $THREADS




# create a fasta index for proper header generation after extracting specific regional reads
# input: $LOCAL1/all-reverse-start.fasta
# output: $LOCAL1/all-reverse-start.fasta.fai
$SAMTOOLS faidx $LOCAL1/all-reverse-start.fasta

# create a fasta index for proper header generation after extracting specific regional reads
# input: $LOCAL1/all-oriented.fasta
# output: $LOCAL1/all-oriented.fasta.fai
$SAMTOOLS faidx $LOCAL1/all-oriented.fasta



# extract the "soft" mapped reads
# input: $LOCAL1/illumina-to-final.sam
# output: $LOCAL1/illumina-to-final-soft.bam
$SAMTOOLS view -@ $THREADS -f 2 -F 256 -bS $LOCAL1/illumina-to-final.sam > $LOCAL1/illumina-to-final-soft.bam

# sort the extracted soft reads
# input: $LOCAL1/illumina-to-final-soft.bam
# output: $LOCAL1/illumina-to-final-soft-sorted.bam
$SAMTOOLS sort $LOCAL1/illumina-to-final-soft.bam -o $LOCAL1/illumina-to-final-soft-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-filtered-to-final-sorted.bam
# output: $LOCAL1/nanopore-filtered-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/illumina-to-final-soft-sorted.bam -@ $THREADS

# use more stringent filtering parameters from the illumina alignment against the final assembly in preparation for depth determining
# input: $LOCAL1/illumina-to-final.sam
# output: $LOCAL1/illumina-to-final-hard.bam 
$SAMTOOLS view -@ $THREADS -bS -f 2 -F 3848 $LOCAL1/illumina-to-final.sam > $LOCAL1/illumina-to-final-hard.bam

# sort the extracted hard reads
# input: $LOCAL1/illumina-to-final-hard.bam
# output: $LOCAL1/illumina-to-final-hard-sorted.bam
$SAMTOOLS sort $LOCAL1/illumina-to-final-hard.bam -o $LOCAL1/illumina-to-final-hard-sorted.bam -@ $THREADS

# index the sorted reads
# input: $LOCAL1/nanopore-filtered-to-final-hard-sorted.bam
# output: $LOCAL1/nanopore-filtered-to-final-sorted.bam.bai
$SAMTOOLS index $LOCAL1/illumina-to-final-hard-sorted.bam -@ $THREADS

# create a fasta index for proper header generation after extracting specific regional reads
# input: $LOCAL1/all-oriented.fasta
# output: $LOCAL1/all-oriented.fasta.fai
$SAMTOOLS faidx $LOCAL1/all-oriented.fasta

for name in $MAGS
do

# collect the relevant illumina reads
# input: $LOCAL1/illumina-to-final-soft-sorted.bam
# input: $LOCAL1/illumina-to-final-soft-sorted.bam.bai
# input: $LOCAL1/all-oriented.fasta.fai
# output: $LOCAL/$name/illumina-to-final-soft.bam
$SAMTOOLS view $LOCAL1/illumina-to-final-soft-sorted.bam $name -o $LOCAL/$name/illumina-to-final-soft.bam -T $LOCAL1/all-oriented.fasta -@ $THREADS


# input: $LOCAL1/nanopore-filteredf-to-final-sorted.bam
# input: $LOCAL1/nanopore-filteredf-to-final-sorted.bam.bai
# input: $LOCAL1/all-oriented.fasta.fai
# output: $LOCAL/$name/nanopore-filtered-to-final-forward.bam
$SAMTOOLS view $LOCAL1/nanopore-filteredf-to-final-sorted.bam $name -o $LOCAL/$name/nanopore-filtered-to-final-forward.bam -T $LOCAL1/all-oriented.fasta -@ $THREADS
# input: $LOCAL1/nanopore-filteredr-to-final-sorted.bam
# input: $LOCAL1/nanopore-filteredr-to-final-sorted.bam.bai
# input: $LOCAL1/all-reverse-start.fasta.fai
# output: $LOCAL/$name/nanopore-filtered-to-final-reverse.bam
$SAMTOOLS view $LOCAL1/nanopore-filteredr-to-final-sorted.bam $name-reversed -o $LOCAL/$name/nanopore-filtered-to-final-reverse.bam -T $LOCAL1/all-reverse-start.fasta -@ $THREADS


# sort the reads by position
# input: $LOCAL/$name/nanopore-filtered-to-final-forward.bam
# output: $LOCAL/$name/nanopore-filtered-to-final-forward-sorted.bam
$SAMTOOLS sort -@ $THREADS $LOCAL/$name/nanopore-filtered-to-final-forward.bam -o $LOCAL/$name/nanopore-filtered-to-final-forward-sorted.bam
# input: $LOCAL/$name/nanopore-filtered-to-final-reverse.bam
# output: $LOCAL/$name/nanopore-filtered-to-final-reverse-sorted.bam
$SAMTOOLS sort -@ $THREADS $LOCAL/$name/nanopore-filtered-to-final-reverse.bam -o $LOCAL/$name/nanopore-filtered-to-final-reverse-sorted.bam


# sort the reads by name
# input: $LOCAL/$name/illumina-to-final-soft.bam
# output: $LOCAL/$name/illumina-to-final-soft-sorted-by-name.bam
$SAMTOOLS sort -n -@ $THREADS $LOCAL/$name/illumina-to-final-soft.bam -o $LOCAL/$name/illumina-to-final-soft-sorted-by-name.bam




# convert the concatenated filtered bam file to fastq
# input: $LOCAL/$name/nanopore-filtered-to-final-forward-sorted.bam
# output: $LOCAL/$name/nanopore-filtered-to-final-forward.fastq
$SAMTOOLS fastq $LOCAL/$name/nanopore-filtered-to-final-forward-sorted.bam > $LOCAL/$name/nanopore-filtered-to-final-forward.fastq -@ $THREADS
# input: $LOCAL/$name/nanopore-filtered-to-final-reverse-sorted.bam
# output: $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq
$SAMTOOLS fastq $LOCAL/$name/nanopore-filtered-to-final-reverse-sorted.bam > $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq -@ $THREADS




# concatenate the forward and reverse nanopore fastq files
# input: $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq
# output: $LOCAL/$name/nanopore-filtered-to-final.fastq
cat $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq > $LOCAL/$name/nanopore-filtered-to-final.fastq

# get the read names mapping
# input: $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq
# output: $LOCAL/$name/unique-read-names.txt
grep "^@" $LOCAL/$name/nanopore-filtered-to-final-reverse.fastq > $LOCAL/$name/unique-read-names.txt

# concatenate reads that are not found in the fastq already to the fastq
# input: $LOCAL/$name/unique-read-names.txt
# input: $LOCAL/$name/nanopore-filtered-to-final-forward.fastq
# output: $LOCAL/$name/nanopore-filtered-to-final.fastq
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $LOCAL/$name/unique-read-names.txt $LOCAL/$name/nanopore-filtered-to-final-forward.fastq | grep -v -A 3 -f - >> $LOCAL/$name/nanopore-filtered-to-final.fastq



# separate the reads into forward and reverse types (could have used bam2fastq, still required a sorted bam file)
# input: $LOCAL/$name/illumina-to-final-soft-sorted-by-name.bam
# output: $LOCAL/$name/illumina-soft-paired1.fq
# output: $LOCAL/$name/illumina-soft-paired2.fq
$SAMTOOLS fastq -1 $LOCAL/$name/illumina-soft-paired1.fq -2 $LOCAL/$name/illumina-soft-paired2.fq -0 /dev/null -s /dev/null -n $LOCAL/$name/illumina-to-final-soft-sorted-by-name.bam -@ $THREADS



done 


echo Finished final mapping














###################################

# begin reassembly attempts

###################################

## Re-assembly and analysis using different tools

for name in $MAGS
do

source /Volumes/data/bin/miniconda3/bin/activate

# reassemble using flye not in meta mode with the subset of nanopore reads filtered from the mapping against the final assembly
# input: $LOCAL/$name/nanopore-filtered-to-final.fastq
# output dir: $LOCAL2/$name/flye/ has multiple files
flye --nano-raw $LOCAL/$name/nanopore-filtered-to-final.fastq --threads $THREADS -g 5m -o $LOCAL2/$name/flye/

# move the contigs fasta file and assembly info file to up two directories
# input: $LOCAL2/$name/flye/assembly.fasta
# output: $LOCAL2/$name-flye-assembly.fasta
mv $LOCAL2/$name/flye/assembly.fasta $LOCAL2/$name-flye-assembly.fasta
# input: $LOCAL2/$name/flye/assembly_info.txt
# output: $LOCAL2/$name-flye.log
mv $LOCAL2/$name/flye/assembly_info.txt $LOCAL2/$name-flye.log

conda deactivate

# assemble with unicycler hybrid assembly
# input: $LOCAL/$name/nanopore-filtered-to-final.fastq
# input: $LOCAL/$name/illumina-soft-paired1.fq
# input: $LOCAL/$name/illumina-soft-paired2.fq
# output dir: $LOCAL2/$name/unicycler/ multiple unicycler files
/usr/bin/python3.6 $UNICYCLER -t $THREADS -1 $LOCAL/$name/illumina-soft-paired1.fq -2 $LOCAL/$name/illumina-soft-paired2.fq -l $LOCAL/$name/nanopore-filtered-to-final.fastq -o $LOCAL2/$name/unicycler/ --no_pilon

# move the relevant files up two directories
# input: $LOCAL2/$name/unicycler/assembly.fasta
# output: $LOCAL2/$name-unicycler-assembly.fasta
mv $LOCAL2/$name/unicycler/assembly.fasta $LOCAL2/$name-unicycler-assembly.fasta
# input: $LOCAL2/$name/unicycler/unicycler.log
# output: $LOCAL2/$name-unicycler.log
mv $LOCAL2/$name/unicycler/unicycler.log $LOCAL2/$name-unicycler.log

# reassemble using unicycler long-read assembly 
# input: $LOCAL/$name/nanopore-filtered-to-final.fastq
# output dir: $LOCAL2/$name/unicycler-long multiple unicycler long read files
/usr/bin/python3.6 $UNICYCLER -t $THREADS -l $LOCAL/$name/nanopore-filtered-to-final.fastq -o $LOCAL2/$name/unicycler-long --no_pilon

# move the relevant files up two directories
# input: $LOCAL2/$name/unicycler-long/assembly.fasta
# output: $LOCAL2/$name-unicycler-long-assembly.fasta
mv $LOCAL2/$name/unicycler-long/assembly.fasta $LOCAL2/$name-unicycler-long-assembly.fasta
# input: $LOCAL2/$name/unicycler-long/unicycler.log
# output: $LOCAL2/$name-unicycler-long.log
mv $LOCAL2/$name/unicycler-long/unicycler.log $LOCAL2/$name-unicycler-long.log

# assemble with spades
# input: $LOCAL/$name/$name-illumina-soft-paired1.fq
# input: $LOCAL/$name/$name-illumina-soft-paired2.fq
# output dir: $LOCAL2/$name/spades multiple spades assembly files
$SPADES -1 $LOCAL/$name/illumina-soft-paired1.fq -2 $LOCAL/$name/illumina-soft-paired2.fq --careful -t $THREADS -o $LOCAL2/$name/spades 

# move the relevant information up two directories
# input: $LOCAL2/$name/spades/contigs.fasta
# output: $LOCAL2/$name-spades-assembly.fasta 
mv $LOCAL2/$name/spades/contigs.fasta $LOCAL2/$name-spades-assembly.fasta 
# input: $LOCAL2/$name/spades/spades.log
# output: $LOCAL2/$name-spades.log
mv $LOCAL2/$name/spades/spades.log $LOCAL2/$name-spades.log


done

echo Finished all re-assemblies 














###################################

# create dot plots for the flye re-assemblies against the validated and corrected assemblies to ensure they are assembled correctly

###################################

for name in $MAGS
do

# create dot plot(s) for the given genome
# input: $LOCAL2/$name-flye-assembly.fasta
# input: $LOCAL/$name/$name-final-assembly.fasta
# output: 
# output dir: $LOCAL4/$name has various files (not used again)
# requires: dot-plot2.py
./dot-plot2.py $name $LOCAL4 $SMALLPLOTS $LOCAL2/$name-flye-assembly.fasta $LOCAL/$name/$name-final-assembly.fasta

mv $LOCAL4/$name $LOCAL/output/dotplots

done

echo Dot plots made



















###################################

# generate preliminary files for use in figure (circos_plot.Rmd)

###################################



source /Volumes/data/bin/miniconda3/bin/activate

# calculate the depth over a 1000 base window using mapped unfiltered and filtered for reverse reads for nanopore
# input: $LOCAL1/$name/nanopore-filtered-to-final-forward-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/nanopore-filteredf1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/nanopore-filteredf1 $LOCAL1/nanopore-filteredf-to-final-sorted.bam
# input: $LOCAL/$name/nanopore-unfiltered-to-final-forward-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/nanopore-unfilteredf1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/nanopore-unfilteredf1 $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam
# input: $LOCAL1/nanopore-filteredr-to-final-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/nanopore-filteredr1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/nanopore-filteredr1 $LOCAL1/nanopore-filteredr-to-final-sorted.bam
# input: $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/nanopore-unfilteredr1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/nanopore-unfilteredr1 $LOCAL1/nanopore-unfilteredr-to-final-sorted.bam

# calculate the depth over a 1000 base window using mapped illumina reads
# input: $LOCAL1/illumina-to-final-soft-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/illumina-soft1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/illumina-soft1 $LOCAL1/illumina-to-final-soft-sorted.bam
# input: $LOCAL1/illumina-to-final-hard-sorted.bam.bai
# output (multiple files; relevant ones provided): $LOCAL1/illumina-hard1.regions.bed.gz
mosdepth -b 1000 $LOCAL1/illumina-hard1 $LOCAL1/illumina-to-final-hard-sorted.bam

# it outputs in a gzipped bed file, so unzip
# input: $LOCAL1/nanopore-filteredf1.regions.bed.gz
# output: $LOCAL1/nanopore-filteredf1.regions.bed
gunzip $LOCAL1/nanopore-filteredf1.regions.bed.gz
# input: $LOCAL1/nanopore-unfilteredf1.regions.bed.gz
# output: $LOCAL1/nanopore-unfilteredf1.regions.bed
gunzip $LOCAL1/nanopore-unfilteredf1.regions.bed.gz
# input: $LOCAL1/nanopore-filteredr1.regions.bed.gz
# output: $LOCAL1/nanopore-filteredr1.regions.bed
gunzip $LOCAL1/nanopore-filteredr1.regions.bed.gz
# input: $LOCAL1/nanopore-unfilteredr1.regions.bed.gz
# output: $LOCAL1/nanopore-unfilteredr1.regions.bed
gunzip $LOCAL1/nanopore-unfilteredr1.regions.bed.gz
# input: $LOCAL1/illumina-soft1.regions.bed.gz
# output: $LOCAL1/illumina-soft1.regions.bed
gunzip $LOCAL1/illumina-soft1.regions.bed.gz
# input: $LOCAL1/illumina-hard1.regions.bed.gz
# output: $LOCAL1/illumina-hard1.regions.bed
gunzip $LOCAL1/illumina-hard1.regions.bed.gz


# apply to each genome
for name in $MAGS
do

# collect only the relevant genome
# input: $LOCAL1/nanopore-filteredf1.regions.bed
# output: $LOCAL/$name/circos/$name-nanopore-filteredf.regions.bed
grep $name $LOCAL1/nanopore-filteredf1.regions.bed > $LOCAL/$name/circos/$name-nanopore-filteredf.regions.bed
# input: $LOCAL1/nanopore-unfilteredf1.regions.bed
# output: $LOCAL/$name/circos/$name-nanopore-unfilteredf.regions.bed
grep $name $LOCAL1/nanopore-unfilteredf1.regions.bed > $LOCAL/$name/circos/$name-nanopore-unfilteredf.regions.bed
# input: $LOCAL1/nanopore-filteredr1.regions.bed
# output: $LOCAL/$name/circos/$name-nanopore-filteredr.regions.bed
grep $name $LOCAL1/nanopore-filteredr1.regions.bed > $LOCAL/$name/circos/$name-nanopore-filteredr.regions.bed
# input: $LOCAL1/nanopore-unfilteredr1.regions.bed
# output: $LOCAL/$name/circos/$name-nanopore-unfilteredr.regions.bed
grep $name $LOCAL1/nanopore-unfilteredr1.regions.bed > $LOCAL/$name/circos/$name-nanopore-unfilteredr.regions.bed

# input: $LOCAL1/illumina-soft1.regions.bed
# output: $LOCAL/$name/circos/$name-illumina-soft.regions.bed
grep $name $LOCAL1/illumina-soft1.regions.bed > $LOCAL/$name/circos/$name-illumina-soft.regions.bed
# input: $LOCAL1/illumina-hard1.regions.bed
# output: $LOCAL/$name/circos/$name-illumina-hard.regions.bed
grep $name $LOCAL1/illumina-hard1.regions.bed > $LOCAL/$name/circos/$name-illumina-hard.regions.bed

# modify the unfiltered.regions.bed and filtered.regions.bed files to include accurate measurements
# input: $LOCAL/$name/circos/$name-nanopore-filteredf.regions.bed
# input: $LOCAL/$name/circos/$name-nanopore-unfilteredf.regions.bed
# input: $LOCAL/$name/circos/$name-nanopore-filteredr.regions.bed
# input: $LOCAL/$name/circos/$name-nanopore-unfilteredr.regions.bed
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL/$name/circos/$name-nanopore-unfiltered.regions.bed
# output: $LOCAL/$name/circos/$name-nanopore-filtered.regions.bed
# requires: modify-nanopore-bed-v2.py
./modify-nanopore-bed-v2.py $LOCAL/$name/circos/$name-nanopore-filtered.regions.bed $LOCAL/$name/circos/$name-nanopore-unfiltered.regions.bed $LOCAL/$name/circos/$name-nanopore-filteredf.regions.bed $LOCAL/$name/circos/$name-nanopore-unfilteredf.regions.bed $LOCAL/$name/circos/$name-nanopore-filteredr.regions.bed $LOCAL/$name/circos/$name-nanopore-unfilteredr.regions.bed $LOCAL/$name/$name-final-assembly-oriented.fasta $name

# generate length of genome
# input: $LOCAL/$name/$name-final-assembly.fasta
# output: $LOCAL/$name/circos/$name-final-length.txt
awk '/^>/{if (l!="") l=0; next}{l+=length($0)}END{print l}' $LOCAL/$name/$name-final-assembly-oriented.fasta > $LOCAL/$name/circos/$name-final-length.txt

# get length of genome
# input: $LOCAL/$name/circos/$name-final-length.txt
# output: $length
length=`cat $LOCAL/$name/circos/$name-final-length.txt`
# generate cytoband file needed for circlize
# input: $length
# output: $LOCAL/$name/circos/$name-cytoband.txt
echo -e "$name\t0\t$length\tnot_real\tnot_real" > $LOCAL/$name/circos/$name-cytoband.txt

# generate the required bed files from the .gff output.

# first two lines are headers
# i'm printing it into bed format, then printing "1" as a value file.

# generate header for BED files for positive strand, negative strand, tRNA and rRNA
# output: $LOCAL/$name/circos/$name-cds-positive.bed
# output: $LOCAL/$name/circos/$name-cds-negative.bed
# output: $LOCAL/$name/circos/$name-cds-trna.bed
# output: $LOCAL/$name/circos/$name-cds-rrna.bed
echo -e "chr\tstart\tend\tcds_pos" > $LOCAL/$name/circos/$name-cds-positive.bed
echo -e "chr\tstart\tend\tcds_neg" > $LOCAL/$name/circos/$name-cds-negative.bed
echo -e "chr\tstart\tend\ttrna" > $LOCAL/$name/circos/$name-cds-trna.bed
echo -e "chr\tstart\tend\trrna" > $LOCAL/$name/circos/$name-cds-rrna.bed

# use tab separation (OFS arg)
# coding sequences
# input: $LOCAL/$name/circos/prokka/$name.gff
# output: $LOCAL/$name/circos/$name-cds-positive.bed
awk -v OFS='\t' '$3 ~ "CDS" && $7 ~ "+" {print $1, $4, $5, "1"}' $LOCAL/$name/circos/prokka/$name.gff >> $LOCAL/$name/circos/$name-cds-positive.bed
# input: $LOCAL/$name/circos/prokka/$name.gff
# output: $LOCAL/$name/circos/$name-cds-negative.bed
awk -v OFS='\t' '$3 ~ "CDS" && $7 ~ "-" {print $1, $4, $5, "1"}' $LOCAL/$name/circos/prokka/$name.gff >> $LOCAL/$name/circos/$name-cds-negative.bed

# get tRNA and rRNA gene loci
# input: $LOCAL/$name/circos/prokka/$name.gff
# output: $name/${name}-cds-trna.bed
awk -v OFS='\t' '$3 ~ "tRNA" {print $1, $4, $5, "1"}' $LOCAL/$name/circos/prokka/$name.gff >> $LOCAL/$name/circos/$name-cds-trna.bed
# input: $LOCAL/$name/circos/prokka/$name.gff
# output: $LOCAL/$name/circos/$name-cds-rrna.bed
awk -v OFS='\t' '$3 ~ "rRNA" {print $1, $4, $5, "1"}' $LOCAL/$name/circos/prokka/$name.gff >> $LOCAL/$name/circos/$name-cds-rrna.bed

# modify the bed files based on the positions of dnaa gene rarrangement and size of the genome
# input: $LOCAL/$name/$name-orientation-shift.txt
# input: $LOCAL/$name/circos/$name-cds-trna.bed
# input: $LOCAL/$name/circos/$name-cds-rrna.bed
# input: $LOCAL/$name/circos/$name-cds-positive.bed
# input: $LOCAL/$name/circos/$name-cds-negative.bed
# output: $LOCAL/$name/$name-orientation-shift.txt
# output: $LOCAL/$name/circos/$name-cds-trna.bed
# output: $LOCAL/$name/circos/$name-cds-rrna.bed
# output: $LOCAL/$name/circos/$name-cds-positive.bed
# output: $LOCAL/$name/circos/$name-cds-negative.bed
# requires: bed-file-orientation.py
./bed-file-orientation.py $LOCAL/$name/$name-orientation-shift.txt $LOCAL/$name/circos/$name-cds-trna.bed $LOCAL/$name/circos/$name-cds-rrna.bed $LOCAL/$name/circos/$name-cds-positive.bed $LOCAL/$name/circos/$name-cds-negative.bed

mv $LOCAL/$name/circos/$name-cds* $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-gc* $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-cytoband.txt $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-illumina-hard.regions.bed $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-illumina-soft.regions.bed $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-nanopore-filtered.regions.bed $LOCAL/output/circos
mv $LOCAL/$name/circos/$name-nanopore-unfiltered.regions.bed $LOCAL/output/circos

done





###################################

# rough annotations and genome investigations

###################################

for name in $MAGS
do

# create a file for illumina and nanopore total bp coverages
# output: $LOCAL3/$name/$name-coverage.txt
echo "$name average genome coverages" > $LOCAL3/$name/output/$name-coverage.txt

# get illumina and nanopore total bp coverages
# input: $LOCAL/$name/circos/$name-illumina-hard.regions.bed
# output: $LOCAL3/$name/output/$name-coverage.txt
awk '{sum+=$4} END { print "Illumina coverage\t",sum}' $LOCAL/output/circos/$name-illumina-hard.regions.bed >> $LOCAL3/$name/output/$name-coverage.txt
# input: $LOCAL/$name/circos/$name-nanopore-filtered.regions.bed
# output: $LOCAL3/$name/output/$name-coverage.txt
awk '{sum+=$4} END { print "Nanopore coverage\t",sum}' $LOCAL/output/circos/$name-nanopore-filtered.regions.bed >> $LOCAL3/$name/output/$name-coverage.txt

# determine the number of tRNA in the genome
# input: $LOCAL/$name/circos/prokka/$name.gff
# output: $LOCAL3/$name/output/tRNA-number.txt
awk -v OFS='\t' '$3 ~ "tRNA" {print $1, $4, $5, "1"}' $LOCAL/$name/circos/prokka/$name.gff | wc -l > $LOCAL3/$name/output/tRNA-number.txt

# start anvio6 virtual environment
source /Volumes/data/virtualenvs/anvio-6/bin/activate

# copy the genome for annotating
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL3/$name/$name-final-assembly-oriented.fasta
cp $LOCAL/$name/$name-final-assembly-oriented.fasta $LOCAL3/$name/$name-final-assembly-oriented.fasta

# generate a contigs database of the final genome assembly using anvio
# input: $LOCAL3/$name/$name-final-assembly-oriented.fasta
# output: $LOCAL3/$name/$name-annotation.db
anvi-gen-contigs-database -f $LOCAL3/$name/$name-final-assembly-oriented.fasta -n $name-validation-annotation -o $LOCAL3/$name/$name-annotation.db

##############
##THE FOLLOWING FUNCTIONS MUST BE RUN IN ORDER
# run hmms on the database (no longer automatic in anvio-6)
# input: $LOCAL/$name/$name-annotation.db
anvi-run-hmms -T $THREADS -c $LOCAL3/$name/$name-annotation.db

# affiliate single copy gene taxonomies 
# input: $LOCAL/$name/$name-annotation.db
anvi-run-scg-taxonomy -T $THREADS -c $LOCAL3/$name/$name-annotation.db

# estimate genome taxonomy 
# input: $LOCAL/$name/$name-annotation.db
# output: $LOCAL3/$name/$name-taxonomies.txt
anvi-estimate-genome-taxonomy -T $THREADS --just-do-it -c $LOCAL3/$name/$name-annotation.db -o $LOCAL3/$name/$name-taxonomies.txt

# estimate genomes completeness
# input: $LOCAL/$name/$name-annotation.db
# output: $LOCAL3/$name/$name-completeness.txt
anvi-estimate-genome-completeness -c $LOCAL3/$name/$name-annotation.db -o $LOCAL3/$name/$name-completeness.txt

# get contig stats for the contig
# input: $LOCAL/$name/$name-annotation.db
# output: $LOCAL3/$name/contigs-stats.txt
anvi-display-contigs-stats $LOCAL3/$name/$name-annotation.db --report-as-text -o $LOCAL3/$name/contigs-stats.txt

done

echo Finished all anvio annotations








###################################

# create the summary TABLE1 for the assemblies

###################################

# create column titles for the TABLE1
# output: $TABLE1
echo "| Contig Number | Predicted Taxonomy | Total Length | GC Content | % Completion | % Redundancy | Illumina Coverage | Nanopore (F) Coverage | Number of rRNAs | Number of tRNAs |" > $TABLE1
# input: $TABLE1
# output: $TABLE1
echo "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|" >> $TABLE1

# create column titles for TABLE2
# output: 
echo "| Predicted Taxonomy | Flye Reassembled | Unicycler Hybrid Reassembled | Miniasm Reassembled | Spades Reassembled |" > $TABLE2
# input: $TABLE2
# output: $TABLE2
echo "|:---:|:---:|:---:|:---:|:---:|" >> $TABLE2

# get information from each genome analysis
for name in $MAGS
do 

# use python script to populate TABLE1
# input: $TABLE1
# input: $LOCAL3/$name/$name-taxonomies.txt
# input: $LOCAL3/$name/contigs-stats.txt
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# input: $LOCAL3/$name/$name-completeness.txt
# input: $LOCAL3/$name/output/$name-coverage.txt
# input: $LOCAL3/$name/output/tRNA-number.txt
# input: $LOCAL2/$name-flye-assembly.fasta
# input: $LOCAL2/$name-flye.log
# input: $LOCAL2/$name-unicycler-assembly.fasta
# input: $LOCAL2/$name-unicycler.log
# input: $LOCAL2/$name-unicycler-long-assembly.fasta
# input: $LOCAL2/$name-unicycler-long.log
# input: $LOCAL2/$name-spades-assembly.fasta 
# output: $TABLE1
# requires: table1.py
./table1.py $name $TABLE1 $LOCAL3/$name/$name-taxonomies.txt $LOCAL3/$name/contigs-stats.txt $LOCAL/$name/$name-final-assembly-oriented.fasta $LOCAL3/$name/$name-completeness.txt $LOCAL3/$name/output/$name-coverage.txt $LOCAL3/$name/output/tRNA-number.txt $LOCAL2/$name-flye-assembly.fasta $LOCAL2/$name-flye.log $LOCAL2/$name-unicycler-assembly.fasta $LOCAL2/$name-unicycler.log $LOCAL2/$name-unicycler-long-assembly.fasta $LOCAL2/$name-unicycler-long.log $LOCAL2/$name-spades-assembly.fasta 


# use python script to populate TABLE2
# input: $TABLE1
# input: $LOCAL3/$name/$name-taxonomies.txt
# input: $LOCAL3/$name/contigs-stats.txt
# input: $LOCAL/$name/$name-final-assembly-oriented.fasta
# input: $LOCAL3/$name/$name-completeness.txt
# input: $LOCAL3/$name/output/$name-coverage.txt
# input: $LOCAL3/$name/output/tRNA-number.txt
# input: $LOCAL2/$name-flye-assembly.fasta
# input: $LOCAL2/$name-flye.log
# input: $LOCAL2/$name-unicycler-assembly.fasta
# input: $LOCAL2/$name-unicycler.log
# input: $LOCAL2/$name-unicycler-long-assembly.fasta
# input: $LOCAL2/$name-unicycler-long.log
# input: $LOCAL2/$name-spades-assembly.fasta 
# output: $TABLE2
# requires: table2.py
./table2.py $name $TABLE2 $LOCAL3/$name/$name-taxonomies.txt $LOCAL3/$name/contigs-stats.txt $LOCAL/$name/$name-final-assembly-oriented.fasta $LOCAL3/$name/$name-completeness.txt $LOCAL3/$name/output/$name-coverage.txt $LOCAL3/$name/output/tRNA-number.txt $LOCAL2/$name-flye-assembly.fasta $LOCAL2/$name-flye.log $LOCAL2/$name-unicycler-assembly.fasta $LOCAL2/$name-unicycler.log $LOCAL2/$name-unicycler-long-assembly.fasta $LOCAL2/$name-unicycler-long.log $LOCAL2/$name-spades-assembly.fasta 

done

# progress update
echo Finished Table 1 and Table 2






###################################

# determine proportion of reads allocated to each genome

###################################

# make the table file
# output: $TABLE3
echo "| Contig Number | Percent of all Nanopore Reads Mapping | Percent of all Illumina Reads Mapping |" > $TABLE3
echo "|:---:|:---:|:---:|" >> $TABLE3

# assign total nanopore reads and total illumina reads variable
TOTALNANOPORE=`grep -c "^@" $LONGREADS`
TOTALILLU=`zcat $FWD | grep -c "^@"`
TOTALILLUMINA=`echo "$TOTALILLU*2" | bc`

# repeat for every genome
for name in $MAGS

do

# get the fraction of nanopore reads
NANO=`grep -c $name $LOCAL1/nanopore-unfilteredf-to-final.sam`
NANFRAC=`echo "scale=4; $NANO/$TOTALNANOPORE*100" | bc`

# get the fraction of illumina reads
ILLU=`grep -c $name $LOCAL1/illumina-to-final.sam`
ILLFRAC=`echo "scale=4; $ILLU/$TOTALILLUMINA*100" | bc`

# print fractions of illumina and nanopore reads to table 2
echo "| $name | $NANFRAC | $ILLFRAC |" >> $TABLE3

done

# progress update
echo Finished Table 3


# for name in $MAGS
# do

# done
