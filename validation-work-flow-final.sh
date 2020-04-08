#!/bin/bash
# Dan Giguere, Alec Bahcheli, Ben Joris
# Before beginning, make sure you have the following:
  # all the forward and reverse short reads are in accessible directories
  # all the long (nanopore) reads are in a directory
  # the path to the original metaFlye assembly with all files maintained there

# set a size cutoff for the mags to meet in order to be analyzed: 300kb eliminates the chloroplast and mitochondria genomes
SIZE=300000

# dot plot parameters: setting to 'Yes' will generate smaller dot plots for every 50kb of the genome in addition to the whole genome dot plot
SMALLPLOTS='Yes'

# choose the number of threads to use for programs capable of multi-threading
THREADS=8

# define the short forward and reverse genes from illumina and the long reads from nanopore
FWD='algae_R1_strict.fastq.gz'
REV='algae_R2_strict.fastq.gz'
LONGREADS='algae_hac_long_reads.fastq'

# define the directory containing the metaFlye assembly
RAWFASTA='flye_assembly_meta_hac'

# directory to work in for the project
MASTERDIR=''

# directory for reassemblies, annotations and dot plots
POLISHING=`echo $MASTERDIR/polishing`
ASSEMBLIES=`echo $MASTERDIR/re-assemblies`
ANVIOANNOTATIONS=`echo $MASTERDIR/contig-annotations`
DOTPLOTS=`echo $MASTERDIR/dot-plots`
TABLE1=`echo $MASTERDIR/output/table1.txt`
TABLE2=`echo $MASTERDIR/output/table2.txt`
TABLE3=`echo $MASTERDIR/output/table3.txt`

# path to necessary tools: bowtie2, bowtie2-build, minimap2, samtools, bcftools, rebaler, reapr, pilon, unicycler, flye, tRNAscan-SE-2.0. this assumes all tools are available in your PATH, or you can input the locations for each
BOWTIE='bowtie2'
BUILD='bowtie2-build'
MAP='minimap2'
SAMTOOLS='samtools'
BCFTOOLS='bcftools'
REBALER='/usr/bin/python3.6 rebaler-runner.py'
REAPR='reapr'
PILON='pilon-1.23.jar'
UNICYCLER='unicycler-runner.py'
FLYE='flye'
SPADES='spades.py'
TRNASCAN='tRNAscan-SE'
# path to python
PYTHON='/usr/bin/python3.6'

# make a directory for the project (as defined in the variable MASTERDIR)
mkdir -p $MASTERDIR $MASTERDIR/raw-contigs $POLISHING/pilon-logs $MASTERDIR/output $MASTERDIR/output/circos $MASTERDIR/output/dotplots

###################################

# extract initial circularized genomes

###################################

# extract circular contigs that are over 300kb in size; input requires the metaFlye directory

# requires: extract-circularized-fasta.py
./extract-circularized-fasta.py $RAWFASTA/assembly_info.txt $RAWFASTA/assembly.fasta $MASTERDIR/raw-contigs/mags-to-analyze.txt $MASTERDIR/raw-contigs $SIZE

# define the mags to analyze from the output of the last script
MAGS=`awk '{print $1}' $MASTERDIR/raw-contigs/mags-to-analyze.txt`

###################################

# validation and polishing processes

###################################

# repeat the following validation and polishing steps for each contig
for name in $MAGS
    do

    # create a directory for the contig validation, polishing, reassembly, etc
    mkdir -p $MASTERDIR/$name $MASTERDIR/$name/output $MASTERDIR/$name/output/fcd_errors $MASTERDIR/$name/output/total_coverage $ASSEMBLIES $ASSEMBLIES/$name $ASSEMBLIES/$name/flye $ASSEMBLIES/$name/unicycler $ASSEMBLIES/$name/unicycler-long $ASSEMBLIES/$name/spades $ANVIOANNOTATIONS $ANVIOANNOTATIONS/$name $ANVIOANNOTATIONS/$name/output $DOTPLOTS $DOTPLOTS/$name $MASTERDIR/$name/circos

    # link the raw fasta file to that directory for analysis
    cp $MASTERDIR/raw-contigs/$name $MASTERDIR/$name/$name-contigs.fa

    # the header line needs to be changed to circular=true for rebaler rebaler
    sed '1 s/\(^>.*$\)/\1circular=true/' $MASTERDIR/$name/$name-contigs.fa > $MASTERDIR/$name/$name-contigs-edited.fa

    # concatenate the genomes into a single fasta file for minimap2 mapping
    cat $MASTERDIR/$name/$name-contigs-edited.fa >> $POLISHING/raw-genomes.fa

done

# map the nanopore long reads to the mag
$MAP -t $THREADS -aQLx map-ont --secondary=no --sam-hit-only $POLISHING/raw-genomes.fa $LONGREADS > $POLISHING/raw-genome-nanopore.sam

# extract read name, bit score, alignment score, and cigar string
# only on non-header lines (lines that don't begin with "@")
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $POLISHING/raw-genome-nanopore.sam > $POLISHING/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
./cigar-parse.py $POLISHING/sam.txt $POLISHING/final-read-names.txt

# filter the nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# get all header lines by searching for lines beginning with "@"
awk '$0 ~ "@"{print $0}' $POLISHING/raw-genome-nanopore.sam > $POLISHING/nanopore-filtered-to-raw.sam

# filter 
# while parsing the first file, save the first field in an array incrementing (++) after each line. Once the first file has been parsed (NOT NR==FNR), stop saving the first field in the array (;next) and evaluate whether the arguments in the array are in a line of the second file. If so, save the whole line
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $POLISHING/final-read-names.txt $POLISHING/raw-genome-nanopore.sam >> $POLISHING/nanopore-filtered-to-raw.sam

# progress
echo Finished filtering nanopore reads against raw assemblies

# activate conda virtual environment to run rebaler
source miniconda3/bin/activate

# loop through all MAGs 
for name in $MAGS
do

    # get sam header
    # get nanopore reads for a given mag by grepping against genome name
    grep "^@" $POLISHING/nanopore-filtered-to-raw.sam > $MASTERDIR/$name/nanopore-filtered-to-raw.sam
    grep $name $POLISHING/nanopore-filtered-to-raw.sam | grep -v "^@" >> $MASTERDIR/$name/nanopore-filtered-to-raw.sam

    # intialize (sort) the sam file using samtools
    $SAMTOOLS sort $MASTERDIR/$name/nanopore-filtered-to-raw.sam -o $MASTERDIR/$name/nanopore-filtered-to-raw-sorted.bam -@ $THREADS

    # convert the sam file into fastq for use in polishing
    $SAMTOOLS fastq $MASTERDIR/$name/nanopore-filtered-to-raw-sorted.bam > $MASTERDIR/$name/nanopore-filtered-to-raw.fastq -@ $THREADS

    # use rebaler to correct the assembly using nanopore
    # changed to the edited file
    $REBALER -t $THREADS $MASTERDIR/$name/$name-contigs-edited.fa $MASTERDIR/$name/nanopore-filtered-to-raw.fastq > $MASTERDIR/$name/rebaler-genome.fasta

    # export the fastq and genomes
    cat $MASTERDIR/$name/rebaler-genome.fasta >> $POLISHING/all-rebaler-genomes.fasta
    cat $MASTERDIR/$name/nanopore-filtered-to-raw.fastq >> $POLISHING/nanopore-filtered.fastq

done

echo Finished rebaler corrections 

# remap the nanopore against the rebaler assemblies
$MAP -t $THREADS -aLQx map-ont --sam-hit-only --secondary=no $POLISHING/all-rebaler-genomes.fasta $POLISHING/nanopore-filtered.fastq > $POLISHING/nanopore-to-rebaler.sam

# make a bowtie2 database to map the illumina reads against the rebaler assemblies
$BUILD $POLISHING/all-rebaler-genomes.fasta $POLISHING/all-rebaler-polished

# map the illumina reads against the rebaler assemblies
# output: $POLISHING/all-illumina-to-rebaler.sam
$BOWTIE -x $POLISHING/all-rebaler-polished -1 $FWD -2 $REV -p $THREADS -S $POLISHING/all-illumina-to-rebaler.sam --no-unal --no-discordant --end-to-end --very-sensitive

# progress update
echo Mapped to rebaler

# sort the sam file
$SAMTOOLS sort -@ $THREADS $POLISHING/nanopore-to-rebaler.sam -o $POLISHING/nanopore-to-rebaler-sorted.bam

# index the bam file 
$SAMTOOLS index -@ $THREADS $POLISHING/nanopore-to-rebaler-sorted.bam

# remove non-primary alignments and convert to bam
$SAMTOOLS view -@ $THREADS -F 256 -bS $POLISHING/all-illumina-to-rebaler.sam > $POLISHING/all-illumina-to-rebaler.bam

# use samtools to sort the mapped reads
$SAMTOOLS sort $POLISHING/all-illumina-to-rebaler.bam -o $POLISHING/all-illumina-to-rebaler-sorted.bam -@ $THREADS

# index the sorted mapped reads creating and index
$SAMTOOLS index -b $POLISHING/all-illumina-to-rebaler-sorted.bam -@ $THREADS

# pilon does not requires subsetting the reads to those specific for a particular contig
# correct the rebaler assembly using both illumina and nanopore reads via pilon. Send logs to log directory. OUTPUTS MULTIPLE FILES
/usr/bin/java -Xmx32G -jar $PILON --genome $POLISHING/all-rebaler-genomes.fasta --frags $POLISHING/all-illumina-to-rebaler-sorted.bam --nanopore $POLISHING/nanopore-to-rebaler-sorted.bam --fix snps,indels --threads $THREADS --verbose --outdir $POLISHING/pilon --changes &> $POLISHING/pilon-logs/pilon.log

# removes "circular-true" from fasta headers
sed 's/circular=true_pilon//g' $POLISHING/pilon/pilon.fasta > $POLISHING/all-pilon-genomes.fasta

# map the filtered nanopore reads to the pilon reassembly
$MAP -t $THREADS -aLQx map-ont --secondary=no $POLISHING/all-pilon-genomes.fasta $POLISHING/nanopore-filtered.fastq > $POLISHING/nanopore-to-pilon.sam

# build a bowtie2 database to map the illumina reads against the pilon polished assemblies
$BUILD $POLISHING/all-pilon-genomes.fasta $POLISHING/all-pilon-genomes

# re-map the short reads against the pilon polished references
$BOWTIE -x $POLISHING/all-pilon-genomes -1 $FWD -2 $REV -p $THREADS -S $POLISHING/illumina-to-all-pilon-genomes.sam --no-unal --no-discordant --end-to-end --very-sensitive

# remove unmapped reads and convert to bam
$SAMTOOLS view -F 4 -bS $POLISHING/nanopore-to-pilon.sam -o $POLISHING/nanopore-to-pilon.bam -@ $THREADS

# use samtools to sort nanpore mapping to the pilon reassembly
$SAMTOOLS sort $POLISHING/nanopore-to-pilon.bam -o $POLISHING/nanopore-to-pilon-sorted.bam -@ $THREADS

# convert the illumina reads to bam and sort them
$SAMTOOLS sort $POLISHING/illumina-to-all-pilon-genomes.sam > $POLISHING/illumina-to-all-pilon-genomes-sorted.bam -@ $THREADS

# progress update
echo Finished mapping to pilon

# index the reads in preparation for reapr
$SAMTOOLS index $POLISHING/illumina-to-all-pilon-genomes-sorted.bam

# use reapr pipeline to check for errors in assembly according to illumina reads
$REAPR pipeline $POLISHING/all-pilon-genomes.fasta $POLISHING/illumina-to-all-pilon-genomes-sorted.bam $POLISHING/reapr

# progress update
echo Finished reapr

# generate mpileup files for variant calling
$BCFTOOLS mpileup --threads $THREADS -o $POLISHING/mpileup.txt -f $POLISHING/all-pilon-genomes.fasta $POLISHING/illumina-to-all-pilon-genomes-sorted.bam

# collect variants from mpileup
$BCFTOOLS call -cv -o $POLISHING/variants.txt --threads $THREADS $POLISHING/mpileup.txt

# filter for variants that have > 50% of reads supporting the call (indels)
$BCFTOOLS filter -i '%MAX(IMF)>0.5' $POLISHING/variants.txt -Ob > $POLISHING/variants-filtered.bcf

# normalization is needed to left align each of the variants
$BCFTOOLS norm -d all -f $POLISHING/all-pilon-genomes.fasta $POLISHING/variants-filtered.bcf -Ob > $POLISHING/variants-norm.bcf

# index these variants
$BCFTOOLS index $POLISHING/variants-norm.bcf -o $POLISHING/variants-norm.bcf.csi

# apply all the variants to a final polished assembly
$BCFTOOLS consensus -f $POLISHING/all-pilon-genomes.fasta $POLISHING/variants-norm.bcf > $POLISHING/final-assemblies.fasta

# extract the final reads from each polsihed assembly
./extract-circularized-fasta.py $POLISHING/final-assemblies.fasta $MASTERDIR

echo Finished all validation and polishing


###################################

# orient genome assemblies on dnaA gene

###################################

for name in $MAGS
    do

    # these steps are done to generate improved nanopore mapping to the final assembly

    # calculating GC content and GC skew

    # from the fasta file itself, use seqinr to calculate the GC content in 1000 base windows that can be added to the files in R

    # activate conda virtual environment to run prokka
    source miniconda3/bin/activate

    # run prokka on all genomes, make sure you have the most recent version of tbl2asn, or this will not work
    prokka --outdir $MASTERDIR/$name/circos/prokka --prefix $name $MASTERDIR/$name/$name-final-assembly.fasta --cpu $THREADS

    # re-orient the genome around the dnaA gene where there is a switch in GC-skew
    grep 'dnaA\|DnaA\|dnaa' $MASTERDIR/$name/circos/prokka/$name.gff > $MASTERDIR/$name/circos/dnaA-loci.txt

    ./orient.py $MASTERDIR/$name/circos/dnaA-loci.txt $MASTERDIR/$name/$name-final-assembly.fasta $MASTERDIR/$name/$name-final-assembly-oriented.fasta $name $MASTERDIR/$name/$name-orientation-shift.txt

    # get the gc content, gc skew and culmulative gc skew in 1000 base windows
    ./circlize_gc_information.R -i $MASTERDIR/$name/$name-final-assembly-oriented.fasta -o $MASTERDIR/$name/circos/$name-gc-info.txt

    # rebuild the final assembly to have 80kb of the beginning of the genome copied to the end of the genome to improve remapping of long reads
    ./modify-genome.py $MASTERDIR/$name/$name-final-assembly-oriented.fasta $MASTERDIR/$name/$name-final-assembly-reverse-start.fasta

    # group all genomes oriented and their reverse pair
    cat $MASTERDIR/$name/$name-final-assembly-reverse-start.fasta >> $POLISHING/all-reverse-start.fasta
    # group all oriented genomes
    cat $MASTERDIR/$name/$name-final-assembly-oriented.fasta >> $POLISHING/all-oriented.fasta

done

# progress update
echo Finished genome orientation

###################################

# Remap the reads to the oriented and reverse-oriented genomes
# Filter according to the previous filtering parameters
# Concatenate the filtered oriented and reverse-oriented mapped reads and convert to fastq for use in reassembly

###################################

# remap unfiltered nanopore reads against the final forward and reverse assemblies
$MAP -t $THREADS -aLQx map-ont --secondary=no --sam-hit-only $POLISHING/all-reverse-start.fasta $LONGREADS > $POLISHING/nanopore-unfilteredr-to-final.sam

# remap unfiltered nanopore reads against the final forward and reverse assemblies
$MAP -t $THREADS -aLQx map-ont --secondary=no --sam-hit-only $POLISHING/all-oriented.fasta $LONGREADS > $POLISHING/nanopore-unfilteredf-to-final.sam

# extract relevant information on nanopore reverse sequence alignments scores
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $POLISHING/nanopore-unfilteredr-to-final.sam > $POLISHING/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
./cigar-parse.py $POLISHING/sam.txt $POLISHING/final-read-names.txt

# filter the reverse nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# this gets the sam headers in the file
awk '$0 ~ "@"{print $0}' $POLISHING/nanopore-unfilteredr-to-final.sam > $POLISHING/nanopore-filteredr-to-final.sam
# while parsing the first file, save the first field in an array incrementing (++) after each line. Once the first file has been parsed (NOT NR==FNR), stop saving the first field in the array (;next) and evaluate whether the arguments in the array are in a line of the second file. If so, save the whole line
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $POLISHING/final-read-names.txt $POLISHING/nanopore-unfilteredr-to-final.sam >> $POLISHING/nanopore-filteredr-to-final.sam

# extract relevant information on nanopore forward sequence alignments scores
awk '$0 !~ "@" {a = $1; b = $2; c = substr($14, 6); d = $6; print a"\t"b"\t"c"\t"d}' $POLISHING/nanopore-unfilteredf-to-final.sam > $POLISHING/sam.txt

# get the length of the sequence and length mapped against the query from the cigar string
./cigar-parse.py $POLISHING/sam.txt $POLISHING/final-read-names.txt

# filter the reverse nanopore mapped reads for a primary alignment (i.e., bit flag in sam output is 0 or 16)
# while parsing the first file, save the first field in an array incrementing (++) after each line. Once the first file has been parsed (NOT NR==FNR), stop saving the first field in the array (;next) and evaluate whether the arguments in the array are in a line of the second file. If so, save the whole line
awk '$0 ~ "@"{print $0}' $POLISHING/nanopore-unfilteredf-to-final.sam > $POLISHING/nanopore-filteredf-to-final.sam
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $POLISHING/final-read-names.txt $POLISHING/nanopore-unfilteredf-to-final.sam >> $POLISHING/nanopore-filteredf-to-final.sam

# make a bowtie2 database for mapping illumina reads to the finally assembly
$BUILD $POLISHING/all-oriented.fasta $POLISHING/genomes-final-assembly

# map the illumina reads against the final assembly
$BOWTIE -x $POLISHING/genomes-final-assembly -1 $FWD -2 $REV -p $THREADS -S $POLISHING/illumina-to-final.sam --no-unal --no-discordant --end-to-end --very-sensitive

# sort the unfiltered-forward nanopore reads
$SAMTOOLS sort $POLISHING/nanopore-unfilteredf-to-final.sam -o $POLISHING/nanopore-unfilteredf-to-final-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/nanopore-unfilteredf-to-final-sorted.bam -@ $THREADS

# sort the unfiltered-reverse nanopore reads
$SAMTOOLS sort $POLISHING/nanopore-unfilteredr-to-final.sam -o $POLISHING/nanopore-unfilteredr-to-final-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/nanopore-unfilteredr-to-final-sorted.bam -@ $THREADS

# sort the filtered-forward nanopore reads
$SAMTOOLS sort $POLISHING/nanopore-filteredf-to-final.sam -o $POLISHING/nanopore-filteredf-to-final-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/nanopore-filteredf-to-final-sorted.bam -@ $THREADS

# sort the filtered-reverse nanopore reads
$SAMTOOLS sort $POLISHING/nanopore-filteredr-to-final.sam -o $POLISHING/nanopore-filteredr-to-final-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/nanopore-filteredr-to-final-sorted.bam -@ $THREADS

# create a fasta index for proper header generation after extracting specific regional reads
$SAMTOOLS faidx $POLISHING/all-reverse-start.fasta

# create a fasta index for proper header generation after extracting specific regional reads
$SAMTOOLS faidx $POLISHING/all-oriented.fasta

# extract the "soft" mapped reads
$SAMTOOLS view -@ $THREADS -f 2 -F 256 -bS $POLISHING/illumina-to-final.sam > $POLISHING/illumina-to-final-soft.bam

# sort the extracted soft reads
$SAMTOOLS sort $POLISHING/illumina-to-final-soft.bam -o $POLISHING/illumina-to-final-soft-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/illumina-to-final-soft-sorted.bam -@ $THREADS

# use more stringent filtering parameters from the illumina alignment against the final assembly in preparation for depth determining
$SAMTOOLS view -@ $THREADS -bS -f 2 -F 3848 $POLISHING/illumina-to-final.sam > $POLISHING/illumina-to-final-hard.bam

# sort the extracted hard reads
$SAMTOOLS sort $POLISHING/illumina-to-final-hard.bam -o $POLISHING/illumina-to-final-hard-sorted.bam -@ $THREADS

# index the sorted reads
$SAMTOOLS index $POLISHING/illumina-to-final-hard-sorted.bam -@ $THREADS

# create a fasta index for proper header generation after extracting specific regional reads
$SAMTOOLS faidx $POLISHING/all-oriented.fasta

for name in $MAGS
do

    # collect the relevant illumina reads
    $SAMTOOLS view $POLISHING/illumina-to-final-soft-sorted.bam $name -o $MASTERDIR/$name/illumina-to-final-soft.bam -T $POLISHING/all-oriented.fasta -@ $THREADS

    $SAMTOOLS view $POLISHING/nanopore-filteredf-to-final-sorted.bam $name -o $MASTERDIR/$name/nanopore-filtered-to-final-forward.bam -T $POLISHING/all-oriented.fasta -@ $THREADS

    $SAMTOOLS view $POLISHING/nanopore-filteredr-to-final-sorted.bam $name-reversed -o $MASTERDIR/$name/nanopore-filtered-to-final-reverse.bam -T $POLISHING/all-reverse-start.fasta -@ $THREADS

    # sort the reads by position
    $SAMTOOLS sort -@ $THREADS $MASTERDIR/$name/nanopore-filtered-to-final-forward.bam -o $MASTERDIR/$name/nanopore-filtered-to-final-forward-sorted.bam

    $SAMTOOLS sort -@ $THREADS $MASTERDIR/$name/nanopore-filtered-to-final-reverse.bam -o $MASTERDIR/$name/nanopore-filtered-to-final-reverse-sorted.bam

    # sort the reads by name
    $SAMTOOLS sort -n -@ $THREADS $MASTERDIR/$name/illumina-to-final-soft.bam -o $MASTERDIR/$name/illumina-to-final-soft-sorted-by-name.bam

    # convert the concatenated filtered bam file to fastq
    $SAMTOOLS fastq $MASTERDIR/$name/nanopore-filtered-to-final-forward-sorted.bam > $MASTERDIR/$name/nanopore-filtered-to-final-forward.fastq -@ $THREADS

    $SAMTOOLS fastq $MASTERDIR/$name/nanopore-filtered-to-final-reverse-sorted.bam > $MASTERDIR/$name/nanopore-filtered-to-final-reverse.fastq -@ $THREADS

    # concatenate the forward and reverse nanopore fastq files
    cat $MASTERDIR/$name/nanopore-filtered-to-final-reverse.fastq > $MASTERDIR/$name/nanopore-filtered-to-final.fastq

    # get the read names mapping
    grep "^@" $MASTERDIR/$name/nanopore-filtered-to-final-reverse.fastq > $MASTERDIR/$name/unique-read-names.txt

    # concatenate reads that are not found in the fastq already to the fastq
    # while parsing the first file, save the first field in an array incrementing (++) after each line. Once the first file has been parsed (NOT NR==FNR), stop saving the first field in the array (;next) and evaluate whether the arguments in the array are in a line of the second file. If so, save the whole line
    awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $MASTERDIR/$name/unique-read-names.txt $MASTERDIR/$name/nanopore-filtered-to-final-forward.fastq | grep -v -A 3 -f - >> $MASTERDIR/$name/nanopore-filtered-to-final.fastq

    # separate the reads into forward and reverse types (could have used bam2fastq, still required a sorted bam file)
    $SAMTOOLS fastq -1 $MASTERDIR/$name/illumina-soft-paired1.fq -2 $MASTERDIR/$name/illumina-soft-paired2.fq -0 /dev/null -s /dev/null -n $MASTERDIR/$name/illumina-to-final-soft-sorted-by-name.bam -@ $THREADS

done 


echo Finished final mapping

###################################

# begin reassembly attempts

###################################

## Re-assembly and analysis using different tools

for name in $MAGS
do

    source miniconda3/bin/activate

    # reassemble using flye not in meta mode with the subset of nanopore reads filtered from the mapping against the final assembly
    flye --nano-raw $MASTERDIR/$name/nanopore-filtered-to-final.fastq --threads $THREADS -g 5m -o $ASSEMBLIES/$name/flye/

    # move the contigs fasta file and assembly info file to up two directories
    mv $ASSEMBLIES/$name/flye/assembly.fasta $ASSEMBLIES/$name-flye-assembly.fasta
    # input: $ASSEMBLIES/$name/flye/assembly_info.txt
    # output: $ASSEMBLIES/$name-flye.log
    mv $ASSEMBLIES/$name/flye/assembly_info.txt $ASSEMBLIES/$name-flye.log

    conda deactivate

    # assemble with unicycler hybrid assembly
    /usr/bin/python3.6 $UNICYCLER -t $THREADS -1 $MASTERDIR/$name/illumina-soft-paired1.fq -2 $MASTERDIR/$name/illumina-soft-paired2.fq -l $MASTERDIR/$name/nanopore-filtered-to-final.fastq -o $ASSEMBLIES/$name/unicycler/ --no_pilon

    # move the relevant files up two directories
    mv $ASSEMBLIES/$name/unicycler/assembly.fasta $ASSEMBLIES/$name-unicycler-assembly.fasta

    mv $ASSEMBLIES/$name/unicycler/unicycler.log $ASSEMBLIES/$name-unicycler.log

    # reassemble using unicycler long-read assembly 
    /usr/bin/python3.6 $UNICYCLER -t $THREADS -l $MASTERDIR/$name/nanopore-filtered-to-final.fastq -o $ASSEMBLIES/$name/unicycler-long --no_pilon

    # move the relevant files up two directories
    mv $ASSEMBLIES/$name/unicycler-long/assembly.fasta $ASSEMBLIES/$name-unicycler-long-assembly.fasta

    mv $ASSEMBLIES/$name/unicycler-long/unicycler.log $ASSEMBLIES/$name-unicycler-long.log

    # assemble with spades
    $SPADES -1 $MASTERDIR/$name/illumina-soft-paired1.fq -2 $MASTERDIR/$name/illumina-soft-paired2.fq --careful -t $THREADS -o $ASSEMBLIES/$name/spades 

    # move the relevant information up two directories
    mv $ASSEMBLIES/$name/spades/contigs.fasta $ASSEMBLIES/$name-spades-assembly.fasta 

    mv $ASSEMBLIES/$name/spades/spades.log $ASSEMBLIES/$name-spades.log

done

echo Finished all re-assemblies 

###################################

# create dot plots for the flye re-assemblies against the validated and corrected assemblies to ensure they are assembled correctly

###################################

for name in $MAGS
do

    # create dot plot(s) for the given genome
    ./dot-plot.py $name $DOTPLOTS $SMALLPLOTS $ASSEMBLIES/$name-flye-assembly.fasta $MASTERDIR/$name/$name-final-assembly.fasta

    mv $DOTPLOTS/$name $MASTERDIR/output/dotplots

done

echo Dot plots made

###################################

# generate  files for circos plots

###################################

source miniconda3/bin/activate

# calculate the depth over a 1000 base window using mapped unfiltered and filtered for reverse reads for nanopore
mosdepth -b 1000 $POLISHING/nanopore-filteredf1 $POLISHING/nanopore-filteredf-to-final-sorted.bam
mosdepth -b 1000 $POLISHING/nanopore-unfilteredf1 $POLISHING/nanopore-unfilteredr-to-final-sorted.bam
mosdepth -b 1000 $POLISHING/nanopore-filteredr1 $POLISHING/nanopore-filteredr-to-final-sorted.bam
mosdepth -b 1000 $POLISHING/nanopore-unfilteredr1 $POLISHING/nanopore-unfilteredr-to-final-sorted.bam

# calculate the depth over a 1000 base window using mapped illumina reads
mosdepth -b 1000 $POLISHING/illumina-soft1 $POLISHING/illumina-to-final-soft-sorted.bam
mosdepth -b 1000 $POLISHING/illumina-hard1 $POLISHING/illumina-to-final-hard-sorted.bam

# it outputs in a gzipped bed file, so unzip
gunzip $POLISHING/nanopore-filteredf1.regions.bed.gz
gunzip $POLISHING/nanopore-unfilteredf1.regions.bed.gz
gunzip $POLISHING/nanopore-filteredr1.regions.bed.gz
gunzip $POLISHING/nanopore-unfilteredr1.regions.bed.gz
gunzip $POLISHING/illumina-soft1.regions.bed.gz
gunzip $POLISHING/illumina-hard1.regions.bed.gz


# apply to each genome
for name in $MAGS
do

    # collect only the relevant genome
    grep $name $POLISHING/nanopore-filteredf1.regions.bed > $MASTERDIR/$name/circos/$name-nanopore-filteredf.regions.bed
    grep $name $POLISHING/nanopore-unfilteredf1.regions.bed > $MASTERDIR/$name/circos/$name-nanopore-unfilteredf.regions.bed
    grep $name $POLISHING/nanopore-filteredr1.regions.bed > $MASTERDIR/$name/circos/$name-nanopore-filteredr.regions.bed
    grep $name $POLISHING/nanopore-unfilteredr1.regions.bed > $MASTERDIR/$name/circos/$name-nanopore-unfilteredr.regions.bed
    grep $name $POLISHING/illumina-soft1.regions.bed > $MASTERDIR/$name/circos/$name-illumina-soft.regions.bed
    grep $name $POLISHING/illumina-hard1.regions.bed > $MASTERDIR/$name/circos/$name-illumina-hard.regions.bed

    # modify the unfiltered.regions.bed and filtered.regions.bed files to include accurate measurements
    ./modify-nanopore-bed.py $MASTERDIR/$name/circos/$name-nanopore-filtered.regions.bed $MASTERDIR/$name/circos/$name-nanopore-unfiltered.regions.bed $MASTERDIR/$name/circos/$name-nanopore-filteredf.regions.bed $MASTERDIR/$name/circos/$name-nanopore-unfilteredf.regions.bed $MASTERDIR/$name/circos/$name-nanopore-filteredr.regions.bed $MASTERDIR/$name/circos/$name-nanopore-unfilteredr.regions.bed $MASTERDIR/$name/$name-final-assembly-oriented.fasta $name

    # generate length of genome
    awk '/^>/{if (l!="") l=0; next}{l+=length($0)}END{print l}' $MASTERDIR/$name/$name-final-assembly-oriented.fasta > $MASTERDIR/$name/circos/$name-final-length.txt

    # get length of genome
    length=`cat $MASTERDIR/$name/circos/$name-final-length.txt`
    # generate cytoband file needed for circlize
    echo -e "$name\t0\t$length\tnot_real\tnot_real" > $MASTERDIR/$name/circos/$name-cytoband.txt

    # generate the required bed files from the .gff output.

    # first two lines are headers
    # i'm printing it into bed format, then printing "1" as a value file.

    # generate header for BED files for positive strand, negative strand, tRNA and rRNA
    echo -e "chr\tstart\tend\tcds_pos" > $MASTERDIR/$name/circos/$name-cds-positive.bed
    echo -e "chr\tstart\tend\tcds_neg" > $MASTERDIR/$name/circos/$name-cds-negative.bed
    echo -e "chr\tstart\tend\ttrna" > $MASTERDIR/$name/circos/$name-cds-trna.bed
    echo -e "chr\tstart\tend\trrna" > $MASTERDIR/$name/circos/$name-cds-rrna.bed

    # use tab separation (OFS arg)
    # coding sequences
    awk -v OFS='\t' '$3 ~ "CDS" && $7 ~ "+" {print $1, $4, $5, "1"}' $MASTERDIR/$name/circos/prokka/$name.gff >> $MASTERDIR/$name/circos/$name-cds-positive.bed

    awk -v OFS='\t' '$3 ~ "CDS" && $7 ~ "-" {print $1, $4, $5, "1"}' $MASTERDIR/$name/circos/prokka/$name.gff >> $MASTERDIR/$name/circos/$name-cds-negative.bed

    # get tRNA and rRNA gene loci
    awk -v OFS='\t' '$3 ~ "tRNA" {print $1, $4, $5, "1"}' $MASTERDIR/$name/circos/prokka/$name.gff >> $MASTERDIR/$name/circos/$name-cds-trna.bed

    awk -v OFS='\t' '$3 ~ "rRNA" {print $1, $4, $5, "1"}' $MASTERDIR/$name/circos/prokka/$name.gff >> $MASTERDIR/$name/circos/$name-cds-rrna.bed

    # modify the bed files based on the positions of dnaa gene rarrangement and size of the genome
    ./bed-file-orientation.py $MASTERDIR/$name/$name-orientation-shift.txt $MASTERDIR/$name/circos/$name-cds-trna.bed $MASTERDIR/$name/circos/$name-cds-rrna.bed $MASTERDIR/$name/circos/$name-cds-positive.bed $MASTERDIR/$name/circos/$name-cds-negative.bed

    mv $MASTERDIR/$name/circos/$name-cds* $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-gc* $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-cytoband.txt $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-illumina-hard.regions.bed $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-illumina-soft.regions.bed $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-nanopore-filtered.regions.bed $MASTERDIR/output/circos
    mv $MASTERDIR/$name/circos/$name-nanopore-unfiltered.regions.bed $MASTERDIR/output/circos

done


###################################

# rough annotations and genome investigations

###################################

for name in $MAGS
do

    # create a file for illumina and nanopore total bp coverages
    echo "$name average genome coverages" > $ANVIOANNOTATIONS/$name/output/$name-coverage.txt

    # get illumina and nanopore total bp coverages
    awk '{sum+=$4} END { print "Illumina coverage\t",sum}' $MASTERDIR/output/circos/$name-illumina-hard.regions.bed >> $ANVIOANNOTATIONS/$name/output/$name-coverage.txt
    awk '{sum+=$4} END { print "Nanopore coverage\t",sum}' $MASTERDIR/output/circos/$name-nanopore-filtered.regions.bed >> $ANVIOANNOTATIONS/$name/output/$name-coverage.txt

    # determine the number of tRNA in the genome
    awk -v OFS='\t' '$3 ~ "tRNA" {print $1, $4, $5, "1"}' $MASTERDIR/$name/circos/prokka/$name.gff | wc -l > $ANVIOANNOTATIONS/$name/output/tRNA-number.txt

    # start anvio6 virtual environment
    source anvio-6/bin/activate

    # copy the genome for annotating
    cp $MASTERDIR/$name/$name-final-assembly-oriented.fasta $ANVIOANNOTATIONS/$name/$name-final-assembly-oriented.fasta

    # generate a contigs database of the final genome assembly using anvio
    anvi-gen-contigs-database -f $ANVIOANNOTATIONS/$name/$name-final-assembly-oriented.fasta -n $name-validation-annotation -o $ANVIOANNOTATIONS/$name/$name-annotation.db

    ##THE FOLLOWING FUNCTIONS MUST BE RUN IN ORDER
    # run hmms on the database (no longer automatic in anvio-6)
    anvi-run-hmms -T $THREADS -c $ANVIOANNOTATIONS/$name/$name-annotation.db

    # affiliate single copy gene taxonomies 
    anvi-run-scg-taxonomy -T $THREADS -c $ANVIOANNOTATIONS/$name/$name-annotation.db

    # estimate genome taxonomy 
    anvi-estimate-genome-taxonomy -T $THREADS --just-do-it -c $ANVIOANNOTATIONS/$name/$name-annotation.db -o $ANVIOANNOTATIONS/$name/$name-taxonomies.txt

    # estimate genomes completeness
    anvi-estimate-genome-completeness -c $ANVIOANNOTATIONS/$name/$name-annotation.db -o $ANVIOANNOTATIONS/$name/$name-completeness.txt

    # get contig stats for the contig
    anvi-display-contigs-stats $ANVIOANNOTATIONS/$name/$name-annotation.db --report-as-text -o $ANVIOANNOTATIONS/$name/contigs-stats.txt

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
    ./table1.py $name $TABLE1 $ANVIOANNOTATIONS/$name/$name-taxonomies.txt $ANVIOANNOTATIONS/$name/contigs-stats.txt $MASTERDIR/$name/$name-final-assembly-oriented.fasta $ANVIOANNOTATIONS/$name/$name-completeness.txt $ANVIOANNOTATIONS/$name/output/$name-coverage.txt $ANVIOANNOTATIONS/$name/output/tRNA-number.txt $ASSEMBLIES/$name-flye-assembly.fasta $ASSEMBLIES/$name-flye.log $ASSEMBLIES/$name-unicycler-assembly.fasta $ASSEMBLIES/$name-unicycler.log $ASSEMBLIES/$name-unicycler-long-assembly.fasta $ASSEMBLIES/$name-unicycler-long.log $ASSEMBLIES/$name-spades-assembly.fasta 


    # use python script to populate TABLE2
    ./table2.py $name $TABLE2 $ANVIOANNOTATIONS/$name/$name-taxonomies.txt $ANVIOANNOTATIONS/$name/contigs-stats.txt $MASTERDIR/$name/$name-final-assembly-oriented.fasta $ANVIOANNOTATIONS/$name/$name-completeness.txt $ANVIOANNOTATIONS/$name/output/$name-coverage.txt $ANVIOANNOTATIONS/$name/output/tRNA-number.txt $ASSEMBLIES/$name-flye-assembly.fasta $ASSEMBLIES/$name-flye.log $ASSEMBLIES/$name-unicycler-assembly.fasta $ASSEMBLIES/$name-unicycler.log $ASSEMBLIES/$name-unicycler-long-assembly.fasta $ASSEMBLIES/$name-unicycler-long.log $ASSEMBLIES/$name-spades-assembly.fasta 

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
    NANO=`grep -c $name $POLISHING/nanopore-unfilteredf-to-final.sam`
    NANFRAC=`echo "scale=4; $NANO/$TOTALNANOPORE*100" | bc`

    # get the fraction of illumina reads
    ILLU=`grep -c $name $POLISHING/illumina-to-final.sam`
    ILLFRAC=`echo "scale=4; $ILLU/$TOTALILLUMINA*100" | bc`

    # print fractions of illumina and nanopore reads to table 2
    echo "| $name | $NANFRAC | $ILLFRAC |" >> $TABLE3

done

# progress update
echo Finished Table 3
