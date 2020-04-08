#!/usr/bin/env Rscript
# Alec Bahcheli
library(optparse)
library(seqinr)

option_list <- list(
    make_option(c("-q","--query"), type="character", default=NULL,
            help="query fasta file for dotplot [default]",
            dest="query_sequence"),
    make_option(c("-s","--sequence"), type="character", default=NULL,
            help="sequences fasta file for dotplot [default]",
            dest="sequence_file"),
make_option(c("-l","--labels"), type="character", default=NULL,
            help="dotplot labels off [default]",
            dest="labels"),
    make_option(c("-w","--winsize"), type="character", default=NULL,
            help="window size",
            dest="window"),
    make_option(c("-t","--steps"), type="character", default=NULL,
            help="step size",
            dest="steps"),
    make_option(c("-m","--matches"), type="character", default=NULL,
            help="minimum number of matches in window",
            dest="matches"),
    make_option(c("-o","--output"), type="character", default=NULL,
            help="path output file [default = ./]",
            dest="output_directory")
            )
    
parser <- OptionParser(usage = "%prog -q query.fasta -s sequence.fasta -o out_directoy [options]",option_list=option_list)
            
opt = parse_args(parser)

#Read in a fasta as a list of vectors of characters - this means that each header will be 1 item in the list
seq1 <- read.fasta(file=opt$query_sequence,as.string=FALSE,seqtype=c("DNA"))
# seq2 <-read.fasta(file=opt$sequence_file,as.string=FALSE,seqtype=c("DNA"))

# generate a dotplot of the first vector in the list against itself without labels on axis
if (opt$labels == 0){
        png(opt$output_directory)
        seqinr::dotPlot(seq1[[1]],seq1[[1]],wsize=as.numeric(opt$window),wstep=as.numeric(opt$steps),nmatch=as.numeric(opt$matches), xlab='', ylab='', xaxt='n', yaxt='n')
        dev.off()
# generate a dotplot of the first vector in the list against itself
} else {
        pdf(opt$output_directory)
        seqinr::dotPlot(seq1[[1]],seq1[[1]],wsize=as.numeric(opt$window),wstep=as.numeric(opt$steps),nmatch=as.numeric(opt$matches), xlab='', ylab='')
        dev.off()
}



