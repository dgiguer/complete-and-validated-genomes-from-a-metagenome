# Visualizing coverage in circos plots

Daniel J Giguere

# Summary 

This file contains all the code to reproduce the circos style visualizations from the data output for the validation script. The required input files are generated in `validation-work-flow-final.sh` under the Circos step. An example data set is provided. 

### Input files

- table of GC content, skew and culmulative skew calculated from `circlize_gc_information.R`
- unfiltered Illumina coverage
- filtered Illumina coverage
- unfiltered nanopore coverage
- filtered nanopore coverage
- coding sequences (positive and negative strand in separate files)
- location of tRNA and rRNA genes
- cytoband file (required for length and genome name)

Below is the R code required to reproduce the image for the *Blastomonas* genome (Supplemental Figure 4). The `genome` variable assigns the relative path to the file. The output file will be a PDF of the genome (it may take 30 seconds or so to plot all the data). Additional information about the `circlize` package can be found [here](https://github.com/jokergoo/circlize). 

```
library(circlize)

# relative path and genome is loaded into this variable 
genome <- "data/contig_125"

# read in the input files
gc_info <- read.table(paste0(genome, "-gc-info.txt"), header=TRUE)
ill_cov <- read.table(paste0(genome, "-illumina-soft.regions.bed"))
ill_cov_filt <- read.table(paste0(genome, "-illumina-hard.regions.bed"))
nano_cov <- read.table(paste0(genome, "-nanopore-unfiltered.regions.bed"))
nano_cov_filt <- read.table(paste0(genome, "-nanopore-filtered.regions.bed"))

# read in coding sequence bed files
cds_pos <- read.table(paste0(genome, "-cds-positive.bed"), header=TRUE)
cds_neg <- read.table(paste0(genome, "-cds-negative.bed"), header=TRUE)
trna <- read.table(paste0(genome, "-cds-trna.bed"), header=TRUE)
rrna <- read.table(paste0(genome, "-cds-rrna.bed"), header=TRUE)
cytoband.df <- read.table(paste0(genome, "-cytoband.txt"), colClasses=c("character", "numeric", "numeric", "character", "character"))

# generate master bed dataframe for coverage plotting
final_bed <- ill_cov
colnames(final_bed) <- c("chr", "start", "end", "illumina_coverage")

# add coverage information to data frame
final_bed$illumina_coverage_filtered <- ill_cov_filt$V4
final_bed$nanopore_coverage <- nano_cov$V4
final_bed$nanopore_coverage_filtered <- nano_cov_filt$V4

# add gc info columns
final_bed$gc_content <- gc_info$gc_content
final_bed$gc_skew <- gc_info$gc_skew
final_bed$gc_culm <- gc_info$gc_culm

# this is what the required bed file input looks like.
# > head(final_bed)
#         chr start  end illumina_coverage illumina_coverage_filtered
#1 contig_125     0 1000            138.72                     138.72
#2 contig_125  1000 2000            128.83                     128.83
#3 contig_125  2000 3000            185.93                     185.93
#4 contig_125  3000 4000            204.40                     204.40
#5 contig_125  4000 5000            136.34                     136.34
#6 contig_125  5000 6000            150.13                     150.13
#  nanopore_coverage nanopore_coverage_filtered gc_content     gc_skew
#1             89.35                      73.12      0.648  0.28657866
#2             92.51                      74.79      0.660  0.29699970
#3             83.96                      68.93      0.674  0.01530687
#4             84.00                      66.35      0.641 -0.16899646
#5             85.32                      67.56      0.669 -0.33155737
#6             89.91                      69.73      0.662 -0.51428315
#        gc_culm
#1  0.0026374353
#2  0.0053707774
#3  0.0055116493
#4  0.0039563442
#5  0.0009049617
#6 -0.0038280794

```

Since duplicated genes exist in the genome, read coverage can be incorrect reported in the thousands due to mis mapping, which would negatively affect the visualization. This isn't a huge issue for this specific genome, but some of the unfiltered nanopore reads and Illumina reads do get modified to reflect this. None of the filtered nanopore coverage is modified since the max coverage value (99) is lower than 1.75 X the mean (117). I arbitrarily chose 1.75 X as the cut off to visualize. 

```
# calculate maximum nanopore coverage to display.
maximum_nanopore <- mean(final_bed$nanopore_coverage_filtered) * 1.75
maximum_illumina <- mean(final_bed$illumina_coverage) * 1.75

# replace each one that is over limit each one with maximum
final_bed$nanopore_coverage[which(final_bed$nanopore_coverage > maximum_nanopore)] <- maximum_nanopore

# adjust maximum for purpose of plotting.
final_bed$illumina_coverage[which(final_bed$illumina_coverage > maximum_illumina)] <- maximum_illumina
final_bed$illumina_coverage_filtered[which(final_bed$illumina_coverage_filtered > maximum_illumina)] <- maximum_illumina

# clear plotting region if circos plot already active.
circos.clear()

# open PDF device. this will output the file in the data folder - modify as you wish. 
pdf(paste0(genome, "-circos.pdf"), width = 12, height=15)

# rotate so the beginning of the genome is at the top
# make space between start and end for labels
# remove as much space as possible between layers.
circos.par(start.degree=75, gap.after = 30, cell.padding=c(0,0,0,0))

# this is what the cytoband.df input needs to look like. V4 and V5 are required only when you are plotting multiple chromosomes (like human data for example)
#> cytoband.df
#          V1 V2      V3       V4       V5
#1 contig_125  0 3836404 not_real not_real

# initialize using the cytoband dataframe but don't plot the default.
circos.initializeWithIdeogram(cytoband.df, plotType=NULL)

# generate track for genome.
circos.track(ylim = c(0,1), panel.fun = function(x,y) {
}, track.height = 0.05, bg.col="#5e6cff")

# plot all coverage tracks
circos.genomicTrackPlotRegion(final_bed, numeric.column=c("illumina_coverage", "illumina_coverage_filtered", "nanopore_coverage", "nanopore_coverage_filtered") , ylim = c(-(maximum_nanopore), maximum_illumina), panel.fun = function(region, value, ...) {

  circos.genomicLines(region, value, numeric.column="illumina_coverage", baseline = 0, col = rgb(0,0,0,0.5), border = NA, area=TRUE, ...)

  circos.genomicLines(region, value, numeric.column="illumina_coverage_filtered", baseline = 0, col = "#9cc7ff", border = NA, area=TRUE, ...)

  circos.genomicLines(region, -value, numeric.column="nanopore_coverage", baseline = 0, col = rgb(0,0,0,0.6), border = NA, area=TRUE, ...)

  circos.genomicLines(region, -value, numeric.column="nanopore_coverage_filtered", baseline = 0, col = "#f5a556", border = NA, area=TRUE, ...)

  }, track.height = 0.2, bg.border=NA)

# add maximum coverage values for nanopore and illumina reads for scale.
circos.text(0, maximum_illumina, adj=c(1.2, 0), labels = as.integer(maximum_illumina), cex = 0.7)
circos.text(0, 0, labels = "0", cex = 0.7, adj=c(1.2, 0))
circos.text(0, -maximum_nanopore, labels = as.integer(maximum_nanopore), adj=c(1.2, 0), cex = 0.7)
 

# plot gc_content track
circos.genomicTrackPlotRegion(final_bed, numeric.column="gc_content" , ylim = c(min(final_bed$gc_content),max(final_bed$gc_content)), panel.fun = function(region, value, ...) {

    circos.genomicLines(region, value, numeric.column="gc_content", baseline = 0, col = rgb(0,0,0,0.7), border = NA, ...)

    },track.height = 0.05)

# gc skew and culm track
# the gc_skew and gc_culm columns need to converted to a proportion (i.e. between 0.0 to 1.0) to be plotted on the same track.

# convert to proportion for plotting purposes. 
final_bed$gc_skew <- final_bed$gc_skew / max(final_bed$gc_skew)
final_bed$gc_culm <- final_bed$gc_culm / max(abs(final_bed$gc_culm))

# plot GC skew and culmulative skew on the same track.
circos.genomicTrackPlotRegion(final_bed, numeric.column="gc_culm", ylim = c(-1, 1), panel.fun = function(region, value, ...) {

  circos.genomicLines(region, value, numeric.column="gc_culm", baseline = 0, col = rgb(0,0,0,0.5), border = NA, area=TRUE,...)

  circos.genomicLines(region, value, numeric.column="gc_skew", baseline = 0, col = rgb(0,0,0,1), border = NA,...)

}, bg.border = NA, track.height = 0.05)

# plot positive coding sequence track
circos.genomicTrackPlotRegion(cds_pos, numeric.column="cds_pos", ylim = c(0, 1), panel.fun = function(region, value, ...) {

  # in every region, print a vertical line at the start region.
  for (k in seq_len(nrow(region))){
    # plot vertical lines for each START of the ORF from bottom to top.
      circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 0.5, straight = TRUE, col = rgb(0,0,0,0.5))
  }

}, bg.border = NA, track.height = 0.05)

# plot negative coding sequence track
circos.genomicTrackPlotRegion(cds_neg, numeric.column="cds_neg", ylim = c(0, 1), panel.fun = function(region, value, ...) {

  # in every region, print a vertical line at the start region.
  for (k in seq_len(nrow(region))){
    # plot vertical lines for each START of the orf from bottom to top.
      circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 0.5, straight = TRUE, col = rgb(0,0,0,0.5))
  }

}, bg.border = NA, track.height = 0.05)

# plot tRNA and rRNA track
circos.genomicTrackPlotRegion(list(trna, rrna), ylim = c(0, 1), panel.fun = function(region, value, ...) {

  # this will iterate through the list of bed dataframes (trna and rrna)
  i = getI(...)

  # plot tRNAs
  if (i == 1){

    for (k in seq_len(nrow(region))){
      # plot vertical lines for each START of the orf from bottom to top.
        circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 2, straight = TRUE, col = rgb(0,0,0,0.2))
    }

  } else {
    # plot rRNAs as red
    for (k in seq_len(nrow(region))){
      # plot vertical lines for each START of the orf from bottom to top.
        circos.lines(rep(region[k, 1], 2), c(0, 1), lwd = 2, straight = TRUE, col = rgb(1,0,0,1))
    }
  }
}, bg.border = NA, track.height = 0.05)

# clear the plotting area
circos.clear()

# plot the labels in the middle
text(0, 0.85, "Illumina coverage", cex = 1)
text(0, 0.75, "Nanopore coverage", cex = 1)
text(0, 0.66, "GC content", cex = 1)
text(0, 0.6, "GC skew", cex = 1)
text(0, 0.52, "CDS (+)", cex = 1)
text(0, 0.45, "CDS (-)", cex = 1)
text(0, 0.38, "tRNA (Black)\n rRNA (red)", cex = 1)
text(0, 0, "Blastomonas", cex = 2)
text(0, -0.1, paste0("Genome size: ", cytoband.df$V3), cex = 1)

dev.off()
```

The output from this script can be visualized [here](data/contig_125-circos.pdf).