#!/usr/bin/python
# Alec Bahcheli
import sys, re, os

# get input parameters from parent shell script (eg contig name, relative path, start and stop nts for the dot plot, size of fragments to be analyzed, windows size for analysis, step size and miminum number of matches)
name=sys.argv[1]
local3=sys.argv[2]
smallplots=sys.argv[3]
flyeassembly = sys.argv[4]
finalassembly = sys.argv[5]
start = 0
genome=open(flyeassembly)

# get information on the genome to be analyzed; if the start and stop parameters are not provided, assume analysis of the while genome
goi=""
header=""
for line in genome:
    if ">" in line:
        header=line
    else:
        goi = goi + line.strip("\n")
start = 0
stop = len(goi)
size = len(goi)
winsize = 50000
steps = 15000
matches = 15000
# evaluate the number of windows to use based on the genome size 
windows = (stop - start) / size + 1
x = 1
# for each window, generate a dot plot
while x < windows:
    interest=open("%s/concat%s.fasta" %(local3, x), "w")
    interest.write(header)
    interest.write(goi[(start + ((x - 1) * size)):(min((start + (x * size)), len(goi)))] + "\n")
    interest.close()
    s1 = str(start + ((x - 1) * size))
    s2 = str(min(start + (x * size), len(goi)))
    cmd = "./dot-plot.R -q %s/concat%s.fasta -s %s -l 0 -o %s/%s-dotplot_full_genome-%s-%s.png -w %s -t %s -m %s" %(local3, x, finalassembly, local3, name, s1, s2, winsize, steps, matches)
    os.system(cmd)
    cmd = "rm %s/concat%s.fasta" %(local3, x)
    os.system(cmd)
    x = x + 1
genome.close()
# determine if small plots are to be made; if yes, do 50kb segments
if smallplots == "Yes":
    winsize = 500
    # using a matches size ~35% of the window size generates reasonably noisy plots that allows you to distinguish repeated regions
    matches = 175
    steps = 200
    size = 50000
    windows = (stop - start) / size + 1
    x = 1
    while x < windows:
        interest=open("%s/concat%s.fasta" %(local3, x), "w")
        interest.write(header)
        interest.write(goi[(start + ((x - 1) * size)):(min((start + (x * size)), len(goi)))] + "\n")
        interest.close()
        s1 = str(start + ((x - 1) * size))
        s2 = str(min(start + (x * size), len(goi)))
        cmd = "./dot-plot.R -q %s/concat%s.fasta -s %s -l 1 -o %s/%s-dotplot-%s-%s.png -w %s -t %s -m %s" %(local3, x, finalassembly, local3, name, s1, s2, winsize, steps, matches)
        os.system(cmd)
        cmd = "rm %s/concat%s.fasta" %(local3, x)
        os.system(cmd)
        x = x + 1
