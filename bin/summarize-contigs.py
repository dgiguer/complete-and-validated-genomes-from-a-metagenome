#!/usr/bin/python3.6
import sys, re
# Alec Bahcheli
#define the list of mags you are interested in summarizing
name = sys.argv[1]
table = open(sys.argv[2], "a")
table.write("| %s | " %name)

flye = open(sys.argv[5])
unicycler = open(sys.argv[6])
miniasm = open(sys.argv[7])
reassembly = "re-assemblies"

taxonomies = open(sys.argv[3]).readlines()
#get predicted taxonomy
for line in taxonomies:
    if "bin_name" not in line:
        taxon = line.split("\t")
        #split the taxonomy prediction from anvio into a list, print the relevant taxonomic classification (if it has one, otherwise print unknown)
        if len(taxon) > 2 and len(taxon[len(taxon) - 1]) > 2:
            table.write(taxon[len(taxon) - 1].strip("\n") + " | ")
        else:
            newt = []
            for element in taxon:
                if len(element) > 4:
                        newt.append(element)
            if len(newt) > 2:
                table.write(newt[len(newt) - 1].strip("\n") + " | ")
            elif len(newt) <= 2:
                table.write("Unknown | ")
info = open(sys.argv[4]).readlines()
#get total genome length and number of rRNAs
rRNA = ""
length = float(0)
for line in info:
    if "Total Length" in line:
        tl = line.split("\t")
        table.write(tl[len(tl) - 1].strip("\n") + " | ")
        length = int(tl[len(tl) - 1].strip("\n"))
    if "Ribosomal_RNAs" in line:
        tl = line.split("\t")
        rRNA = tl[len(tl) - 1].strip("\n")
contig = open(sys.argv[5]).readlines()
#get gc content
gc = float(0)
for line in contig:
    if ">" not in line:
        line = line.strip("\n")
        for element in line:
            if element == "G" or element == "C":
                gc = gc + 1
table.write(str(round(gc / length * 100, 1)) + " | ")
#get completeness and redundancy
details = open(sys.argv[6]).readlines()
for line in details:
    if name in line:
        dets = line.split("\t")
        table.write(dets[3] + " | " + dets[4] + " | ")
#get nanopore and illumina coverage
coverage = open(sys.argv[7]).readlines()
for line in coverage:
    if "Illumina" in line:
        ilcov = line.split("\t")
        table.write(str(int(ilcov[1]) / length) + " | ")
    if "Nanopore" in line:
        nancov = line.split("\t")
        table.write(str(int(nancov[1]) / length) + " | ")
table.write(rRNA + " | ")
#get number of tRNAs
trna = open(sys.argv[8]).readlines()
table.write(trna[0] + " | ")
#check flye reassembly
def printconditions(x, y):
    if x:
        table.write("1 circular contig | ")
    else:
        ycount = 0
        for line in y:
            if ">" in line:
                ycount = ycount + 1
        if ycount == 1:
            table.write("1 linear contig | ")
        elif ycount > 1:
            table.write("%s linear contigs | " %(str(ycount)))
fly = open(sys.argv[9]).readlines()
flyeinfo = open(sys.argv[10]).readlines()
fl = False
for line in flyeinfo:
    if "#" not in line:
        line = line.split("\t")
        if line[3] == "+":
            for line in flye:
                if name in line:
                    fl = True
printconditions(fl, fly)
#check unicycler assembly
uni = open(sys.argv[11]).readlines()
with open(sys.argv[12]) as file1:
    unicyclerinfo = file1.read()
un = False
if re.search("circular", unicyclerinfo):
    for line in unicycler:
        if name in line:
            un = True
printconditions(un, uni)
#check miniasm reassembly
mini = open(sys.argv[13]).readlines()
with open(sys.argv[14]) as file2:
    miniasminfo = file2.read()
mi = False
if re.search("1 circular unitig", miniasminfo):
    mi = True
printconditions(mi,mini)
# check spades assembly
spades = open(sys.argv[15]).readlines()
# spades will never close a genome / will never specify the genome has been closed
sp = False
printconditions(sp,spades)
# insert spot for graph
table.write("~dot-plots/%s/%s-dotplot_full_genome-0-XXXXXXX.png | " %(name, name))
table.write("\n")
print("Finished %s" %name)
table.close()
