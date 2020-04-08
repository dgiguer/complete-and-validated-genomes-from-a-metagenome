#!/usr/bin/python3.6
import sys, re
# Alec Bahcheli
#define the list of mags you are interested in summarizing
name = sys.argv[1]
table = open(sys.argv[2], "a")
table.write("| %s | " %name)

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
        table.write(str(float(ilcov[1].strip("\n")[1:]) / length * 1000) + " | ")
    if "Nanopore" in line:
        nancov = line.split("\t")
        table.write(str(float(nancov[1].strip("\n")[1:]) / length * 1000) + " | ")
table.write(rRNA + " | ")
#get number of tRNAs
trna = open(sys.argv[8]).readlines()
table.write(trna[0].strip("\n") + " | ")
table.write("\n")
print("Finished %s" %name)
table.close()
