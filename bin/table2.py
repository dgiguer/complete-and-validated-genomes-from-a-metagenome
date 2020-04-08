#!/usr/bin/python3.6
import sys, re
# Alec Bahcheli
#define the list of mags you are interested in summarizing
name = sys.argv[1]
table = open(sys.argv[2], "a")
table.write("| %s | " %name)

reassembly = "re-assemblies"

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
            fl = True
printconditions(fl, fly)
#check unicycler assembly
uni = open(sys.argv[11]).readlines()
with open(sys.argv[12]) as file1:
    unicyclerinfo = file1.read()
un = False
if re.search("circular", unicyclerinfo):
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
table.write("\n")
print("Finished %s" %name)
table.close()
