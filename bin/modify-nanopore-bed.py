#!/usr/bin/python3.6
# Alec Bahcheli
# modify the bed files to correctly include read depth

import sys, math

results1 = sys.argv[1]
results2 = sys.argv[2]
forfiltered = sys.argv[3]
forunfiltered = sys.argv[4]
revfiltered = sys.argv[5]
revunfiltered = sys.argv[6]
orientedgenome = sys.argv[7]
name = sys.argv[8]
print(name)
genome = open(orientedgenome)
filterresults = open(results1, "w")
unfilterresults = open(results2, "w")
bedff = open(forfiltered).readlines()
beduf = open(forunfiltered).readlines()
bedfr = open(revfiltered).readlines()
bedur = open(revunfiltered).readlines()
seq = ""
# get genome length
for line in genome:
    if ">" not in line:
        seq = seq + line.strip("\n")
# get genome length to nearest (greater) thousand
genomelength = math.ceil(len(seq) / 1000)
# determine where to cut genome (half of genome to lowest thousand)
cut = math.floor((math.floor(len(seq) / 1000)) / 2)
# define approximately a quarter of the genome
fragment = math.floor(cut / 2)
startforward = cut - fragment
stopforward = cut + fragment
startopposite = stopforward - cut
end = genomelength - cut
stopopposite = startforward + end
# create a list for the new bed file information
filtered = []
x = 0
# correct the length of the genome position

# get coverage from the beginning to the first quarter of the genome (using the reverse / opposite oriented genome)
for line in bedfr[end:stopopposite]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str((1 + x) * 1000)
    filtered.append('\t'.join(original))
    x += 1
# get coverage from the first quarter to the third quarter of the genome (using the forward oriented genome)
for line in bedff[startforward:stopforward]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str((1 + x) * 1000)
    filtered.append('\t'.join(original))
    x += 1
# get coverage from the third quarter to the end of the genome (using the reverse / opposite oriented genome)
for line in bedfr[startopposite:end]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str(min(((1 + x) * 1000), len(seq)))
    if ((1 + x) * 1000) > len(seq):
        original[3] = str(round(float(original[3]) * ((1 + x) / (len(seq) / 1000)), 2))
    filtered.append('\t'.join(original))
    x += 1
filterresults.write(''.join(filtered))
filterresults.close()
# repeat for unfiltered genome
unfiltered = []
x = 0
# get coverage from the beginning to the first quarter of the genome (using the reverse / opposite oriented genome)
for line in bedur[end:stopopposite]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str((1 + x) * 1000)
    unfiltered.append('\t'.join(original))
    x += 1
# get coverage from the first quarter to the third quarter of the genome (using the forward oriented genome)
for line in beduf[startforward:stopforward]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str((1 + x) * 1000)
    unfiltered.append('\t'.join(original))
    x += 1
# get coverage from the third quarter to the end of the genome (using the reverse / opposite oriented genome)
for line in bedur[startopposite:end]:
    original = line.split("\t")
    original[0] = name
    original[1] = str(x * 1000)
    original[2] = str(min(((1 + x) * 1000), len(seq)))
    if ((1 + x) * 1000) > len(seq):
        original[3] = str(round(float(original[3]) * ((1 + x) / (len(seq) / 1000)), 2))
    unfiltered.append('\t'.join(original))
    x += 1
unfilterresults.write(''.join(unfiltered))
unfilterresults.close()
genome.close()
