#!/usr/bin/python3.6
# Alec Bahcheli

import sys, re, math
dnaa = sys.argv[1]
finalassembly = sys.argv[2]
output = sys.argv[3]
name = sys.argv[4]
shiftvalue = open(sys.argv[5], "w")

dnaa=open(dnaa).readlines()
genome = open(finalassembly)
# open a file for outputting the oriented genome
orientedgenome = open(output, "w")
# spots will be the potential ori-starts
spot = []
# define the raw genome
raw = ""
# get the dnaa gene loci
for line in dnaa:
    if name in line:
        line = line.split("\t")
        spot.append(int(line[3]) - 1)

if len(spot) != 0:
    for line in genome:
        if ">" in line:
            line = re.sub("circular=true_pilon",'',line)
            orientedgenome.write(line.strip("\n") + "\n")
        else:
            raw = raw + line.strip("\n")
    # there are usually 60nts per line in fasta sequences
    numlines = math.ceil(len(raw) / 60)
    result = raw[int(spot[0]):] + raw[:int(spot[0])]
    for i in range(numlines):
        orientedgenome.write(result[(i * 60): min(((i + 1) * 60), len(raw))] + "\n")
    orientedgenome.close()
else:
    for line in genome:
        if ">" in line:
            line = re.sub("circular=true_pilon",'',line)
            orientedgenome.write(line.strip("\n") + "\n")
        else:
            raw = raw + line.strip("\n")
    # there are usually 60nts per line in fasta sequences
    numlines = math.ceil(len(raw) / 60)
    result = raw[int(spot[0]):] + raw[:int(spot[0])]
    for i in range(numlines):
        orientedgenome.write(result[(i * 60): min(((i + 1) * 60), len(raw))] + "\n")
    orientedgenome.close()
genome.close()
shiftvalue.write(str(spot[0]) + "\t" + str(len(raw)) + "\n")
shiftvalue.close()












# # if there are multiple dnaa genes check gc content and find the greatest change in gcskew near the dnaa gene
# if len(spot) > 1:
#     gcskew = []
#     gc = gc[1:]
#     locations = []
#     locus = {}
#     sp = 100
#     # get the 
#     for line in gc:
#         thing = line.split()
#         gcskew.append(float(thing[1]))
#     # for each potential spot, check the surrounding areas plus or minus 10000kb
#     for element in spot:
#         locations.append(math.floor(element / 1000))
#         locus[math.floor(element / 1000)] = element
#     for element in locations:
#         gccon1 = sum(gcskew[(element - 100):element]) / 101
#         gccon2 = sum(gcskew[element:(element + 100)]) / 101
#         locus[abs(gccon1 - gccon2)] = element
#     for key in locus.keys():
#         sp = min(key, sp)
#     start = locus[locus[sp]]
#     for line in genome:
#         if ">" in line:
#             line = line.strip("circular=true_pilon")
#             orientedgenome.write(line.strip("\n"))
#         else:
#             raw = raw + line.strip("\n")
#     # there are usually 60nts per line in fasta sequences
#     numlines = math.ceil(len(raw) / 60)
#     result = raw[start:] + raw[:start]
#     for i in range(numlines):
#         orientedgenome.write(result[(i * 60): min(((i + 1) * 60), len(raw))] + "\n")
#     orientedgenome.close()
# if there is only 1 dnaa gene, orient at that gene
## Dan recommended orienting on an arbitrary dnaa gene initially then modifying after