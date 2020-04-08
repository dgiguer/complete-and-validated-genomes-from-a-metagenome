#!/usr/bin/python3.6
# Alec Bahcheli

import sys, re

# get the sam file concatenated and open a file to store the filtered reads
samfile = open(sys.argv[1]).readlines()
results = open(sys.argv[2], "w")

# parse each sam entry in the concatenated file
for line in samfile:
    line = line.strip("\n")
    line = line.split("\t")
    length = 0
    query = 0
    cigar = re.findall("([0-9]*[MISH])", line[3])
    for element in cigar:
        if re.search("[ISH]", element):
            length += int(element.strip("[ISH]")) 
        elif re.search("M", element):
            length += int(element.strip("M"))
            query += int(element.strip("M"))
    line[3] = length
    line.append(query)
    if int(line[2]) > 1 and (int(line[1]) == 0 or int(line[1]) == 16) and line[3] > 1000 and (line[4] / line[3]) > 0.9 and (line[3] / int(line[2])) < 0.9:
        results.write(line[0] + "\n")
results.close()

