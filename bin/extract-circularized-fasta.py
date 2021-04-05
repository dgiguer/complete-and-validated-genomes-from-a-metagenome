#!/usr/bin/python3.6
# Alec Bahcheli
# intended to extract sequences of interest


import sys, re

if len(sys.argv) == 6:
    circ_contigs=[]
    contoanalyze = open(sys.argv[3], "w")
    infofh = open(sys.argv[1])
    for line in infofh:
        if "#" not in line:
            info = line.split("\t")
            if (info[3] == "Y" or info[3] == "+") and int(info[1]) > int(sys.argv[5]):
                contoanalyze.write(info[0] + " \n")
                circ_contigs.append(info[0])
    infofh.close()
    contoanalyze.close()
    # Write fasta entry to new file if circularized
    writeflag=False
    with open(sys.argv[2]) as allfh:
        for line in allfh:
            if ">" in line:
                writeflag=False
            if line[1:].strip() in circ_contigs:
                writeflag=True
                filename = line[1:].strip()
            if writeflag:
                temp = sys.argv[4]
                with open("%s/%s" %(temp, filename), "a") as circfh:
                    circfh.write(line)

elif len(sys.argv) == 7:
    circ_contigs=[]
    contoanalyze = open(sys.argv[3], "w")
    infofh = open(sys.argv[1])
    for line in infofh:
        if "#" not in line:
            info = line.split("\t")
            if info[3] == "+" and int(info[1]) <= int(sys.argv[5]):
                contoanalyze.write(info[0] + " \n")
                circ_contigs.append(info[0])
    infofh.close()
    contoanalyze.close()
    # Write fasta entry to new file if circularized
    writeflag=False
    with open(sys.argv[2]) as allfh:
        for line in allfh:
            if ">" in line:
                writeflag=False
            if line[1:].strip() in circ_contigs:
                writeflag=True
                filename = line[1:].strip()
            if writeflag:
                temp = sys.argv[4]
                with open("%s/%s" %(temp, filename), "a") as circfh:
                    circfh.write(line)

# check if the input is a certain length to use this function
elif len(sys.argv) == 4:
    # get the name of the fasta sequence of interest
    circ_contigs = sys.argv[3]
    # check if there are multiple names as input
    if len(circ_contigs.split("\n")) > 1:
        # create a list if there are multiple names
        circ_contigs = circ_contigs.split("\n")
    # Write fasta entry to new file if circularized
    writeflag=False
    # open the original fasta file with many sequences
    with open(sys.argv[1]) as allfh:
        # parse file
        for line in allfh:
            # check headers
            if ">" in line:
                # don't write unless name in header of sequence
                writeflag=False
                # modify line so it's readable as a list
                line1 = line.strip(">").split(" ")
                # check if the names of interest are a list or a single string
                if isinstance(circ_contigs, list):
                    # parse names if list
                    for thing in circ_contigs:
                        # parse the list of elements in the fasta line
                        for spot in line1:
                            # if name is of interest
                            if re.search(thing, spot):
                                # allow the script to write the header and following lines into output file
                                writeflag=True
                                # get the name for part of the file name
                                filename = thing
                        # other possibility for checking if name is in sequence header: if the name is in the beginning after the ">"
                        if line[1:].strip() == thing or line1[0] == thing:
                            # write if name of interest
                            writeflag=True
                            filename = thing
                # if name of interest is alone and not a list
                else:
                    # check if the beginning (after the ">") if the name of interest
                    if line[1:].strip() == circ_contigs or line1[0] == circ_contigs:
                        writeflag=True
                        filename = circ_contigs
                    # otherwise check each element of the line in list format
                    else:
                        for spot in line1:
                            if spot.strip(",") == circ_contigs:
                                print("yes")
                                writeflag=True
                                filename = circ_contigs
            # if the name of interest was in the header, write the header and sequence
            if writeflag:
                temp = sys.argv[2]
                with open("%s/%s" %(temp, filename), "a") as circfh:
                    circfh.write(line)

else:
    with open(sys.argv[1]) as allfh:
        writeflag=True
        for line in allfh:
            if ">" in line:
                filename = line[1:].strip()
            if writeflag:
                directory = sys.argv[2]
                with open("%s/%s/%s-final-assembly.fasta" %(directory, filename, filename), "a") as circfh:
                    circfh.write(line)



