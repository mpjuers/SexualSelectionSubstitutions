#!/usr/bin/env python3
# Example usage: overlapDetector.py -i matingBehavior.interest.txt -w 15000
import sys

# This deals with the command line input from the user
for x in range(1, len(sys.argv)):
    if sys.argv[x] == "-i":
        input_file = str(sys.argv[x + 1])
        view = open(input_file)
        sys.stderr.write(
            "Your input file is: " + str(sys.argv[x + 1])
        )  # Prints this to Standard error rather than Standard Out
    elif sys.argv[x] == "-w":
        windowlen = int(sys.argv[x + 1])
        sys.stderr.write(
            "You want to group positions within "
            + str(sys.argv[x + 1])
            + " bp of each other; "
        )

"""
#This is the hardcoded filename and settings used for while writing the code
view = open("matingBehavior.interest.txt")
#This is the hardcoded window length for position overlap
windowlen = 15000
"""

# Read the different lines of the file
lines = view.readlines()
pos = []
for line in lines:
    line = line.strip()  # removes the /n at the end of the line
    line = line.partition("_")  # uses _ as the seperator in the line
    pos.append(
        [line[0], line[2]]
    )  # Appends each position and chromosome to the object pos

# Create a dictionary of all the chromosomes as keys, with all positions a values
chrom_dict = dict()
for position in pos:
    if position[0] in chrom_dict:
        # append the new number to the existing array at this slot
        chrom_dict[position[0]].append(int(position[1]))
    else:
        # create a new array in this slot
        chrom_dict[position[0]] = [int(position[1])]
# https://stackoverflow.com/questions/3199171/append-multiple-values-for-one-key-in-a-dictionary

posGroupList = []  # Creates an empty list for all our sections of interest
for chrom in chrom_dict:
    newlist = chrom_dict[
        chrom
    ]  # converts the list of position values into an actual list
    newlist.sort()  # sort the list of positions in ascending order
    i = 0
    chromList = []
    # print(newlist)
    for position in newlist:
        if len(chromList) == 0:
            chromList.append([chrom, position])
            i = i + 1
        else:
            if (position - newlist[i - 1]) <= windowlen:
                chromList[i - 1].append(position)
            else:
                chromList.append([chrom, position])
                i = i + 1
    # print(chromList) #for testing purposes, to check if the within-chromosme lists are correct
    posGroupList.extend(chromList)

# Here is the code for printing the grouped positions
# print(posGroupList)
print("chromosome,start,end")  # Prints a header column, "chromosome,start,end"
for position in posGroupList:  # Prints to StdOut in csv format
    if len(position) == 2:
        start = position[1] - int(windowlen / 2)
        end = position[1] + int(windowlen / 2)
        print(str(position[0]) + "," + str(start) + "," + str(end))
    else:
        index = len(position)
        start = position[1] - int(windowlen / 2)
        end = position[-1] + int(windowlen / 2)
        print(str(position[0]) + "," + str(start) + "," + str(end))
