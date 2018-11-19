#!/usr/bin/env python

from __future__ import print_function  #python3 like print()
from __future__ import absolute_import
import string, sys
import numpy
import math
import subprocess
import argparse
from argparse import RawTextHelpFormatter  #needed to go next line in the help text
import os
import re  #re.split('(\d+)',"O23") = ['O', '23', '']
from six.moves import range
import subprocess

parser = argparse.ArgumentParser(
    description=
    'Program that split a tragectory into different files to be converted with\
    manage_crystal and then glues them back(by Daniele Ongari)',
    formatter_class=RawTextHelpFormatter)

parser.add_argument(
    "inputfile",
    type=str,
    help="path to the input file to read\n" +
    "IMPLEMENTED: pdb")

parser.add_argument(
    "output",
    type=str,
    help="Output filename.extension or just the extension\n" +
    "IMPLEMENTED: xyz")

parser.add_argument(
    "-sss",
    "-savesnapshots",
    action="store_true",
    dest="sss",
    default=False,
    help="Don't clear the single snapshots, in the input and output format")

args = parser.parse_args()

# Open input file and
if not os.path.isfile(args.inputfile):
    sys.exit("ERROR: The file %s doesn't exist!" % args.inputfile)
inputfilename = os.path.splitext(args.inputfile)[0]
inputformat = os.path.splitext(args.inputfile)[1][1:]
inpfile = open("{}.{}".format(inputfilename,inputformat), 'r')

# Parse output filename
if len(args.output.split(".")) > 1:  #output defined as name.format
    outputfilename = os.path.splitext(args.output)[0]
    outputformat = os.path.splitext(args.output)[1][1:]
else:  #output defined as format
    outputfilename = inputfilename
    outputformat = args.output
outputfile = "{}.{}".format(outputfilename,outputformat)

# Split and convert snapshots from the trajectory
if inputformat == 'pdb':
    ss = 0 #snapshot count
    line = None
    # Split snapshots
    while True:
        line = inpfile.readline()
        data = line.split()
        if line == "":
            print("Extracted %d snapshots" %ss)
            break
        elif len(data)>0 and data[0]=="MODEL":
            outfile = open("{0}-ss{1:05d}.pdb".format(outputfilename,ss), 'w')
        elif len(data)>0 and data[0]=="ENDMDL":
            print("END", end="", file=outfile)
            outfile.close()
            ss+=1
        else:
            print(line.strip(), end="\n", file=outfile)
    inpfile.close()
    # Convert snapshots
    print("Converting %s to %s " %(inputformat,outputformat), end="")
    for iss in range(ss):
        print(".", end="")
        sys.stdout.flush()
        cmd = ["../../manage_crystal.py",
               "{0}-ss{1:05d}.pdb".format(outputfilename,iss),
               "-o",
               "{0}-ss{1:05d}.{2}".format(outputfilename,iss,outputformat),
               "-silent",
              ]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output, error = process.communicate()
    print(" all done!")

# Finally combine the temporary files
if outputformat=="xyz":
    with open(outputfile, 'w') as outfile:
        for iss in range(ss):
            fname = "{0}-ss{1:05d}.{2}".format(outputfilename,iss,outputformat)
            with open(fname) as inpfile:
                outfile.write(inpfile.read())

# Delete snapshots
if not args.sss:
    print("Deleting snaphsots (use -sss to save them)")
    for iss in range(ss):
        os.remove("{0}-ss{1:05d}.{2}".format(outputfilename,iss,inputformat))
        os.remove("{0}-ss{1:05d}.{2}".format(outputfilename,iss,outputformat))
else:
    print("Saved snaphsots")
