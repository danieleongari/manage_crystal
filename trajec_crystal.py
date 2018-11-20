#!/usr/bin/env python

from __future__ import print_function  #python3 like print()
from __future__ import absolute_import
import string, sys
import numpy
import math
import subprocess
import argparse
from argparse import RawTextHelpFormatter # to go next line in the help text
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
    "inpfile",
    type=str,
    help="path to the input file to read\n" +
    "IMPLEMENTED: pdb")

parser.add_argument(
    "-o",
    "--output",
    type=str,
    default=None,
    help="Output filename.extension or just the extension\n" +
    "IMPLEMENTED: xyz, axsf")

parser.add_argument(
    "-sss",
    "-savesnapshots",
    action="store_true",
    dest="sss",
    default=False,
    help="Don't clear the single snapshots, in the input and output format")

args = parser.parse_args()

################################################## Manage input/output filenames
# Open input file
if not os.path.isfile(args.inpfile):
    sys.exit("ERROR: The file %s doesn't exist!" % args.inpfile)
inpfilename = os.path.splitext(args.inpfile)[0]
inpformat = os.path.splitext(args.inpfile)[1][1:]
inpfile = "{}.{}".format(inpfilename,inpformat)
inpfiler = open(inpfile, 'r')

# Parse output filename
if args.output == None:
    sys.exit("ERROR: An -o output is mandatory!")
if len(args.output.split(".")) > 1:  #output defined as name.format
    outfilename = os.path.splitext(args.output)[0]
    outformat = os.path.splitext(args.output)[1][1:]
else:  #output defined as format
    outfilename = inpfilename
    outformat = args.output
outfile = "{}.{}".format(outfilename,outformat)

################################ Split and convert snapshots from the trajectory
if inpformat == 'pdb':
    ss = 0 #snapshot count
    line = None
    # Split snapshots
    while True:
        line = inpfiler.readline()
        data = line.split()
        if line == "":
            print("Extracted %d snapshots" %ss)
            break
        elif len(data)>0 and data[0]=="MODEL":
            ssfile = "{0}-ss{1:05d}.pdb".format(outfilename,ss)
            ssfilew = open(ssfile, 'w')
        elif len(data)>0 and data[0]=="ENDMDL":
            print("END", end="", file=ssfilew)
            ssfilew.close()
            ss+=1
        else:
            print(line.strip(), end="\n", file=ssfilew)
    inpfiler.close()
    # Convert snapshots
    print("Converting %s to %s " %(inpformat,outformat), end="")
    for iss in range(ss):
        print(".", end="")
        sys.stdout.flush()
        cmd = ["manage_crystal.py",
               "{0}-ss{1:05d}.pdb".format(outfilename,iss),
               "-o",
               "{0}-ss{1:05d}.{2}".format(outfilename,iss,outformat),
               "-silent",
              ]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output, error = process.communicate()
    print(" all done!")

############################################ Finally combine the temporary files
if outformat=="xyz":
    with open(outfile, 'w') as outfilew:
        for iss in range(ss):
            ssfile = "{0}-ss{1:05d}.{2}".format(outfilename,iss,outformat)
            with open(ssfile) as ssfiler:
                outfilew.write(ssfiler.read())

if outformat=="axsf":
    with open(outfile, 'w') as outfile:
        print("ANIMSTEPS %s" %ss, file=outfile)
        print("CRYSTAL %s" %ss, file=outfile)
        for iss in range(ss):
            print("PRIMVEC %s" %(iss+1), file=outfile)
            fname = "{0}-ss{1:05d}.{2}".format(outfilename,iss,outformat)
            with open(fname) as inpfiler:
                for i, line in enumerate(inpfiler): #skip header
                    if i>2 and i<6:
                        print(line.strip(), file=outfile)
            print("PRIMCOORD %s" %(iss+1), file=outfile)
            with open(fname) as inpfiler:
                for i, line in enumerate(inpfiler): #skip header
                    if i>6:
                        print(line.strip(), file=outfile)

############################################################### Delete snapshots
if not args.sss:
    print("Deleting snaphsots (use -sss to save them)")
    for iss in range(ss):
        os.remove("{0}-ss{1:05d}.{2}".format(outfilename,iss,inpformat))
        os.remove("{0}-ss{1:05d}.{2}".format(outfilename,iss,outformat))
else:
    print("Saved snaphsots")
