#!/usr/bin/env python

# This script is using ASE and opanbabel (see Drive for the installation)
# to apply the symmetries of a cif file making a cif P1 file
# which can be read by manage_crystal.py
# Since the parser of ASE is crappy I need obabel to make a nice cif file
# but I needed to modify:
# /usr/lib/python2.7/dist-packages/ase/io/cif.py (line: 233)
#     elif '_space_group_name_h-m_alt' in tags:
#        symbolHM = tags['_space_group_name_h-m_alt']

# USAGE:
# cif2p1.py input.cif flag 
# output: input_flag.cif

# PARTICULAR USAGE:
# for f in *; do cif2p1.py $f untouch; done

import ase.io
import string,sys,os
import subprocess


#subprocess.call("cp "+sys.argv[1]+" "+sys.argv[1]+"_copy", shell=True)

openbabel="/usr/local/bin/obabel"
subprocess.call(openbabel+" "+sys.argv[1]+" -O cif2pytmp1.cif >> cif2pytmp2.verbose", shell=True)

asecrystal=ase.io.read("cif2pytmp1.cif",0) #Raads both cif and pdb
 
ase.io.write("cif2pytmp3.pdb",asecrystal,"pdb")

manage_crystal="/home/daniele/Programs/python_my/manage_crystal.py"

outname = sys.argv[1].split(".")[-2]+"_"+sys.argv[2]+".cif"
subprocess.call(manage_crystal+" cif2pytmp3.pdb "+outname+" >> cif2pytmp4.verbose", shell=True)

subprocess.call("rm -f cif2pytmp*", shell=True)


