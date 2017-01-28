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


openbabel="/usr/local/bin/obabel"
manage_crystal="/home/daniele/Programs/python_my/manage_crystal.py"


# improve the quality of the cif file
subprocess.call(openbabel+" "+sys.argv[1]+" -O cif2pytmp1.cif >> cif2pytmp2.verbose", shell=True)

# change the cif to forces ASE to read "_symmetry_equiv_pos_as_xyz" instead of other stuff: obabel should print the list of operations to do correctly
subprocess.call("sed -i '/_space_group_symop_operation_xyz/c\_symmetry_equiv_pos_as_xyz'         cif2pytmp1.cif", shell=True)
subprocess.call("sed -i '/_symmetry_space_group_name_H-M/c\_symmetry_space_group_name_H-M  P 1'  cif2pytmp1.cif", shell=True)
subprocess.call("sed -i '/_space_group_name_H-M_alt/c\_space_group_name_H-M_alt  P 1'            cif2pytmp1.cif", shell=True)

# import in ASE and print in pdb (it should be able to recognise and remove overlaps)
asecrystal=ase.io.read("cif2pytmp1.cif",0) #Reads both cif and pdb
ase.io.write("cif2pytmp3.pdb",asecrystal,"pdb")

# use manage_crystal to convert from pdb to cif in the way I want
outname = sys.argv[1].split(".")[-2]+"_"+sys.argv[2]+".cif"
subprocess.call(manage_crystal+" cif2pytmp3.pdb "+outname+" >> cif2pytmp4.verbose", shell=True)

subprocess.call("rm -f cif2pytmp*", shell=True)


