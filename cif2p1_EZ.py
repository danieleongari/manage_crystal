#!/usr/bin/python2.7

# Daniele ongari 15jan2017. EZ is a modification that uses only obabel and ase (not manage_crystal)

# This script is using ASE and opanbabel (see Drive for the installation)
# to apply the symmetries of a cif file making a cif P1 file
# which can be read by manage_crystal.py
# Since the parser of ASE is crappy I need obabel to make a nice cif file
# but I needed to modify:
# /usr/lib/python2.7/dist-packages/ase/io/cif.py (line: 233)
#add     elif '_space_group_name_h-m_alt' in tags:
#add       symbolHM = tags['_space_group_name_h-m_alt']

# USAGE:
# cif2p1_EZ.py input.cif  
# output: input_P1.cif

# PARTICULAR USAGE:
# for f in *; do cif2p1_EZ.py $f ; done

import ase.io
import string,sys,os
import subprocess


openbabel="/usr/local/bin/obabel" #be sure you have a recent version of obabel (_atom_site_fract_x instead of _atom_site_Cartn_x)
#manage_crystal="/home/daniele/Programs/python_my/manage_crystal.py"


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
outname = sys.argv[1].split(".")[-2]+"_P1.cif"
subprocess.call(openbabel+" cif2pytmp3.pdb -O "+outname+" >> cif2pytmp4.verbose", shell=True)

subprocess.call("rm -f cif2pytmp*", shell=True)


