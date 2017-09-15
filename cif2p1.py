#!/usr/bin/python2.7

# Daniele ongari 11Sep2017. Using just ASE v3.14 and manage crystal

# This script is using ASE: pip install ase

# USAGE:
# cif2p1_EZ.py input.cif  
# output: input_P1.cif

# PARTICULAR USAGE:
# for f in *; do cif2p1_EZ.py $f ; done

import ase
import ase.io
import ase.build
import sys

inpname = sys.argv[1]
outname = sys.argv[1].split(".")[-2]+"_P1.cif"

# import in ASE and print in pdb (it should be able to recognise and remove overlaps)
cif_original  = ase.io.read(inpname, format='cif') 
cif_unwrapped = ase.build.make_supercell(cif_original, [[1,0,0], [0,1,0], [0,0,1]])
ase.io.write(outname,cif_unwrapped,format='cif')

