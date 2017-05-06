#!/usr/bin/python2.7

import os
import re
import sys
import subprocess

manage_crystal="/home/daniele/Programs/python_my/manage_crystal.py"
pb="/home/daniele/Documents/Project6_zeopp_volume/poreblazer_v3.0.2_d/poreblazer.exe" 

if len(sys.argv)==1: sys.exit("Usage: runporeblazer yourcrystal.xxx. Output: ./yourcrystal_PB/")

for i in range(1,len(sys.argv)):
        filename=sys.argv[i]
        name=filename.split(".")[0]

        #Read the cell dimensions from the cif file
	file = open(filename, "r")
	name=filename.split(".")[0]
        print " %d %s" %(i,name)
 	for line in file:
            if re.search("_cell_length_a", line):
               	A=line.split()[1]
            if re.search("_cell_length_b", line):
               	B=line.split()[1]
            if re.search("_cell_length_c", line):
               	C=line.split()[1]
            if re.search("_cell_angle_alpha", line):
               	a=line.split()[1]
            if re.search("_cell_angle_beta", line):
               	b=line.split()[1]
            if re.search("_cell_angle_gamma", line):
               	c=line.split()[1]
        file.close()

        #make the xyz file for pb
	subprocess.call(manage_crystal+" "+filename+" -o xyz -silent", shell=True)
        
        #make the input file with the info for pb 
	pbinp= open(name+".inp","w")
 	print >> pbinp, "%s.xyz"   %name
 	print >> pbinp, "%s %s %s" %(A,B,C)
 	print >> pbinp, "%s %s %s" %(a,b,c)
        pbinp.close()
 
        #run pb
	subprocess.call(pb+" < "+name+".inp >> "+name+".out", shell=True)
	subprocess.call("mkdir "+name+"_PB", shell=True)
	subprocess.call("mv new.xyz "+name+".inp " +name+".out "+name+".xyz "+name+"_PB", shell=True)
	
	
	






