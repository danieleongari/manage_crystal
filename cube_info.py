#!/usr/bin/env python
"""/
Python program to print information from a cube file. Daniele Ongari 7/11/16
"""

import string,sys
import numpy
import math
import subprocess

############################################################################# HELP 
if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='-help' or sys.argv[1]=='help':
        print
        print '####################################################################################'
	print '#  Python program to read a cube file and print infos:'
	print '#'
	print '#  $ %s inputfile.cube' % (sys.argv[0])
        print '####################################################################################'
	print
	sys.exit()

############################################################################# utilities
ANGS2BOHR=1.88973 
BOHR2ANGS=1./1.88973 

############################################################################## INPUT

#reading input file: name and format
for nfile in range(1,len(sys.argv)):
        print
        print sys.argv[nfile]
	file = open(sys.argv[nfile],'r')
	header1 = file.readline()
	header2 = file.readline()
	data = file.readline().split()
	nat = int(data[0])
	orig = [ float(data[1])/ANGS2BOHR, float(data[2])/ANGS2BOHR, float(data[3])/ANGS2BOHR ]
	data = file.readline().split() # X
	data = file.readline().split() # Y
	data = file.readline().split() # Z

	for i in range(1,nat):
	  data = file.readline().split()

	gridValues=[]

	while True:
	 data = file.readline().split()
	 if len(data)==0:           #if the file is finished stop  
	    break
	 gridValues.extend(data)

	gridValues=[float(i) for i in gridValues]

	print "MIN: %.3f " %min(gridValues)
        print "MAX: %.3f " %max(gridValues)

        file.close()





