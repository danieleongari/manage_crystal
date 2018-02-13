#!/usr/bin/python3
"""/
Compute radius from cube density (B.Wells, eq 18). Tested for Gaussian.
"""

import string,sys
import numpy
import math
import subprocess
import matplotlib.pyplot as plt

############################################################################# HELP 
if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='-help' or sys.argv[1]=='help':
        print()
        print('####################################################################################')
        print('#  Python program to read a cube file and print infos:')
        print('#')
        print('#  $ %s inputfile.cube' % (sys.argv[0]))
        print('####################################################################################')
        print()
        sys.exit()

############################################################################# utilities
ANGS2BOHR=1.88973 
BOHR2ANGS=1./1.88973 

############################################################################## INPUT


npoints=[0]*3
binsize=[0]*3 #works only with orthogonal cells 

#reading input file: name and format
for nfile in range(1,len(sys.argv)):
	print("Reading: %s" %sys.argv[nfile])
	file = open(sys.argv[nfile],'r')
	header1 = file.readline()
	header2 = file.readline()

	data = file.readline().split()
	nat = int(data[0])
	orig = [ float(data[1]), float(data[2]), float(data[3]) ] #(bohr)

	data = file.readline().split() # X
	npoints[0]=float(data[0])
	binsize[0]=float(data[1])      #(bohr)                       
	data = file.readline().split() # Y
	npoints[1]=float(data[0])
	binsize[1]=float(data[2])      #(bohr)
	data = file.readline().split() # Z
	npoints[2]=float(data[0])
	binsize[2]=float(data[3])      #(bohr)

	if nat>1:
	  print("This script works only with one atom: here you have more!")
	  sys.exit()
	else:
	#for i in range(0,nat): 
	  data = file.readline().split()
	  atpos = [ float(data[2]), float(data[3]), float(data[4]) ] #(bohr) REMEBER the format of the line: AN AN X Y Z

	gridValues=[]
	count=-1
	x=0
	y=0
	z=0
	rxrho=0
	rho=0

	while True:
	 data = file.readline().split()
	 if len(data)==0:           #if the file is finished stop  
	    break
	 for gridvalue in data:
	   count=count+1;
	   z=orig[2]+(count % npoints[2])*binsize[2]
	   y=orig[1]+math.floor(count/npoints[2]) % npoints[1] * binsize[1]
	   x=orig[0]+math.floor(count/npoints[1]/npoints[2]) % npoints[0] * binsize[0]
 
	   r=math.sqrt((x-atpos[0])**2+(y-atpos[1])**2+(z-atpos[2])**2)
	   #print(x,y,z,r)
	   rxrho=rxrho+r*float(gridvalue)
	   rho=rho+float(gridvalue)

	print("radius (bohr)= %f" %(rxrho/rho))
