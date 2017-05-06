#!/usr/bin/python2.7
"""/
Python program to print information from a cube file. Daniele Ongari 7/11/16
"""

import string,sys
import numpy
import math
import subprocess
import matplotlib.pyplot as plt

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

shift  =[0]*3
npoints=[0]*3
binsize=[0]*3 #works only with orthogonal cells 


#reading input file: name and format
for nfile in range(1,len(sys.argv)):
        print
        print "Reading: %s" %sys.argv[nfile]
	file = open(sys.argv[nfile],'r')
	header1 = file.readline()
	header2 = file.readline().split()
        #shift[0]=float(header2[22])
	#shift[1]=float(header2[23])
	#shift[2]=float(header2[24])

	data = file.readline().split()
	nat = int(data[0])
	orig = [ float(data[1])/ANGS2BOHR, float(data[2])/ANGS2BOHR, float(data[3])/ANGS2BOHR ]

	data = file.readline().split() # X
        npoints[0]=float(data[0])
        binsize[0]=float(data[1])/ANGS2BOHR                          
	data = file.readline().split() # Y
        npoints[1]=float(data[0])
        binsize[1]=float(data[2])/ANGS2BOHR
	data = file.readline().split() # Z
        npoints[2]=float(data[0])
        binsize[2]=float(data[3])/ANGS2BOHR

	for i in range(0,nat):
	  data = file.readline().split()

	gridValues=[]
        count=-1
        x=0
        y=0
        z=0

	while True:
	 data = file.readline().split()
	 if len(data)==0:           #if the file is finished stop  
	    break
         for gridvalue in data:
           count=count+1;
           z=orig[2]+(count % npoints[2])*binsize[2]
           y=orig[1]+math.floor(count/npoints[2]) % npoints[1] * binsize[1]
           x=orig[0]+math.floor(count/npoints[1]/npoints[2]) % npoints[0] * binsize[0]
 
           if (x>orig[0]+shift[0]) and (y>orig[1]+shift[1]) and (z>orig[2]+shift[2]):
           	gridvalue=float(gridvalue)
           	#print "%f %f %f %f" %(x,y,z,gridvalue)
	   	gridValues.append(gridvalue)		#(kJ)

	gridValues=[float(i) for i in gridValues]

	print "MIN value = %.3f " %min(gridValues)
        print "MAX value = %.3f " %max(gridValues)

        #Compute Henry coefficient
        
        #T= 298      					#K
        #R= 0.008314  					#kJ/mol/K
        #beta=1/R/T;  					#mol/kJ

        #boltz_coeff=[ math.exp(-beta*U) for U in gridValues]
        #boltz_ave=sum(boltz_coeff)/len(gridValues)
        #mu_ex=-1/beta*math.log(boltz_ave)
        #K_H=beta*boltz_ave                          	#(mol/kJ)
        #print
        #print "Thermodynamic data at %d K:" %T
        #print "mu_ex (kJ/mol) = %.3f "                         %mu_ex
        #print "K_H (mol/Pa/kg) = %.3f (mol/J) / rho_framework (kg/m3)"  %(K_H/1000)
	

        #print
        #print "Printing histogram of values"   
        #f, ax = plt.subplots()   
        #ax.hist(gridValues, 1000, normed=1)
        #plt.yscale('log')
        #ax.set_title(sys.argv[nfile])
        #ax.set_xlabel('Values')
        #ax.set_ylabel('Normalized histogram')
	#plt.show()
	
	
        file.close()





