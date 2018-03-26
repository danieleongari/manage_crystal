#!/usr/bin/env python2
"""/
Python program to obtain guestdensity from a raspa Movie pdb. Daniele Ongari 11/12/16

Future idea:
- choose how many models do and the step 
- insert a gaussian occupancy instead of the spherical step
"""

import string,sys
import numpy
import math
import subprocess
import matplotlib.pyplot as plt
import argparse 
from argparse import RawTextHelpFormatter #needed to go next line in the help text

parser = argparse.ArgumentParser(description='Python program to make a CUBE file for the occupancy, given a PDB movie.', formatter_class=RawTextHelpFormatter)

parser.add_argument("PDB", 
                       type=str,
                       help="path to the PDB file with the movie")

parser.add_argument("outputfile",
                             type=str,
                             help="Name of the output CUBE file")

parser.add_argument("-b",
                             type=float,
                          default=0.1,
                             dest="binlengthinput",
                             help="Binsize of the grid (Angstrom)\n" + 
                                  "Default: 0.1")

parser.add_argument("-r",
                             type=float,
                          default=1.0,
                             dest="radius",
                             help="Radius of the occupation (Angstrom)\n" +
                                  "[Use r<b, if you want to count +1 only for the\n"+ 
                                  " closest gridpoint. It doesn't matter the value!]\n"+  
                                  "[It correspond to 3*sigma if 'gauss' method is used]\n"+        #+/- 3*sigma contains the 99.7% of the gaussian                       
                                  "Default: 1.0")

parser.add_argument("-a",
                             type=str,
                          default="readAll",
                             dest="atomName",                             
                             help="Name of the atom to consider\n" +
                                  "Default: readAll")

parser.add_argument("-method",
                             type=str,
                          default="unit",
                             dest="method",                             
                             help="Method used: 'unit' or 'gauss' integration.\n" +
                                  "[unit: all the area around for a radius r will be counted as +1]\n"+
                                  "[gauss gaussian normalized value around the particles]\n"+
                                  "Default: unit")

parser.add_argument("-norm",
                             type=str,
                          default="atom",
                             dest="norm",                             
                             help="Normalization  used: 'model' or 'atom'.\n" +
                                  "Default: 'model'")

args = parser.parse_args()
args.outputfile+=".cube"

############################################################################# variables
ANGS2BOHR=1.88973 

############################################################################## INPUT reading
ABC=[0]*3           #UC dimensions
abc=[0]*3           #UC angles
ABCbins=[0]*3       #number of bins in each direction
binLenReal=[0]*3    #bin length
dist=[0]*3
occ=[]              #occupation, not normalized 
radiusBin=[0]*3
sigma=args.radius/3 #sigma for "gauss" method

print
print "Reading: %s"    %args.PDB
print "will print: %s" %args.outputfile
file = open(args.PDB,'r')

model_count=0
atom_count=0

skip = file.readline()
data = file.readline().split()
ABC[0]=float(data[1])         
ABC[1]=float(data[2])  
ABC[2]=float(data[3])     
abc[0]=math.radians(float(data[4]))     
abc[1]=math.radians(float(data[5]))     
abc[2]=math.radians(float(data[6]))                              #only orthogonal for now
        

ABCbins[0]=int(math.floor(ABC[0]/args.binlengthinput));
ABCbins[1]=int(math.floor(ABC[1]/args.binlengthinput)); 
ABCbins[2]=int(math.floor(ABC[2]/args.binlengthinput));   
   
binLenReal[0]=ABC[0]/float(ABCbins[0]-1)                  #rounding the specified binLength to the UC
binLenReal[1]=ABC[1]/float(ABCbins[1])   
binLenReal[2]=ABC[2]/float(ABCbins[2])   
           
radiusBin[0]=int(math.floor(args.radius/binLenReal[0]))+1 #specified integration radius, in bins
radiusBin[1]=int(math.floor(args.radius/binLenReal[1]))+1 
radiusBin[2]=int(math.floor(args.radius/binLenReal[2]))+1

for x in range(0,ABCbins[0]):
  for y in range(0,ABCbins[1]):
   for z in range(0,ABCbins[2]):
             occ.append(0)

while True:
                data = file.readline().split()
        	if len(data)==0:           #if the file is finished stop  
	            break 	
        	elif data[0]=='MODEL':     
                    model_count+=1      
                    #print "model %d" %model_count
	           
        	elif data[0]=='ATOM' and (data[2]==args.atomName or args.atomName=="readAll" ):  #READ only the atoms specified with [-a] or all tyhe atom if not specified    
                    atom_count+=1     
	            pos=[float(data[4]),float(data[5]),float(data[6])]                                                                            #position of the read atom
                    posGrid=[int(math.floor(pos[0]/binLenReal[0])), int(math.floor(pos[1]/binLenReal[1])), int(math.floor(pos[2]/binLenReal[2]))] #floor bin corresponding to the position

                    #Now I consider the square around the atomic position (next lines are considering the PBC)
                    for x in range(posGrid[0]-radiusBin[0],posGrid[0]+radiusBin[0]):    
                      if   (x<0): 
                           x=x+ABCbins[0]
                      elif (x>=ABCbins[0]): 
                           x=x-ABCbins[0]

                      for y in range(posGrid[1]-radiusBin[1],posGrid[1]+radiusBin[1]):
                        if   (y<0): 
                             y=y+ABCbins[1]
                        elif (y>=ABCbins[1]): 
                             y=y-ABCbins[1]

	                for z in range(posGrid[2]-radiusBin[2],posGrid[2]+radiusBin[2]):
                          if  (z<0): 
                              z=z+ABCbins[2]-1
                          elif (z>=ABCbins[2]): 
                               z=z-ABCbins[2]+1 

                          #Now I compute the distance so that I can take only the sphere (inside of the square) around the atom				
                          dist[0]=math.fabs(x*binLenReal[0]-pos[0])+0.0001                    
                          dist[1]=math.fabs(y*binLenReal[1]-pos[1])+0.0001 
                          dist[2]=math.fabs(z*binLenReal[2]-pos[2])+0.0001 

                          dist[0]=dist[0]-ABC[0]*int(round(dist[0]/ABC[0])) #pbc
                          dist[1]=dist[1]-ABC[1]*int(round(dist[1]/ABC[1])) #pbc
                          dist[2]=dist[2]-ABC[2]*int(round(dist[2]/ABC[2])) #pbc
                          
                          dist3d=math.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
                          
                          if (dist3d<args.radius):
                             #print "%d %d %d" %(x,y,z) 
                             index=x*ABCbins[1]*ABCbins[2]+y*ABCbins[2]+z
                             if (args.method=='unit'):  occ[index]+=1
                             if (args.method=='gauss'): occ[index]+=1/math.sqrt(((2*math.pi)**3)*(sigma**6))*math.exp(-dist3d**2/(2*sigma**2)) #https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function
                             #print "%d" %occ[index]

file.close()

       
if (args.norm=='model'): occ_norm=[ (float(i)/float(model_count)) for i in occ] #nurmalizing over the number of "models" (=snapshot in Raspa)
if (args.norm=='atom'):  occ_norm=[ (float(i)/float(atom_count))  for i in occ] #normalizing over the number of "atoms"

############################################################################## OUTPUT writing

ofile=open(args.outputfile, 'w+')
print >> ofile, "*** Computed from: %s"                      %args.PDB
print >> ofile, "*** using binlenght: %f ,radius: %f ,atomName: %s" %(args.binlengthinput,args.radius, args.atomName)

  
print >> ofile,"%6d %12.6f %12.6f %12.6f" %(1,0,0,0); #natoms, origin x,y,z

print >> ofile,"%6i %12.6f %12.6f %12.6f" %(ABCbins[0],binLenReal[0]*ANGS2BOHR,0,0) #gridpoints x
print >> ofile,"%6i %12.6f %12.6f %12.6f" %(ABCbins[1],0,binLenReal[1]*ANGS2BOHR,0) #gridpoints y
print >> ofile,"%6i %12.6f %12.6f %12.6f" %(ABCbins[2],0,0,binLenReal[2]*ANGS2BOHR) #gridpoints z
   
print >> ofile,"%6i %12.6f %12.6f %12.6f" %(40,0.0,0.0,0.0),  #one fake atom needed in the cube file to be read by VMD

index=0
for x in range(0,ABCbins[0]):
  for y in range(0,ABCbins[1]):
    column=1
    for z in range(0,ABCbins[2]):
              if (column==1):
                 print >> ofile,"\n%.6f" %occ_norm[index],
                 column+=1;
              elif (column==6):
                 print >> ofile," %.6f" %occ_norm[index],
                 column=1;
              else:
                 print >> ofile," %.6f" %occ_norm[index],
                 column+=1;
              index+=1




