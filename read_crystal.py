#!/usr/bin/env python
"""/
Python program to read coordinates from a file and handle them. Daniele Ongari 7/11/16

inputfilename          name of the input filename.inpformat
inputformat	       input format

natoms                 number of atoms
ABC[012]               A B C in angstrom
abc[012]               alpha,beta,gamma in radians
cell(3x3)              unit cell matrix in angstrom


  atom[index]          atomic name
    an[index]          atomic number
   xyz[index][x,y,z]   cartesian coordinates in angstrom
 fract[index][x,y,z]   fractional coordinates
charge[index]          partial charge

"""

import string,sys
import numpy
import math
import subprocess

############################################################################# HELP 
if sys.argv[1]=='-h' or sys.argv[1]=='-help' or sys.argv[1]=='help':
        print
        print '####################################################################################'
	print '#  Python program to read coordinates from a file and handle them:'
	print '#'
	print '#  $ %s inputfile.xxx outputfile.yyy z' % (sys.argv[0])
	print '#'
	print '#  xxx=xyz(CELL),pdb,cssr          (next: cp2k-restart, xsf, pwo, pwi, gaussian, dcd+atoms)'
	print '#  yyy=cif,pdb,cssr,xyz(CELL)'
	print '#  z=f,l (for the first or the last coordinate in a dcd or pwo or axsf or log)'
        print '####################################################################################'
	print
	sys.exit()

############################################################################# STANDARD INFOS about the periodic table
atom_tuple=('unassigned','H','He',
            'Li','Be','B','C','N','O','F','Ne',
            'Na','Mg','Al','Si','P','S','Cl','Ar',
            'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
            'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo')

atom_count=[0]*119 #anything assigned to 0, H_index=1, He_index=2, ...

############################################################################## INPUT

#reading input file: name and format
inputfilename = sys.argv[1].split(".")[0]
inputformat= sys.argv[1].split(".")[1]
file = open('./'+inputfilename+'.'+inputformat,'r')

if inputformat=='dcd':
	from pwtools import dcd
	cc,co = dcd.read_dcd_data('./'+inputfilename+'.'+inputformat)
	ABC=[ cc[-1][0], cc[-1][1], cc[-1][2] ]
	abc=[ math.acos(cc[-1][3]), math.acos(cc[-1][4]), math.acos(cc[-1][5]) ]
	cell=numpy.matrix([[                 ABC[0],                    0.0,                                                                           0.0],
			   [ABC[1]*math.cos(abc[2]),ABC[1]*math.sin(abc[2]),                                                                           0.0],
			   [ABC[2]*math.cos(abc[1]),ABC[2]*math.cos(abc[0]),math.sqrt(ABC[2]**2-(ABC[2]*math.cos(abc[1]))**2-(ABC[2]*math.cos(abc[0]))**2)]]) 
	xyz=co[-1]
	file.close()
	user_input = raw_input("Give me a file with the atoms: ")
	file = open(user_input,'r')
	natoms=int(file.readlines()[0].split()[0])
	atom=[]
	an=[]
	#xyz=[]
	charge=[]
	file = open(user_input,'r')
	for i in range(0,natoms):
		line = file.readlines()[i+2]
		data = line.split()
		atom.append(data[0])	
		#an.append(atom_tuple.index(atom[i]))
                #atom_count[an[i]]+=1
		an.append(1)
                atom_count[an[i]]+=1
		#xyz.append([float(data[1]), float(data[2]), float(data[3])])
		charge.append(0.000)                                               #fix

#if inputformat=='pwo':
#	pwo2scf=subprocess.call("which pwo2xsf", shell=True)
#	if sys.argv[3]=="f":
#		subprocess.call(pwoscf+"--inicoor    "+sys.argv[1]" > "inputfilename".xsf", shell=True)
#	elif sys.argv[3]=="l":
#		subprocess.call(pwoscf+"--latestcoor "+sys.argv[1]" > "inputfilename".xsf", shell=True)
#	else:
#		print "********* use f for first and l for last coordinate*****"
#	file.close()
#	file = open('./'+inputfilename+'.xsf','r')
#	#now the xsf file still need to be processed
	

#if inputformat=='xsf' or inputformat=='pwo':

if inputformat=='pdb':
	while True:
		line = file.readline()
		if line.split()[0]=='CRYST1':
			break
	ABC=[float(line[06:15]),float(line[15:24]),float(line[24:33])]
	abc=[math.radians(float(line[33:40])),math.radians(float(line[40:47])),math.radians(float(line[47:54]))]
	#read atom[index]
        atom=[]
        an=[]
        xyz=[]
        charge=[]
	i=0
        while True:
         line = file.readline()
	 data = line.split()
	 if len(data)==0:           #if the file is finished stop  
		break 			
 	 elif data[0]=='END':        #if the file says "END" 
		break
	 else:		
                data = line.split()
                atom.append(data[2])
                an.append(atom_tuple.index(atom[i]))
                atom_count[an[i]]+=1
		xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                charge.append(0.000)
		i=i+1	
	 natoms=i

if inputformat=='cssr':
	line = file.readline()
 	celldim=line.split( )
	ABC=[float(celldim[0]),float(celldim[1]),float(celldim[2])]
	line = file.readline()
 	celldim=line.split( )
	abc=[math.radians(float(celldim[0])),math.radians(float(celldim[1])),math.radians(float(celldim[2]))]
	line = file.readline()
        natoms=int(line.split( )[0])
        line = file.readline() 
        atom=[]
        an=[]
        fract=[]
        charge=[]
	i=0
	for i in range(0,natoms):
		line = file.readline()
		data = line.split( )
		atom.append(data[1])	
		an.append(atom_tuple.index(atom[i]))
                atom_count[an[i]]+=1
		fract.append([float(data[2]), float(data[3]), float(data[4])])
		charge.append(0.000)

if inputformat=='xyz':
	#read number of atoms
	line = file.readline()
	natoms=int(line.split( )[0])

	#read cell in my way of writing it as a comment of xyz
	line = file.readline()
	celldim=line.split( )
	#debug: celldim=['cell:',1,2,3,4,5,6,7,8,9,10,11]

	if celldim[0]=='CELL:' or celldim[0]=='CELL':
		ABC=[float(celldim[1]),float(celldim[2]),float(celldim[3])]
		abc=[math.radians(float(celldim[4])),math.radians(float(celldim[5])),math.radians(float(celldim[6]))]

	elif celldim[0]=='cell:' or celldim[0]=='cell':
		cell=numpy.matrix([[celldim[1],celldim[2],celldim[3]],
		                   [celldim[4],celldim[5],celldim[6]],
		                   [celldim[7],celldim[8],celldim[9]]])
	#read atom[index]
	atom=[]
	an=[]
	xyz=[]
	charge=[]
	for i in range(0,natoms):
		line = file.readline()
		data = line.split( )
		atom.append(data[0])	
		an.append(atom_tuple.index(atom[i]))
                atom_count[an[i]]+=1
		xyz.append([float(data[1]), float(data[2]), float(data[3])])
		charge.append(0.000)

############################################################################# DO SOMETHING

if 'cell' in locals():   #make uc ABC+abc if it was read in cell
  print
  print " ...converting cell (matrix) to CELL (ABCabc)"
  next=todo
elif 'ABC' in locals():  #make uc matrix if it was read in ABC+abc
  print
  print " ...converting CELL (ABCabc) to cell (matrix) "
  #alpha=B^C, beta=A^C, gamma=A^B
  cell=numpy.matrix([[                ABC[0],                     0.0,                                                                           0.0],
		     [ABC[1]*math.cos(abc[2]),ABC[1]*math.sin(abc[2]),                                                                           0.0],
		     [ABC[2]*math.cos(abc[1]),ABC[2]*math.cos(abc[0]),math.sqrt(ABC[2]**2-(ABC[2]*math.cos(abc[1]))**2-(ABC[2]*math.cos(abc[0]))**2)]]) #check this part



if 'fract' in locals(): #convert in cartesian
  print
  print " ...converting fractional coordinates in cartesian"
  xyz=[] 
  for i in range(0,natoms):
	x=fract[i][0]*cell.item((0,0))+fract[i][1]*cell.item((1,0))+fract[i][2]*cell.item((2,0))
	y=fract[i][1]*cell.item((0,1))+fract[i][1]*cell.item((1,1))+fract[i][2]*cell.item((2,1))
	z=fract[i][2]*cell.item((0,2))+fract[i][1]*cell.item((1,2))+fract[i][2]*cell.item((2,2))
	xyz.append([x,y,z])
elif 'xyz' in locals(): #convert in fractionals
  print
  print " ...converting cartesian coordinates in fractional"
  from numpy.linalg import inv
  invcell=inv(cell)
  fract=[]
  for i in range(0,natoms):
	x=xyz[i][0]*invcell.item((0,0))+xyz[i][1]*invcell.item((1,0))+xyz[i][2]*invcell.item((2,0))
	y=xyz[i][1]*invcell.item((0,1))+xyz[i][1]*invcell.item((1,1))+xyz[i][2]*invcell.item((2,1))
	z=xyz[i][2]*invcell.item((0,2))+xyz[i][1]*invcell.item((1,2))+xyz[i][2]*invcell.item((2,2))
	fract.append([x,y,z])


#reading what to do
outputfilename = sys.argv[2].split(".")[0]
outputformat= sys.argv[2].split(".")[1]

############################################################################## OUTPUT INFO
print
print "***************************************************"
print "  keep calm: I am converting %s to %s" % (inputformat, outputformat)
print "***************************************************"

for i in range(1,len(atom_count)):
	if atom_count[i] != 0:
		print('{0:>5} {1:3} atoms'.format(atom_count[i],atom_tuple[i]))

print " ---- --- ----- "
print('{0:>5} {1:3} atoms'.format(natoms,'tot'))
print
############################################################################## OUTPUT FILE


ofile=open('./'+outputfilename+'.'+outputformat, 'w+')

#writing a CIF file
if outputformat=="cif":
	print >> ofile, "data_crystal"
	print >> ofile, " "
	print >> ofile, "_cell_length_a    %.3f" %ABC[0]
	print >> ofile, "_cell_length_b    %.3f" %ABC[1]
	print >> ofile, "_cell_length_c    %.3f" %ABC[2]
	print >> ofile, "_cell_angle_alpha %.3f" %math.degrees(abc[0])
	print >> ofile, "_cell_angle_beta  %.3f" %math.degrees(abc[1])
	print >> ofile, "_cell_angle_gamma %.3f" %math.degrees(abc[2])
	print >> ofile, " "
	print >> ofile, "_symmetry_space_group_name_Hall 'P 1'"
	print >> ofile, "_symmetry_space_group_name_H-M  'P 1'"
	print >> ofile, " "
	print >> ofile, "loop_"
	print >> ofile, "_symmetry_equiv_pos_as_xyz"
	print >> ofile, " 'x,y,z' "
	print >> ofile, " "
	print >> ofile, "loop_"
	print >> ofile, "_atom_site_label"
	print >> ofile, "_atom_site_type_symbol"
	print >> ofile, "_atom_site_fract_x"
	print >> ofile, "_atom_site_fract_y"
	print >> ofile, "_atom_site_fract_z"
	print >> ofile, "_atom_site_charge"

	for i in range(0,natoms):	
        	label=atom[i]    #removed: label=atom[i]+"_"+str(i+1)
		print >> ofile, ('{0:10} {1:5} {2:>9.3f} {3:>9.3f} {4:>9.3f} {5:>9.5f}'.format(label,  atom[i], fract[i][0], fract[i][1], fract[i][2], charge[i]))


#writing a PDB file
if outputformat=="pdb":
	print >> ofile, ('CRYST1{0:>9.3f}{1:>9.3f}{2:>9.3f}{3:>7.2f}{4:>7.2f}{5:>7.2f} P 1           1'.format( ABC[0],ABC[1],ABC[2],math.degrees(abc[0]),math.degrees(abc[1]),math.degrees(abc[2]) ))
	for i in range(0,natoms):
		print >> ofile, "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" %("ATOM", i+1, atom[i],"", "XXX", "X", 1,"",fract[i][0],fract[i][1],fract[i][2],1.00, 0.00, atom[i], "")


#writing a CSSR file
if outputformat=="cssr":
   	print >> ofile, "                               %.3f  %.3f  %.3f"                    %(ABC[0],ABC[1],ABC[2])
	print >> ofile, "                %.3f   %.3f   %.3f   SPGR =  1 P 1         OPT = 1" %(math.degrees(abc[0]),math.degrees(abc[1]),math.degrees(abc[2]))
	print >> ofile, "%d   0"							     %(natoms)
	print >> ofile, "0 %s       : %s"                                                    %(inputfilename,inputfilename)
	for i in range(0,natoms):
		print >> ofile, "%4d %3s %8.3f %8.3f %8.3f    0  0  0  0  0  0  0  0  0.000"    %(i+1, atom[i], fract[i][0],fract[i][1],fract[i][2])
 

if outputformat=="xyz":
   	print >> ofile, "%d"   %(natoms)
	print >> ofile, "CELL: %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %(ABC[0],ABC[1],ABC[2],math.degrees(abc[0]),math.degrees(abc[1]),math.degrees(abc[2]))

	for i in range(0,natoms):
		print >> ofile, "%3s %8.3f %8.3f %8.3f "  %(atom[i], xyz[i][0],xyz[i][1],xyz[i][2])


