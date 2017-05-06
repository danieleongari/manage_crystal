#!/usr/bin/python2.7

import string,sys
import numpy
import math
import subprocess

# Python program to compute the shorter distance between two points in a triclinic cell

############################################ input
A=[0.3,0.1,0.2]
B=[1.9,2.8,7.8]

cell=numpy.matrix([[4.0,0  ,0  ],
		   [3.0,7.0,0  ],
		   [2.0,3.0,8.0]])

############################################ initialize
d=[0]*3 #distance in xyz
s=[0]*3 #distance in abc
t=[0]*3 #distance in abc, with PBC

############################################ start
for i in range (0,3):
 d[i]=A[i]-B[i]

mindist=math.sqrt(d[0]**2+d[1]**2+d[2]**2)
print "NO PDB: dx,dy,dz= %8.3f %8.3f %8.3f | mindist=%8.3f" %(d[0],d[1],d[2],mindist)

from numpy.linalg import inv
invcell=inv(cell)

s[0]=invcell.item((0,0))*d[0]+invcell.item((1,0))*d[1]+invcell.item((2,0))*d[2]
s[1]=invcell.item((0,1))*d[0]+invcell.item((1,1))*d[1]+invcell.item((2,1))*d[2]
s[2]=invcell.item((0,2))*d[0]+invcell.item((1,2))*d[1]+invcell.item((2,2))*d[2]

t[0]=s[0]-int(round(s[0]))
t[1]=s[1]-int(round(s[1]))
t[2]=s[2]-int(round(s[2]))

d[0]=cell.item((0,0))*t[0]+cell.item((1,0))*t[1]+cell.item((2,0))*t[2]
d[1]=cell.item((0,1))*t[0]+cell.item((1,1))*t[1]+cell.item((2,1))*t[2]
d[2]=cell.item((0,2))*t[0]+cell.item((1,2))*t[1]+cell.item((2,2))*t[2]

mindist=math.sqrt(d[0]**2+d[1]**2+d[2]**2)
print "W/ PDB: dx,dy,dz= %8.3f %8.3f %8.3f | mindist=%8.3f" %(d[0],d[1],d[2],mindist)


################################################ print
ofile=open("mindist.xsf", 'w+')

print >> ofile, "CRYSTAL"
print >> ofile, "PRIMVEC 1"
print >> ofile, "     %8.5f %8.5f %8.5f"    %(cell.item((0,0)),cell.item((0,1)),cell.item((0,2)))
print >> ofile, "     %8.5f %8.5f %8.5f"    %(cell.item((1,0)),cell.item((1,1)),cell.item((1,2)))
print >> ofile, "     %8.5f %8.5f %8.5f"    %(cell.item((2,0)),cell.item((2,1)),cell.item((2,2)))
print >> ofile, "PRIMCOORD 1"
print >> ofile, "%d 1"                      %(2)

print >> ofile, "%3s %8.3f %8.3f %8.3f "  %('H ', A[0], A[1], A[2])
print >> ofile, "%3s %8.3f %8.3f %8.3f "  %('He', B[0], B[1], B[2])

ofile.close()

subprocess.call("/home/daniele/Programs/VESTA-x86_64/VESTA mindist.xsf &", shell=True)
subprocess.call("sleep 1; rm mindist.xsf", shell=True)

############################################################## computing the sphere overlap volume
#http://mathworld.wolfram.com/Sphere-SphereIntersection.html
#WE exclude that a sphere overlap with itself! (cell dimension < sphere diameter)

rA=1.0
rB=0.2

V_A=4/3*math.pi*rA**3
V_B=4/3*math.pi*rB**3
if (rA+rB<mindist): V_ovlp=0
else: V_ovlp=math.pi*(rA+rB-mindist)**2 * (mindist**2+2*mindist*rB-3*rB**2+2*mindist*rA+6*rB*rA-3*rA**2) / (12*mindist)

print "V_A, V_B, V_ovlp:  %11.7f %11.7f %11.7f" %(V_A,V_B,V_ovlp)

"""/ BUG:
A=[2.5,0.5,0.0]
B=[6.0,1.5,0.0]
cell=numpy.matrix([[4.0,0.0,0.0],
		   [3.0,2.0,0.0],
                   [0.0,0.0,1.0]])

# It is not working because it is the limiting case where the distance point
# is exactly 0.5 fractional coordinates and it is not shifted.
"""




""" From RASPA potentials.h:
static inline VECTOR ApplyBoundaryCondition(VECTOR dr)
{
  VECTOR s,t;

  switch(BoundaryCondition[CurrentSystem])
  {
    case FINITE:
      break;
    case RECTANGULAR:
    case CUBIC:
      dr.x-=Box[CurrentSystem].ax*(REAL)NINT(dr.x*InverseBox[CurrentSystem].ax);
      dr.y-=Box[CurrentSystem].by*(REAL)NINT(dr.y*InverseBox[CurrentSystem].by);
      dr.z-=Box[CurrentSystem].cz*(REAL)NINT(dr.z*InverseBox[CurrentSystem].cz);
      break;
    case TRICLINIC:
      // convert from xyz to abc
      s.x=InverseBox[CurrentSystem].ax*dr.x+InverseBox[CurrentSystem].bx*dr.y+InverseBox[CurrentSystem].cx*dr.z;
      s.y=InverseBox[CurrentSystem].ay*dr.x+InverseBox[CurrentSystem].by*dr.y+InverseBox[CurrentSystem].cy*dr.z;
      s.z=InverseBox[CurrentSystem].az*dr.x+InverseBox[CurrentSystem].bz*dr.y+InverseBox[CurrentSystem].cz*dr.z;

      // apply boundary condition
      t.x=s.x-(REAL)NINT(s.x);
      t.y=s.y-(REAL)NINT(s.y);
      t.z=s.z-(REAL)NINT(s.z);

      // convert from abc to xyz
      dr.x=Box[CurrentSystem].ax*t.x+Box[CurrentSystem].bx*t.y+Box[CurrentSystem].cx*t.z;
      dr.y=Box[CurrentSystem].ay*t.x+Box[CurrentSystem].by*t.y+Box[CurrentSystem].cy*t.z;
      dr.z=Box[CurrentSystem].az*t.x+Box[CurrentSystem].bz*t.y+Box[CurrentSystem].cz*t.z;
      break;
    default:
      fprintf(stderr,"Error: Unkown boundary condition....\n");
      exit(0);
      break;
  }
  return dr;
}
"""
