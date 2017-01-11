#!/usr/bin/env python

import math
import numpy as np

#xyz= [[x_atom1, y_atom1, z_atom1],
#      [x_atom2, y_atom2, z_atom2],
#       ...
#      [x_atomN, y_atomN, z_atomN]]  

#NB: all matrix and vectors are numpy.matrix and numpy.array. THEY ARE NOT PYTHON LISTs OR ARRAYs

# translate the system xyz to the origin. (if you want to traslate the origin to the point P, use origin=-P)
def transl(xyz,orig):            
	orig1=[]
	for i in range(0,xyz.shape[0]):
		orig1.append(orig)
	return xyz-orig

# make first a rotation of the xyz system in the x direction by xdeg, then y by ydeg, then z by zdeg. (right hand rule for the sense of rotation)
# the center of the rotation is center=0 for the origin (rotate the whole system).
# to rotate a OCO molecule keeping C as center of rotation: rotate_body(co2_xyz,co2xyz[1,:],xdeg,ydeg,zdeg) 
def rotate(xyz,center,xdeg,ydeg,zdeg):	  
	if np.size==1:
		center=np.array([0.,0.,0.])
		print "true"
	transl1=transl(xyz,center)
	xrad=math.radians(xdeg)
   	yrad=math.radians(ydeg)
   	zrad=math.radians(zdeg)
	Qx=np.matrix([[1.            ,0.             ,0.            ],[0.            ,math.cos(xrad),-math.sin(xrad)],[0.             ,math.sin(xrad),math.cos(xrad)]])
	Qy=np.matrix([[math.cos(yrad),0.             ,math.sin(yrad)],[0.            ,1.            ,0.             ],[-math.sin(yrad),0.            ,math.cos(yrad)]])
	Qz=np.matrix([[math.cos(zrad),-math.sin(zrad),0.            ],[math.sin(zrad),math.cos(zrad),0              ],[0.             ,0.            ,1.            ]])
	rotx  =np.transpose(Qx*np.transpose(transl1))
	rotxy =np.transpose(Qy*np.transpose(rotx   ))
	rotxyz=np.transpose(Qz*np.transpose(rotxy  ))
	xyznew=transl(rotxyz,-center)
	return xyznew

a=np.matrix([[0.00000, 0.00000,            3.616354],
	     [0.00000, 0.00000,       3.616354+1.17],
	     [0.00000, 0.00000, 3.616354+1.17+1.17]])

print " 30"
print " "

for i in range(0,91,10):
	#print i
	e=rotate(a,a[0,:],0,-i,45)
	e=transl(e,-1*np.array([0.385400,0.385400,0.]))
	#print "dist 0C = %f" %math.sqrt((e[0,0]-e[1,0])**2+(e[0,1]-e[1,1])**2+(e[0,2]-e[1,2])**2)

	print 'O %f %f %f ' %(e[0,0],e[0,1],e[0,2])
	print 'C %f %f %f ' %(e[1,0],e[1,1],e[1,2])
	print 'O %f %f %f ' %(e[2,0],e[2,1],e[2,2])






