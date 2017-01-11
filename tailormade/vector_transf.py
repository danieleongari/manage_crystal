#!/usr/bin/env python

import math
import numpy as np

def rotate(xyz,xdeg,ydeg,zdeg):	
	xrad=math.radians(xdeg)
   	yrad=math.radians(ydeg)
   	zrad=math.radians(zdeg)
	Qx=np.matrix([[1.            ,0.             ,0.            ],[0.            ,math.cos(xrad),-math.sin(xrad)],[0.             ,math.sin(xrad),math.cos(xrad)]])
	Qy=np.matrix([[math.cos(yrad),0.             ,math.sin(yrad)],[0.            ,1.            ,0.             ],[-math.sin(yrad),0.            ,math.cos(yrad)]])
	Qz=np.matrix([[math.cos(zrad),-math.sin(zrad),0.            ],[math.sin(zrad),math.cos(zrad),0              ],[0.             ,0.            ,1.            ]])
	xyz1=np.transpose(Qx*np.transpose(xyz))
	xyz2=np.transpose(Qy*np.transpose(xyz1))
	xyz3=np.transpose(Qz*np.transpose(xyz2))
	return xyz3


def transl(orig, xyz):
	orig1=[]
	for i in range(0,xyz.shape[0]):
		orig1.append(orig)
	return xyz-orig
  

a=[[0.00000, 0.00000,3.679062],[0.00000, 0.00000, 4.849062],[0.00000, 0.00000, 6.019062]]
a=np.matrix(a)

orig=a[0,:]

for i in 30, 45, 60, 90, 120:
	b= transl(orig,a)
	c=rotate(b,i,0,45)
	e= transl(-orig,c)
	#rint (e[0,0]-e[1,0])**2+(e[0,1]-e[1,1])**2+(e[0,2]-e[1,2])**2

	print 'O %f %f %f ' %(e[0,0],e[0,1],e[0,2])
	print 'C %f %f %f ' %(e[1,0],e[1,1],e[1,2])
	print 'O %f %f %f ' %(e[2,0],e[2,1],e[2,2])


for i in 30, 45, 60, 90:
	b= transl(orig,a)
	c=rotate(b,i,0,0)
	e= transl(-orig,c)
	#rint (e[0,0]-e[1,0])**2+(e[0,1]-e[1,1])**2+(e[0,2]-e[1,2])**2

	print 'O %f %f %f ' %(e[0,0],e[0,1],e[0,2])
	print 'C %f %f %f ' %(e[1,0],e[1,1],e[1,2])
	print 'O %f %f %f ' %(e[2,0],e[2,1],e[2,2])

def transl(orig, xyz):
	orig1=[]
	for i in range(0,xyz.shape[0]):
		orig1.append(orig)
	return xyz-orig
  

a=[[0.00000, 0.00000,3.679062+0.4],[0.00000, 0.00000, 4.849062+0.4],[0.00000, 0.00000, 6.019062+0.4]]
a=np.matrix(a)

orig=a[0,:]

for i in 30, 45, 60, 90, 120:
	b= transl(orig,a)
	c=rotate(b,i,0,45)
	e= transl(-orig,c)
	#rint (e[0,0]-e[1,0])**2+(e[0,1]-e[1,1])**2+(e[0,2]-e[1,2])**2

	print 'O %f %f %f ' %(e[0,0],e[0,1],e[0,2])
	print 'C %f %f %f ' %(e[1,0],e[1,1],e[1,2])
	print 'O %f %f %f ' %(e[2,0],e[2,1],e[2,2])




