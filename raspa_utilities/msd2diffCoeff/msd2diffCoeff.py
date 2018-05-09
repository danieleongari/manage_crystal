#!/usr/bin/env python3

import math 
import scipy 
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import argparse 
from argparse import RawTextHelpFormatter #needed to go next line in the help text

parser = argparse.ArgumentParser(description='Program to compute the diffusion coefficient from Raspa\'s MSD (by Daniele Ongari)', formatter_class=RawTextHelpFormatter)

parser.add_argument("inputfile", 
                      type=str,
                      help="MSD file msd_self_molecule_n.dat")

parser.add_argument("-timerange", 
                      action="store", 
                      type=float,
                      nargs=2, 
                      dest="timerange",
                      default=[1.0,16.0],
                      help="tmin tmax [ps] for the fitting")

args = parser.parse_args()

data=np.genfromtxt(args.inputfile, delimiter="", comments="#",usecols = (0,1,2,3,4))

[rows,columns]=np.shape(data)

time_fit=[]
msd_fit=[]

for i in range(0,rows):
	if args.timerange[0] <= data[i][0] <= args.timerange[1]:
		time_fit.append(data[i][0])
		msd_fit.append(data[i][1])


#Stuff to make the arrays sklearn wants
time_fit = np.asarray(time_fit).reshape(-1, 1)
msd_fit  = np.asarray(msd_fit).reshape(-1, 1)

regr_msd = linear_model.LinearRegression()
regr_msd.fit(time_fit, msd_fit)
msd_pred = regr_msd.predict(time_fit)
msd_r2 = r2_score(msd_fit, msd_pred) 

diffcoeff=regr_msd.coef_ /(2*3)     #angs^2/ps
diffcoeff=diffcoeff*(1e-12)/(1e-20) #m2/s 1A^2=1e-20m2, 1ps=1e-12s


print('Diffusion_coefficient_(m^2/s): \t %.3E' %diffcoeff)
print('R^2_score: \t \t \t %.3f' %msd_r2)



