#!/usr/bin/python2.7
"""/
Python program to print thermodynamic infos.
"""
import string,sys
import numpy
import math
import subprocess
import matplotlib.pyplot as plt

#Constants
h=6.62606957e-34      # J*s       J= 1kg*m^2/s^2
hbar=h/(2.0*math.pi)  # J*s
R=8.314472            # J/K/mol
Navo=6.022e+23        # molec/mol
kb=R/Navo             # J/K/molec

mass=17.0             # g/mol
mass=mass/Navo/1000.  # Kg/molec

T=298.0      #K 
P=1.000      #bar
P=P*1e+5     #Pa = 1kg/m^2

phase="liq"
conc_liq=0.1 #mol/liter


if phase=="id_gas":
   conc=P/(kb*T)  #molec/m^3 ***ideal gas
elif phase=="liq":
   conc=conc_liq*Navo*1000 #molec/m^3


q_trasl= ((2.*math.pi*mass*kb*T/h**2)**(3./2.)) * (1/conc)    # adim        ref:Mendeley/Thermochemistry in Gaussian
S_trasl= R*(math.log(q_trasl)+1+3./2.)                    # J/K/mol     ref:Mendeley/Thermochemistry in Gaussian
ST_trasl= S_trasl*T                                             #J/mol
