#! /usr/bin/python2.7

# domod: downloaded from http://wiki.wpi.edu/deskinsgroup/Density_of_States_-_CP2K on 12/6/17

#---------------------------------------------------
# get-smearing-pdos.py: read one or a pair alpha,  
# beta spin files with the cp2k pdos format and
# return a file "smeared.dat" with the Smeared DOS
#---------------------------------------------------
# Usage: ./get-smearing-pdos.py ALPHA.pdos BETA.pdos
#        or
#        ./get-smearing-pdos.py file.pdos 
#
# Output: 
#         smeared.dat: smeared DOS
#---------------------------------------------------
# Todo:
# - Atomatic name generation of output file
# - Move the algorithm to the module pdos
# - Implement printing of d orbitals
# - ...
#---------------------------------------------------
# Author: Juan Garcia e-mail: jcgarcia [at] wpi.edu
# Date:   11-12-2012
#---------------------------------------------------
import sys
from math import pi, sqrt	
import numpy as np

class pdos:
    """ Projected electronic density of states from CP2K output files

        Attributes
        ----------
        atom: str 
            the name of the atom where the DoS is projected
        iterstep: int
            the iteration step from the CP2K job
        efermi: float
            the energy of the Fermi level [a.u]
        e: float
            (eigenvalue - efermi) in eV
        occupation: int
            1 for occupied state or 0 for unoccupied
        pdos: nested list of float
            projected density of states on each orbital for each eigenvalue
            [[s1, p1, d1,....], [s2, p2, d2,...],...]
            s: pdos in s orbitals
            p: pdos in p orbitals
            d: pdos in d orbitals
            .
            .
            .
        tpdos: list of float
            sum of all the orbitals PDOS
            
        Methods
        -------
        smearing(self,npts, width)
            return the smeared tpdos 
    """
    
    def __init__(self, infilename):
        """Read a CP2K .pdos file and build a pdos instance

        Parameters
        ----------
        infilename: str
            pdos output from CP2K. 

        """    
        input_file = open(infilename, 'r')

        firstline  = input_file.readline().strip().split()
        secondline = input_file.readline().strip().split()


        if 
        # Kind of atom
        self.atom = firstline[6]   #domod: 7 instead of 6
        self.list = firstline[6]
        #iterationstep
        self.iterstep = int(firstline[12][:-1]) #[:-1] delete "," #domod: 14 instead of 12
        # Energy of the Fermi level
        self.efermi = float(firstline[15])   #domod: 17 instead 15

        # it keeps just the orbital names
        secondline[0:5] = []
        self.orbitals = secondline 

        lines = input_file.readlines()
   
        eigenvalue = []
        self.occupation = []
        data = []
        self.pdos = []
        for index, line in enumerate(lines):
            data.append(line.split())
            data[index].pop(0)
            eigenvalue.append(float(data[index].pop(0)))
            self.occupation.append(int(float(data[index].pop(0))))
            self.pdos.append([float(i) for i in data[index]])

        self.e = [(x-self.efermi)*27.211384523 for x in eigenvalue] 

        self.tpdos = []
        for i in self.pdos:
           self.tpdos.append(sum(i))

    def __add__(self, other):
        """Return the sum of two PDOS objects"""
        sumtpdos = [i+j for i,j in zip(self.tpdos,other.tpdos)]
        return sumtpdos

    def delta(self,emin,emax,npts,energy,width):
        """Return a delta-function centered at energy
        
        Parameters
        ----------
        emin: float
            minimun eigenvalue
        emax: float
            maximun eigenvalue
        npts: int
            Number of points in the smeared pdos
        energy: float
            energy where the gaussian is centered
        width: float
            dispersion parameter

        Return 
        ------
        delta: numpy array
            array of delta function values

        """
        
        energies = np.linspace(emin, emax, npts)
        x = -((energies - energy) / width)**2
        return np.exp(x) / (sqrt(pi) * width)

    def smearing(self,npts, width,):
        """Return a gaussian smeared DOS"""

        d = np.zeros(npts)
        emin = min(self.e)
        emax = max(self.e)
        for e, pd in zip(self.e,self.tpdos):
            d +=pd*self.delta(emin,emax,npts,e,width)

        return d
    
def sum_tpdos(tpdos1, tpdos2):
    """Return the sum of two PDOS"""
    return [i+j for i,j in zip(tpdos1,tpdos2)]


############################################################################


if len(sys.argv) == 2:

    infilename = sys.argv[1]

    alpha = pdos(infilename)
    npts = len(alpha.e)
    alpha_smeared = alpha.smearing(npts,0.2)
    eigenvalues = np.linspace(min(alpha.e), max(alpha.e),npts)
    
    g = open('smeared.dat','w')
    for i,j in zip(eigenvalues, alpha_smeared):
        t = str(i).ljust(15) + '     ' + str(j).ljust(15) + '\n'
        g.write(t)

elif len(sys.argv) == 3:

    infilename1 = sys.argv[1]
    infilename2 = sys.argv[2]

    alpha = pdos(infilename1)
    beta = pdos(infilename2)
    npts = len(alpha.e)
    alpha_smeared = alpha.smearing(npts,0.2)
    beta_smeared = beta.smearing(npts,0.2)
    totalDOS = sum_tpdos(alpha_smeared, beta_smeared)
    
    eigenvalues = np.linspace(min(alpha.e), max(alpha.e),npts)

    g = open('smeared.dat','w')
    for i,j in zip(eigenvalues, totalDOS):
        t = str(i).ljust(15) + '     ' + str(j).ljust(15) + '\n'
        g.write(t)

else:
    print '  Wrong number of arguments!'
    print '  usage:'
    print '  ./get-smearing-pdos.py ALPHA.pdos'
    print '  ./get-smearing-pdos.py ALPHA.pdos BETA.pdos'






