import numpy as np
import re  #re.split(r'(\d+)',"Cu23") = ['Cu', '23', '']
import math
from collections import Counter #makes a dictionary
from manage_crystal.atomic_data import atomic_symbol

class Crys:
    """ Crystal object """
    def __init__(self):
        # Initialize the molecule,
        self.natom = 0       # number of atoms
        self.atom_type = []  # atomic type: e.g., Cu1
        self.atom_element = []
        self.inp_xyz = False
        self.inp_fract = False
        self.atom_xyz = []   # cartesian coordinates
        self.atom_fract = []
        self.atom_charge = []
        self.tot_charge = []
        self.inp_lengths_angles = False
        self.inp_matrix = False
        self.length = np.zeros(3)
        self.angle_deg = np.zeros(3)
        self.angle_rad = np.zeros(3)
        self.matrix = np.zeros((3,3))
    def after_parse(self):
        self.natom = len(self.atom_type)
        self.atom_element = [ re.split(r'(\d+)',x)[0] for x in self.atom_type ]
        self.element_count = Counter(self.atom_type)
        self.atom_atnum = [ atomic_symbol.index(x) for x in self.atom_element ]
        if all(x ==0 for x in self.angle_rad):
            self.angle_rad = [math.radians(i) for i in self.angle_deg]
        else:
            self.angle_deg = [math.degrees(i) for i in self.angle_rad]

    def compute_la_from_matrix(self):
        self.length[0] = math.sqrt(self.matrix[0][0]*self.matrix[0][0]+self.matrix[0][1]*self.matrix[0][1]+self.matrix[0][2]*self.matrix[0][2])
        self.length[1] = math.sqrt(self.matrix[1][0]*self.matrix[1][0]+self.matrix[1][1]*self.matrix[1][1]+self.matrix[1][2]*self.matrix[1][2])
        self.length[2] = math.sqrt(self.matrix[2][0]*self.matrix[2][0]+self.matrix[2][1]*self.matrix[2][1]+self.matrix[2][2]*self.matrix[2][2])
        self.angle_rad[0] = math.acos((self.matrix[1][0]*self.matrix[2][0]+self.matrix[1][1]*self.matrix[2][1]+self.matrix[1][2]*self.matrix[2][2])/self.length[1]/self.length[2]) #alpha=B^C
        self.angle_rad[1] = math.acos((self.matrix[0][0]*self.matrix[2][0]+self.matrix[0][1]*self.matrix[2][1]+self.matrix[0][2]*self.matrix[2][2])/self.length[0]/self.length[2]) #beta=A^C
        self.angle_rad[2] = math.acos((self.matrix[0][0]*self.matrix[1][0]+self.matrix[0][1]*self.matrix[1][1]+self.matrix[0][2]*self.matrix[1][2])/self.length[0]/self.length[1]) #gamma=A^B
        self.angle_deg = [math.degrees(i) for i in self.angle_rad]
    def compute_matrix_from_la(self):
        self.matrix[0][0] = self.length[0]
        self.matrix[0][1] = 0.0
        self.matrix[0][2] = 0.0
        self.matrix[1][0] = self.length[1]*math.cos(self.angle_rad[2])
        self.matrix[1][1] = self.length[1]*math.sin(self.angle_rad[2])
        self.matrix[1][2] = 0.0
        self.matrix[2][0] = self.length[2]*math.cos(self.angle_rad[1])
        self.matrix[2][1] = self.length[2]*(math.cos(self.angle_rad[0])-math.cos(self.angle_rad[2])*math.cos(self.angle_rad[1]))/math.sin(self.angle_rad[2])
        self.matrix[2][2] = self.length[2]*math.sqrt(1-(math.cos(self.angle_rad[1]))**2-((math.cos(self.angle_rad[0])-math.cos(self.angle_rad[2])*math.cos(self.angle_rad[1]))/math.sin(self.angle_rad[2]))**2)
    def compute_all_info(self):
        # Compute several thing after the unit cell is specified
        if self.inp_lengths_angles == True:
            self.compute_matrix_from_la()
        if self.inp_matrix == True:
            self.compute_matrix_from_la()
        self.invmatrix=inv(self.matrix)
