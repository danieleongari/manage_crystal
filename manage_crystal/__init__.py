from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import re  #re.split(r'(\d+)',"Cu23") = ['Cu', '23', '']
import math
from collections import Counter  #makes a dictionary
from manage_crystal.atomic_data import atomic_symbol
from numpy.linalg import inv
from six.moves import range


class Crys:
    """ Crystal object """

    def __init__(self):
        # Initialize the molecule,
        self.natom = 0  # number of atoms
        self.atom_type = []  # atomic type: e.g., Cu1
        self.atom_element = []
        self.inp_xyz = False
        self.inp_fract = False
        self.atom_xyz = []  # cartesian coordinates
        self.atom_fract = []
        self.atom_charge = []
        self.tot_charge = []
        self.inp_lengths_angles = False
        self.inp_matrix = False
        self.length = [0.0] * 3
        self.angle_deg = [0.0] * 3
        self.angle_rad = [0.0] * 3
        self.matrix = [[0.0] * 3 for i in range(3)]

    def after_parse(self):
        # Input data check (coordinates and cell)
        if len(self.atom_xyz) > 0: self.inp_xyz = True
        if len(self.atom_fract) > 0: self.inp_fract = True
        if self.inp_xyz and self.inp_fract:
            sys.exit("WARNING: the input contains both Cartesian and \
                      fractional coordinates. EXIT.")
        elif not self.inp_xyz and not self.inp_fract:
            sys.exit("WARNING: no input coordinates. EXIT.")
        if all(x != 0 for x in self.length): self.inp_lengths_angles = True
        if not (x == 0 for x in self.matrix[0]): self.matrix = True
        if self.inp_lengths_angles and self.inp_matrix:
            sys.exit("WARNING: the input contains both lengths & angles \
                      and cell matrix. EXIT.")
        elif not self.compute_la_from_matrix and not self.inp_matrix:
            print("WARNING: no input cell.")
        # Computing stuff
        self.natom = len(self.atom_type)
        self.atom_element = [re.split(r'(\d+)', x)[0] for x in self.atom_type]
        self.element_count = Counter(self.atom_type)
        self.atom_atnum = [atomic_symbol.index(x) for x in self.atom_element]
        if all(x == 0 for x in self.angle_rad):
            self.angle_rad = [math.radians(i) for i in self.angle_deg]
        else:
            self.angle_deg = [math.degrees(i) for i in self.angle_rad]
        if len(self.atom_charge) == 0:
            self.atom_charge = [0] * self.natom

    def compute_la_from_matrix(self):
        self.length[0] = math.sqrt(self.matrix[0][0] * self.matrix[0][0] +
                                   self.matrix[0][1] * self.matrix[0][1] +
                                   self.matrix[0][2] * self.matrix[0][2])
        self.length[1] = math.sqrt(self.matrix[1][0] * self.matrix[1][0] +
                                   self.matrix[1][1] * self.matrix[1][1] +
                                   self.matrix[1][2] * self.matrix[1][2])
        self.length[2] = math.sqrt(self.matrix[2][0] * self.matrix[2][0] +
                                   self.matrix[2][1] * self.matrix[2][1] +
                                   self.matrix[2][2] * self.matrix[2][2])
        self.angle_rad[0] = math.acos(
            (self.matrix[1][0] * self.matrix[2][0] + self.matrix[1][1] *
             self.matrix[2][1] + self.matrix[1][2] * self.matrix[2][2]) /
            self.length[1] / self.length[2])  #alpha=B^C
        self.angle_rad[1] = math.acos(
            (self.matrix[0][0] * self.matrix[2][0] + self.matrix[0][1] *
             self.matrix[2][1] + self.matrix[0][2] * self.matrix[2][2]) /
            self.length[0] / self.length[2])  #beta=A^C
        self.angle_rad[2] = math.acos(
            (self.matrix[0][0] * self.matrix[1][0] + self.matrix[0][1] *
             self.matrix[1][1] + self.matrix[0][2] * self.matrix[1][2]) /
            self.length[0] / self.length[1])  #gamma=A^B
        self.angle_deg = [math.degrees(i) for i in self.angle_rad]

    def compute_matrix_from_la(self):
        #Copied from Raspa>framework.c>UnitCellBox
        self.matrix[0][0] = self.length[0]
        self.matrix[0][1] = 0.0
        self.matrix[0][2] = 0.0
        self.matrix[1][0] = self.length[1] * math.cos(self.angle_rad[2])
        self.matrix[1][1] = self.length[1] * math.sin(self.angle_rad[2])
        self.matrix[1][2] = 0.0
        self.matrix[2][0] = self.length[2] * math.cos(self.angle_rad[1])
        self.matrix[2][1] = self.length[2] * (
            math.cos(self.angle_rad[0]) - math.cos(self.angle_rad[2]) *
            math.cos(self.angle_rad[1])) / math.sin(self.angle_rad[2])
        self.matrix[2][2] = self.length[2] * math.sqrt(
            1 - (math.cos(self.angle_rad[1]))**2 -
            ((math.cos(self.angle_rad[0]) - math.cos(self.angle_rad[2]) *
              math.cos(self.angle_rad[1])) / math.sin(self.angle_rad[2]))**2)

    def compute_both_cell(self):
        if self.inp_lengths_angles: self.compute_matrix_from_la()
        elif self.inp_matrix: self.compute_la_from_matrix()
        self.invmatrix = inv(self.matrix)

    def compute_fract(self):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_fract = [[0.0] * 3 for i in range(self.natom)]
        for i in range(self.natom):
            for j in range(3):
                self.atom_fract[i][j] = self.atom_xyz[i][0] * self.invmatrix[0][j] + \
                                        self.atom_xyz[i][1] * self.invmatrix[1][j] + \
                                        self.atom_xyz[i][2] * self.invmatrix[2][j]

    def compute_xyz(self):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_xyz = [[0.0] * 3 for i in range(self.natom)]
        for i in range(self.natom):
            for j in range(3):
                self.atom_xyz[i][j] = self.atom_fract[i][0] * self.matrix[0][j] + \
                                      self.atom_fract[i][1] * self.matrix[1][j] + \
                                      self.atom_fract[i][2] * self.matrix[2][j]

    def compute_both_coord(self):
        if self.inp_xyz: self.compute_fract()
        elif self.inp_fract: self.compute_xyz()
