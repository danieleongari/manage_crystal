from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import re  #re.split(r'(\d+)',"Cu23") = ['Cu', '23', '']
import math
from collections import Counter  #makes a dictionary
from manage_crystal.periodic_table import ptab_atnum, ptab_mass
from numpy.linalg import inv
from six.moves import range

AVOGCONST = 6.022E+23
M3TOANG3 = 1e10**3
GTOKG = 1 / 1000


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
        self.angle_deg = [0.0] * 3  # Main
        self.angle_rad = [0.0] * 3  # Use ONLY when necessary to simplify math
        self.matrix = [[0.0] * 3 for i in range(3)]

    def check_parse(self):
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

    def compute_atom_count(self):
        self.natom = len(self.atom_type)
        self.atom_element = [re.split(r'(\d+)', x)[0] for x in self.atom_type]
        self.element_count = Counter(self.atom_type)
        self.element = sorted(self.element_count)
        self.nelement = len(self.element)
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
        self.angle_deg = [math.degrees(x) for x in self.angle_rad]

    def compute_matrix_from_la(self):
        #Copied from Raspa>framework.c>UnitCellBox
        self.angle_rad = [math.radians(x) for x in self.angle_deg]
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

    def compute_fract_from_xyz(self):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_fract = [[0.0] * 3 for i in range(self.natom)]
        for i in range(self.natom):
            for j in range(3):
                self.atom_fract[i][j] = self.atom_xyz[i][0] * self.invmatrix[0][j] + \
                                        self.atom_xyz[i][1] * self.invmatrix[1][j] + \
                                        self.atom_xyz[i][2] * self.invmatrix[2][j]

    def compute_xyz_from_fract(self):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_xyz = [[0.0] * 3 for i in range(self.natom)]
        for i in range(self.natom):
            for j in range(3):
                self.atom_xyz[i][j] = self.atom_fract[i][0] * self.matrix[0][j] + \
                                      self.atom_fract[i][1] * self.matrix[1][j] + \
                                      self.atom_fract[i][2] * self.matrix[2][j]

    def compute_both_coord(self):
        if self.inp_xyz: self.compute_fract_from_xyz()
        elif self.inp_fract: self.compute_xyz_from_fract()

    def transl_coord(self, transl):
        """ Translate the coordinates of the crystal by a defined amount """
        for i in range(self.natom):
            self.atom_xyz[i][0] += transl[0]
            self.atom_xyz[i][1] += transl[1]
            self.atom_xyz[i][2] += transl[2]
        self.compute_fract_from_xyz()

    def randomize_coord(self, delta):
        """ Randomize the atomic coordinates, by a normal distribution """
        for i in range(self.natom):
            self.atom_xyz[i][0] += np.random.normal(0, delta, 1)
            self.atom_xyz[i][1] += np.random.normal(0, delta, 1)
            self.atom_xyz[i][2] += np.random.normal(0, delta, 1)
        self.compute_fract_from_xyz()

    def rotate_axis(self, up):
        if up:
            self.length[0], self.length[1], self.length[2] = \
             self.length[2], self.length[0], self.length[1]
            self.angle_deg[0], self.angle_deg[1], self.angle_deg[2] = \
             self.angle_deg[2], self.angle_deg[0], self.angle_deg[1]
            for i in range(self.natom):
                self.atom_xyz[i][0], self.atom_xyz[i][1], \
                 self.atom_xyz[i][2] = self.atom_xyz[i][2], \
                  self.atom_xyz[i][0], self.atom_xyz[i][1]
                self.atom_fract[i][0], self.atom_fract[i][1], \
                 self.atom_fract[i][2] = self.atom_fract[i][2], \
                  self.atom_fract[i][0], self.atom_fract[i][1]
        else:
            self.length[0], self.length[1], self.length[2] = \
             self.length[1], self.length[2], self.length[0]
            self.angle_deg[0], self.angle_deg[1], self.angle_deg[2] = \
             self.angle_deg[1], self.angle_deg[2], self.angle_deg[0]
            for i in range(self.natom):
                self.atom_xyz[i][0], self.atom_xyz[i][1], \
                 self.atom_xyz[i][2] = self.atom_xyz[i][1], \
                  self.atom_xyz[i][2], self.atom_xyz[i][0]
                self.atom_fract[i][0], self.atom_fract[i][1], \
                 self.atom_fract[i][2] = self.atom_fract[i][1], \
                  self.atom_fract[i][2], self.atom_fract[i][0]

    def compute_perp_width(self):
        """ Compute perpendicular widths in the cell """
        ax = crys.matrix[0][0]
        ay = crys.matrix[0][1]
        az = crys.matrix[0][2]
        bx = crys.matrix[1][0]
        by = crys.matrix[1][1]
        bz = crys.matrix[1][2]
        cx = crys.matrix[2][0]
        cy = crys.matrix[2][1]
        cz = crys.matrix[2][2]
        # calculate vector products of cell vectors
        axb1 = ay * bz - az * by
        axb2 = az * bx - ax * bz
        axb3 = ax * by - ay * bx
        bxc1 = by * cz - bz * cy
        bxc2 = bz * cx - bx * cz
        bxc3 = bx * cy - by * cx
        cxa1 = cy * az - ay * cz
        cxa2 = ax * cz - az * cx
        cxa3 = ay * cx - ax * cy
        # calculate volume of cell
        V = math.fabs(ax * bxc1 + ay * bxc2 + az * bxc3)
        # calculate cell perpendicular widths
        self.perp_width = [0.0] * 3
        self.perp_width[0] = V / math.sqrt(bxc1**2 + bxc2**2 + bxc3**2)
        self.perp_width[1] = V / math.sqrt(cxa1**2 + cxa2**2 + cxa3**2)
        self.perp_width[2] = V / math.sqrt(axb1**2 + axb2**2 + axb3**2)

    def dist_ij(self, i, j):
        """ Compute the distance (Å) between the atoms i and j """
        dist_fract = [(self.atom_fract[i][k] - self.atom_fract[j][k])
                      for k in range(3)]
        dist_fract_pbc = [ (dist_fract[k] - int(round(dist_fract[k]))) \
                        for k in range(3) ]
        dist_cart_pbc = [ ( self.matrix[0][k] * dist_fract_pbc[0] + \
                            self.matrix[1][k] * dist_fract_pbc[1] + \
                            self.matrix[2][k] * dist_fract_pbc[2] ) \
                            for k in range(3) ]
        dist = math.sqrt(dist_cart_pbc[0]**2 + \
                         dist_cart_pbc[1]**2 + \
                         dist_cart_pbc[2]**2)
        return dist

    def expand_k_dir(self, k, n):
        """ Expand and replicate the cell n times in k direction """
        for i in range(1, n):
            for j in range(self.natom):
                self.atom_type.append(self.atom_type[j])
                self.atom_charge.append(self.atom_charge[j])
                self.atom_xyz.append([
                    self.atom_xyz[j][0] + i * self.matrix[k][0],
                    self.atom_xyz[j][1] + i * self.matrix[k][1],
                    self.atom_xyz[j][2] + i * self.matrix[k][2]
                ])
        self.length[k] *= n
        for kk in range(3):
            self.matrix[k][kk] *= n
        self.compute_fract_from_xyz()
        self.compute_atom_count()

    def compute_volume_from_la(self):
        """ Compute cell Volume from lengths and angles """
        # www.fxsolver.com/browse/formulas/Triclinic+crystal+system+(Unit+cell's+volume)
        self.angle_rad = [math.radians(x) for x in self.angle_deg]
        vol = self.length[0] * self.length[1] * self.length[2] * \
              math.sqrt(1 - math.cos(self.angle_rad[0])**2 - \
                            math.cos(self.angle_rad[1])**2 - \
                            math.cos(self.angle_rad[2])**2 + \
                            2 * math.cos(self.angle_rad[0]) * \
                                math.cos(self.angle_rad[1]) * \
                                math.cos(self.angle_rad[2])
                       )
        return vol

    def compute_weight_density(self):
        """ Compute crystal density """
        weight = 0  #g/mol_uc
        for element in self.element_count:
            weight += self.element_count[element] * ptab_mass[element]
        vol = self.compute_volume_from_la()
        rho_kgm3 = weight / vol / AVOGCONST * M3TOANG3 * GTOKG
        return weight, rho_kgm3

    def compute_nelectron(self):
        """ Compute the number of electrons """
        nelectron = 0
        for element in self.element_count:
            nelectron += self.element_count[element] * ptab_atnum[element]
        return nelectron
