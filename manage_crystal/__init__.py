import numpy as np
import re  #re.split(r'(\d+)',"Cu23") = ['Cu', '23', '']
import math
import os
import sys
from collections import Counter  #makes a dictionary
from manage_crystal.periodic_table import ptab_atnum, ptab_atnum_inv, ptab_mass, ptab_vdw_uff
from numpy.linalg import inv

AVOGCONST = 6.022E+23
M3TOANG3 = 1e10**3
GTOKG = 0.001


class Crys:
    """ Crystal object """

    def __init__(self):
        # Initialize the molecule,
        self.natom = 0  # number of atoms
        self.atom_type = []  # atomic type: e.g., Cu1
        self.inp_xyz = False
        self.inp_fract = False
        self.atom_xyz = []  # cartesian coordinates
        self.atom_fract = []
        self.atom_charge = []
        self.tot_charge = None
        self.inp_lengths_angles = False
        self.inp_matrix = False
        self.length = [0.0] * 3
        self.angle_deg = [0.0] * 3  # Main
        self.angle_rad = [0.0] * 3  # Use ONLY when necessary to simplify math
        self.matrix = [[0.0] * 3 for i in range(3)]

    def check_matrix_allign(self):
        """ Checks if the cell matrix is alligned with xyz axis """
        # This is important because in an MD the cell can get disalligned and,
        # if atomic coordinates are cartesian (xyz), cell and coordinates would
        # be referred to difference reference systems.
        if all(x == 0 for x in [self.matrix[0][1], self.matrix[0][2], self.matrix[1][2]]):
            self.matrix_alligned = True
        else:
            self.matrix_alligned = False

    def check_parse(self):
        """ Checks input coordinates and cell"""
        # Check atomic coordinates
        if len(self.atom_xyz) > 0:
            self.inp_xyz = True
        if len(self.atom_fract) > 0:
            self.inp_fract = True
        if self.inp_xyz and self.inp_fract:
            sys.exit("WARNING: the input contains both Cartesian and " + "fractional coordinates. EXIT.")
        elif not self.inp_xyz and not self.inp_fract:
            sys.exit("WARNING: no input coordinates. EXIT.")
        # Check input cell
        if any(x != 0 for x in self.length):
            self.inp_lengths_angles = True
            self.matrix_alligned = True
            if any(x == 0 for x in self.length):
                sys.exit("WARNING: found cell length equal to zero. EXIT.")
        if any(x != 0 for x in self.matrix[0]):
            self.inp_matrix = True
            self.check_matrix_allign()
        if self.inp_lengths_angles and self.inp_matrix:
            sys.exit("WARNING: the input contains both lengths & angles " + "and cell matrix. EXIT.")
        elif not self.compute_la_from_matrix and not self.inp_matrix:
            print("WARNING: no input cell.")

    def compute_atom_count(self):
        """ Extract basic info from the atom_type list:
        natoms: int, number of atoms
        atom_element: str, element of the atom type (stripped from numbers)
        element_count: dict, element : number of atoms for that element
        element_list: list, list of elements, sorted by atomic number
        nelement: int, number of different elements

        Also, assign 0 charges to the atoms if not already assigned.
        """

        self.natom = len(self.atom_type)
        self.atom_element = [re.split(r'(\d+)', x)[0] for x in self.atom_type]
        self.element_count = Counter(self.atom_element)
        self.element_list = []
        for an in range(1, 118):
            element = ptab_atnum_inv[an]
            if element in self.element_count:
                self.element_list.append(element)
        self.nelement = len(self.element_list)
        if len(self.atom_charge) == 0:
            self.atom_charge = [0] * self.natom

    def remove_atoms(self, atoms_list):
        """ Remove atoms, given a list of indexes"""
        # exclude duplicates
        atoms_list = list(dict.fromkeys(atoms_list))

        # Remove atoms from list
        for i in sorted(atoms_list, reverse=True):
            del self.atom_type[i]
            del self.atom_xyz[i]
            del self.atom_fract[i]
            del self.atom_charge[i]

        # Remove atoms from sum and recompute counting
        self.compute_atom_count()

    def compute_la_from_matrix(self):
        self.length[0] = math.sqrt(self.matrix[0][0] * self.matrix[0][0] + self.matrix[0][1] * self.matrix[0][1] +
                                   self.matrix[0][2] * self.matrix[0][2])
        self.length[1] = math.sqrt(self.matrix[1][0] * self.matrix[1][0] + self.matrix[1][1] * self.matrix[1][1] +
                                   self.matrix[1][2] * self.matrix[1][2])
        self.length[2] = math.sqrt(self.matrix[2][0] * self.matrix[2][0] + self.matrix[2][1] * self.matrix[2][1] +
                                   self.matrix[2][2] * self.matrix[2][2])
        self.angle_rad[0] = math.acos(
            (self.matrix[1][0] * self.matrix[2][0] + self.matrix[1][1] * self.matrix[2][1] +
             self.matrix[1][2] * self.matrix[2][2]) / self.length[1] / self.length[2])  #alpha=B^C
        self.angle_rad[1] = math.acos(
            (self.matrix[0][0] * self.matrix[2][0] + self.matrix[0][1] * self.matrix[2][1] +
             self.matrix[0][2] * self.matrix[2][2]) / self.length[0] / self.length[2])  #beta=A^C
        self.angle_rad[2] = math.acos(
            (self.matrix[0][0] * self.matrix[1][0] + self.matrix[0][1] * self.matrix[1][1] +
             self.matrix[0][2] * self.matrix[1][2]) / self.length[0] / self.length[1])  #gamma=A^B
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
        self.matrix[2][1] = self.length[2] * (math.cos(self.angle_rad[0]) - math.cos(self.angle_rad[2]) *
                                              math.cos(self.angle_rad[1])) / math.sin(self.angle_rad[2])
        self.matrix[2][2] = self.length[2] * math.sqrt(1 - (math.cos(self.angle_rad[1]))**2 - (
            (math.cos(self.angle_rad[0]) - math.cos(self.angle_rad[2]) * math.cos(self.angle_rad[1])) /
            math.sin(self.angle_rad[2]))**2)

    def compute_fract_from_xyz(self):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_fract = [[0.0] * 3 for i in range(self.natom)]
        self.invmatrix = inv(self.matrix)
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

    def fix_cell_notalligned(self):
        """ Fix the problem of not alligned cell and cartesian coordinates """
        # matrix > la > matrix(alligned)
        # atom_fract are still valid
        # atom_xyz should be overwritten by the atomic_fract in the new cell
        if self.inp_matrix and self.inp_fract and not self.matrix_alligned:
            self.compute_matrix_from_la()
            self.compute_xyz_from_fract()

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
            self.atom_xyz[i][0] += np.random.normal(0, delta, 1)[0]
            self.atom_xyz[i][1] += np.random.normal(0, delta, 1)[0]
            self.atom_xyz[i][2] += np.random.normal(0, delta, 1)[0]
        self.compute_fract_from_xyz()

    def rotate_axis(self, up):
        """ Perform a rotation of the axes:
        'up' (up==True, XYZ to ZXY) or 'down' (up == False, XYZ to YZX) """
        if up:
            self.length[0], self.length[1], self.length[2] = \
             self.length[2], self.length[0], self.length[1]
            self.angle_deg[0], self.angle_deg[1], self.angle_deg[2] = \
             self.angle_deg[2], self.angle_deg[0], self.angle_deg[1]
            for i in range(self.natom):
                self.atom_xyz[i][0], self.atom_xyz[i][1], \
                 self.atom_xyz[i][2] = self.atom_xyz[i][2], \
                  self.atom_xyz[i][0], self.atom_xyz[i][1]

        else:
            self.length[0], self.length[1], self.length[2] = \
             self.length[1], self.length[2], self.length[0]
            self.angle_deg[0], self.angle_deg[1], self.angle_deg[2] = \
             self.angle_deg[1], self.angle_deg[2], self.angle_deg[0]
            for i in range(self.natom):
                self.atom_xyz[i][0], self.atom_xyz[i][1], \
                 self.atom_xyz[i][2] = self.atom_xyz[i][1], \
                  self.atom_xyz[i][2], self.atom_xyz[i][0]
        self.compute_matrix_from_la()
        self.compute_fract_from_xyz()

    def compute_perp_width(self):
        """ Compute perpendicular widths in the cell """
        ax = self.matrix[0][0]
        ay = self.matrix[0][1]
        az = self.matrix[0][2]
        bx = self.matrix[1][0]
        by = self.matrix[1][1]
        bz = self.matrix[1][2]
        cx = self.matrix[2][0]
        cy = self.matrix[2][1]
        cz = self.matrix[2][2]
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
        dist_fract = [(self.atom_fract[i][k] - self.atom_fract[j][k]) for k in range(3)]
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
                    self.atom_xyz[j][0] + i * self.matrix[k][0], self.atom_xyz[j][1] + i * self.matrix[k][1],
                    self.atom_xyz[j][2] + i * self.matrix[k][2]
                ])
        self.length[k] *= n
        self.compute_matrix_from_la()
        self.compute_atom_count()
        self.compute_fract_from_xyz()

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

    def afterparse_basic(self):
        """ Compute all the missing attributes of Crys,
        to be ready for any conversion to another format """
        self.check_parse()
        self.compute_atom_count()
        if self.inp_matrix:
            self.compute_la_from_matrix()
        elif self.inp_lengths_angles:
            self.compute_matrix_from_la()
        if self.inp_xyz:
            self.compute_fract_from_xyz()
        elif self.inp_fract:
            self.compute_xyz_from_fract()
        if not self.matrix_alligned:
            self.fix_cell_notalligned()

    def clear_cell_and_coord(self):
        self.inp_xyz = False
        self.inp_fract = False
        self.atom_xyz = []
        self.atom_fract = []
        self.inp_lengths_angles = False
        self.inp_matrix = False
        self.length = [0.0] * 3
        self.angle_deg = [0.0] * 3
        self.angle_rad = [0.0] * 3
        self.matrix = [[0.0] * 3 for i in range(3)]

    def print_bondvalence(self):

        print()
        print('Bond Valence Sum method: working only for Cu 1/2 as described in Shields2000')
        print('Read the code to see all the many caveats used to define bondings.')

        # read bvparm
        #bvparmf = open(os.path.dirname(__file__) + "/../data/bvparm2016.cif", 'r')
        bvparmf = open(os.path.dirname(__file__) + "/../data/bvparm2016_Cu_l.cif", 'r')
        lines = bvparmf.readlines()
        bvparmf.close()
        bvparm = []
        #for i in range(172, 2089):
        for i in range(13, 39):
            l = lines[i].split()
            bvparm.append([
                l[0],  #0_valence_param_atom_1
                int(l[1]),  #1_valence_param_atom_1_valence
                l[2],  #2_valence_param_atom_2
                int(l[3]),  #3_valence_param_atom_2_valence
                float(l[4]),  #4_valence_param_Ro
                float(l[5]),  #5_valence_param_B
                l[6],  #6_valence_param_ref_id
                lines[i].split("'")[1],  #7_valence_param_details
            ])

        ptab_ox = {'Cu': [1, 2]}
        ptab_organic = ['C', 'N', 'O', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'I']  #Note: I exclude all the others, incl. H

        # compute the BVsum for the first metal in the crys
        for i in range(self.natom):
            ielem = self.atom_element[i]

            if ielem == 'Cu':
                print('Atom: {} {}'.format(ielem, i))
                nneig = 0
                bond_elem_list = []
                bond_dist_list = []
                for j in range(self.natom):
                    jelem = self.atom_element[j]
                    if j != i and jelem in ptab_organic:
                        dist_thr = 0.8 * (ptab_vdw_uff[ielem] + ptab_vdw_uff[jelem])
                        dist = self.dist_ij(i, j)
                        if dist < dist_thr:
                            nneig += 1
                            # if N check further coordination
                            if jelem == 'N':
                                Nneig = 0
                                for k in range(self.natom):
                                    kelem = self.atom_element[j]
                                    if k != j and kelem in ptab_organic:
                                        kdist_thr = 0.5 * (ptab_vdw_uff[jelem] + ptab_vdw_uff[kelem])
                                        kdist = self.dist_ij(j, k)
                                        if kdist < kdist_thr:
                                            Nneig += 1
                                if Nneig == 0 or Nneig > 3:  # due to disorder or "one N atom" solvent
                                    Nneig = ''  # use general nitrogen coefficient
                                jelem += str(Nneig)
                            print('Neighbour {}: {}-{} dist= {:.3f} Angs'.format(nneig, ielem, jelem, dist))
                            bond_elem_list.append(jelem)
                            bond_dist_list.append(dist)
                print()
                ox_diff = {}
                for ox in ptab_ox[ielem]:
                    bvsum = 0

                    # find the appropriate parameters (stop at the first)
                    for j in range(len(bond_elem_list)):
                        jelem = bond_elem_list[j]
                        dist = bond_dist_list[j]
                        for p in bvparm:
                            if p[0] == ielem and p[1] == ox and p[2] == jelem:
                                R0 = p[4]
                                B = p[5]
                                break
                        else:
                            print('Most likely oxidation state: N/A (parameter missing!)')
                            sys.exit()
                        bv = math.exp((R0 - dist) / B)
                        print("Using:", *p, "(BV contrib: {:.2f} )".format(bv))
                        bvsum += bv
                    print('Hypothesys {}({}), BVSum= {:.2f}, diff= {:.3f}'.format(ielem, ox, bvsum, bvsum - ox))
                    ox_diff[ox] = abs(bvsum - ox)
                print('Most likely oxidation state: {}'.format(min(ox_diff, key=ox_diff.get)))
                sys.exit()
