"""Parsers for the different file formats.
Functions are sorted in alphabetic order.
"""

from manage_crystal import Crys
from manage_crystal.periodic_table import ptab_atnum_inv
from manage_crystal.utils import is_number
import numpy as np
import os
import sys
import re

ANGS2BOHR = 1.88973


def parsefloat(string):
    """ Robustly parse floats in files, such as:
    '1.234(5)' = 1.234, note that it means 1.234+/-0.005
    """
    string = float(re.sub(r'\([^)]*\)', '', string))
    return float(string)


def parse_from_filepath(filepath, tm=0):
    """ Utility that takes the filepath, recognise the file format type
    and return a Crys """
    if not os.path.isfile(filepath):
        sys.exit("ERROR: The file %s doesn't exist!" % filepath)
    inputformat = os.path.splitext(filepath)[1][1:]
    file = open(filepath, 'r')
    if inputformat in ["axsf", "xsf"]:
        crys = parse_axsf(file)
    elif inputformat == "car":
        crys = parse_car(file)
    elif inputformat == "cif":
        crys = parse_cif(file)
    elif inputformat in ["subsys", "inp", "restart"]:
        crys = parse_cp2k(file)
    elif inputformat == "cssr":
        crys = parse_cssr(file)
    elif inputformat == "cube":
        crys = parse_cube(file)
    elif inputformat == "pdb":
        crys = parse_pdb(file)
    elif inputformat == "POSCAR":
        crys = parse_poscar(file)
    elif inputformat in ["pwo", "pwi"]:
        crys = parse_pwo(file)
    elif inputformat == "xyz" and tm == 0:
        crys = parse_xyz(file)
    elif inputformat == "xyz" and tm == 3:  # B. Wells Qeq program
        crys = parse_xyz_tm3(file)
    else:
        sys.exit("WARNING: Input file format not implemented. EXIT.")

    file.close()
    crys.check_parse()
    crys.compute_atom_count()
    return crys


def parse_axsf(file):
    ''' Parse .axsf and .xsf files and return a Crys object '''
    c = Crys()
    while True:
        if file.readline().split()[0] == 'PRIMVEC':
            break
    for i in range(3):
        data = file.readline().split()
        for j in range(3):
            c.matrix[i][j] = float(data[j])
    while True:
        if file.readline().split()[0] == 'PRIMCOORD':
            break
    c.natom = int(file.readline().split()[0])
    for i in range(c.natom):
        data = file.readline().split()
        # In this format the atom type can be given as element or atomic number
        if is_number(data[0]):
            # convert from atomic number to element
            c.atom_type.append(ptab_atnum_inv[data[0]])
        else:
            c.atom_type.append(data[0])
        c.atom_xyz.append([float(data[1]), float(data[2]), float(data[3])])
    return c


def parse_car(file):
    ''' Parse .car file and return Crys object.
    TESTED just for Chongli Zhong's hCOF database.'''
    c = Crys()
    read_coord = False
    while True:
        line = file.readline()
        data = line.split()
        if line == "":
            break
        elif len(data) > 0 and data[0] == 'end':
            break
        elif len(data) > 0 and data[0] == 'PBC':
            c.length = [float(data[1]), float(data[2]), float(data[3])]
            c.angle_deg = [float(data[4]), float(data[5]), float(data[6])]
            read_coord = True
        elif read_coord:
            c.atom_type.append(data[0])
            c.atom_xyz.append([float(data[1]), float(data[2]), float(data[3])])
    return c


def parse_cif(file):
    ''' Parse .cif file and return a Crys object.
    Constraints:
    - only valid for P1 symmetry
    - cell data can be specified before/after the atom coordinates
    - it stops reading coordinates when it find a line with less than 4 columns (so no interruptions must exist)
    '''
    c = Crys()
    reading_coord = False
    while True:
        line = file.readline()
        data = line.split()
        if line == "":  #EOF
            break
        if len(data) > 0 \
           and len(data[0].split("_")) > 1 \
           and 'cell' in data[0].split("_"):
            reading_coord = False
            if data[0] == "_cell_length_a":
                c.length[0] = parsefloat(data[1])
            elif data[0] == "_cell_length_b":
                c.length[1] = parsefloat(data[1])
            elif data[0] == "_cell_length_c":
                c.length[2] = parsefloat(data[1])
            elif data[0] == "_cell_angle_alpha":
                c.angle_deg[0] = parsefloat(data[1])
            elif data[0] == "_cell_angle_beta":
                c.angle_deg[1] = parsefloat(data[1])
            elif data[0] == "_cell_angle_gamma":
                c.angle_deg[2] = parsefloat(data[1])
        # if the "_atom_site_***" section starts, remember the order
        if len(data) > 0 \
           and len(data[0].split("_")) > 1 \
           and data[0].split("_")[1] == "atom":
            data_order_dic = {}
            order = 0
            while len(data[0].split("_")) > 1 \
                  and data[0].split("_")[1] == "atom":
                data_order_dic[data[0]] = order
                line = file.readline()
                data = line.split()
                order += 1
            reading_coord = True  #entering in 'reading_coordinate MODE'

        # if less than 'Atom x y z' acan not be a coordinate!
        if len(data) < 4:
            reading_coord = False
        if reading_coord:
            # looks for "type_symbol" before and, if missing for "label"
            if "_atom_site_type_symbol" in data_order_dic:
                c.atom_type.append(data[data_order_dic["_atom_site_type_symbol"]])
            elif "_atom_site_label" in data_order_dic:
                c.atom_type.append(data[data_order_dic["_atom_site_label"]])
            else:
                sys.exit("EXIT: in cif missing type_symbol and label")
            if "_atom_site_fract_x" in data_order_dic:
                c.atom_fract.append([
                    float(data[data_order_dic["_atom_site_fract_x"]]),
                    float(data[data_order_dic["_atom_site_fract_y"]]),
                    float(data[data_order_dic["_atom_site_fract_z"]]),
                ])
            elif "_atom_site_Cartn_x" in data_order_dic:
                c.atom_xyz.append([
                    float(data[data_order_dic["_atom_site_Cartn_x"]]),
                    float(data[data_order_dic["_atom_site_Cartn_y"]]),
                    float(data[data_order_dic["_atom_site_Cartn_z"]]),
                ])
            else:
                sys.exit("EXIT: in cif missing fract_ and Cartn_ coordinates")
            if "_atom_site_charge" in data_order_dic:
                c.atom_charge.append(float(data[data_order_dic["_atom_site_charge"]]))
    return c


def parse_cp2k(file):
    ''' Parse any CP2K input file and return a Crys object.
    Constraints:
    - &CELL should be before &COORD
    - can read both A B C (cell) and ABC ALPHA_BETA_GAMMA (CELL)
    - it complains and exit if units are not [angstrom] and [deg]
    - if SCALED coord, the flag should be before the fract coordinates
    '''
    c = Crys()
    while True:
        data = file.readline().split()
        # Read cell: A B C
        cell_dict = {"A": 0, "B": 1, "C": 2}
        if len(data) > 0 and data[0] in cell_dict:
            if data[1][0] != "[":  # No unit specified. Default: Angstrom.
                shift = 0
            elif data[1].lower() == "[angstrom]":
                shift = 1
            else:
                sys.exit("WARNING: in parsing CP2K, weird units. EXIT")
            for i in range(3):
                c.matrix[cell_dict[data[0]]][i] = float(data[1 + i + shift])
        # Read CELL: ABC
        if len(data) > 3 and data[0] == 'ABC':
            if data[1][0] != "[":  # No unit specified. Default: Angstrom.
                shift = 0
            elif data[1].lower() == "[angstrom]":
                shift = 1
            else:
                sys.exit("WARNING: in parsing CP2K, weird units. EXIT")
            c.length[0] = float(data[1 + shift])
            c.length[1] = float(data[2 + shift])
            c.length[2] = float(data[3 + shift])
        # Read CELL: ALPHA_BETA_GAMMA
        if len(data) > 3 and data[0] == 'ALPHA_BETA_GAMMA':
            if data[1][0] != "[":  # No unit specified. Default: deg.
                shift = 0
            elif data[1].lower() == "[deg]":
                shift = 1
            else:
                sys.exit("WARNING: in parsing CP2K, weird units. EXIT")
            c.angle_deg[0] = float(data[1 + shift])
            c.angle_deg[1] = float(data[2 + shift])
            c.angle_deg[2] = float(data[3 + shift])
        if len(data) > 0 and (data[0] == "&COORD"):
            break
    scaled_coord = False  #Default
    while True:
        data = file.readline().split()
        if data[0] == "SCALED" \
         and (len(data)==1 or data[1].lower() in ["t", "true", ".true."]):
            scaled_coord = True
        elif data[0] == "SCALED" \
         and data[1].lower() in ["f", "false", ".false."]:
            scaled_coord = False
        elif data[0] == "&END":  # End of &COORD section
            break
        elif len(data) > 0:
            c.atom_type.append(data[0])
            if scaled_coord:
                c.atom_fract.append([float(data[1]), float(data[2]), float(data[3])])
            else:
                c.atom_xyz.append([float(data[1]), float(data[2]), float(data[3])])
    return c


def parse_cssr(file):
    ''' Parse .cssr file and return a Crys object '''
    # File format description: http://www.chem.cmu.edu/courses/09-560/docs/msi/modenv/D_Files.html#944777
    c = Crys()
    data = file.readline().split()
    c.length = [float(data[0]), float(data[1]), float(data[2])]
    data = file.readline().split()
    c.angle_deg = [float(data[0]), float(data[1]), float(data[2])]
    c.natom = int(file.readline().split()[0])
    junk = file.readline()
    for i in range(c.natom):
        data = file.readline().split()
        c.atom_type.append(data[1])
        c.atom_fract.append([float(data[2]), float(data[3]), float(data[4])])
        if len(data) == 14:
            c.atom_charge.append(float(data[13]))
    return c


def parse_cube(file):
    ''' Parse .cube file and return a Crys object '''
    c = Crys()
    junk = file.readline()  #header1
    junk = file.readline()  #header2
    data = file.readline().split()
    c.natom = int(data[0])
    for i in range(3):
        data = file.readline().split()
        for j in range(3):
            c.matrix[i][j] = float(data[0]) * float(data[j]) / ANGS2BOHR
    for i in range(c.natom):
        data = file.readline().split()
        # convert from atomic number to element
        c.atom_type.append(ptab_atnum_inv[int(data[0])])
        c.atom_xyz.append([float(data[2]) / ANGS2BOHR, float(data[3]) / ANGS2BOHR, float(data[4]) / ANGS2BOHR])
    return c


def parse_dcd_header(file):
    ''' Parse the dcd header '''
    data_dtype = np.dtype([('h01', 'i4', 1), ('h02', 'S4', 1), ('h03', 'i4', 9), ('h04', 'f4', 1), ('h05', 'i4', 10),
                           ('h06', 'i4', 1), ('h07', 'i4', 1), ('h08', 'i4', 1), ('h09', 'S80', 1), ('h10', 'S80', 1),
                           ('h11', 'i4', 1), ('h12', 'i4', 1), ('natoms', 'i4', 1), ('h13', 'i4', 1)])
    data = np.fromfile(file, data_dtype, 1)
    return data


def parse_dcd_snapshot(file, c):
    ''' Parse the dcd snapshot, updatyng the Crys cell and coordinates '''
    data_dtype = np.dtype([('junk1', 'i4', 1), ('len_ang', 'f8', 6), ('junk2', 'i4', 1), ('junk3', 'i4', 1),
                           ('coord_x', 'f4', c.natom), ('junk4', 'i4', 1),
                           ('junk5', 'i4', 1), ('coord_y', 'f4', c.natom), ('junk6', 'i4', 1), ('junk7', 'i4', 1),
                           ('coord_z', 'f4', c.natom), ('junk8', 'i4', 1)])
    data = np.fromfile(file, data_dtype, 1)
    if len(data) == 0:
        EOF = True
    else:
        EOF = False
        # Clear data to avoid conflicts
        c.clear_cell_and_coord()
        # Parsing the cell (be carefull to the order!)
        c.length[0] = data['len_ang'][0][0]
        c.length[1] = data['len_ang'][0][2]
        c.length[2] = data['len_ang'][0][5]
        c.angle_deg[0] = data['len_ang'][0][4]
        c.angle_deg[1] = data['len_ang'][0][3]
        c.angle_deg[2] = data['len_ang'][0][1]
        # Reset the coordinates and store the new ones
        c.atom_xyz = [[0.0] * 3 for i in range(c.natom)]
        for iatom in range(c.natom):
            c.atom_xyz[iatom][0] = data['coord_x'][0][iatom]
            c.atom_xyz[iatom][1] = data['coord_y'][0][iatom]
            c.atom_xyz[iatom][2] = data['coord_z'][0][iatom]
    return EOF


def parse_pdb(file):
    ''' Parse .pdb file and return Crys object '''
    c = Crys()
    while True:
        line = file.readline()
        data = line.split()
        if line == "":
            break
        elif len(data) > 0 and (data[0] == 'END' or data[0] == 'ENDMDL'):
            break
        elif len(data) > 0 and data[0] == 'CRYST1':
            c.length = [float(line[0o6:15]), float(line[15:24]), float(line[24:33])]
            c.angle_deg = [float(line[33:40]), float(line[40:47]), float(line[47:54])]
        elif len(data) > 0 and (data[0] == "ATOM" or data[0] == "HETATM"):
            c.atom_type.append(data[2])  #maybe read to data[-1]
            c.atom_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return c


def parse_poscar(file):
    ''' Parse Vasp's POSCAR file and return a Crys object '''
    c = Crys()
    junk_title = file.readline()
    scaling_factor = float(file.readline().split([0]))
    for i in range(3):
        data = file.readline().split()
        for j in range(3):
            c.matrix[i][j] = float(data[j]) * scaling_factor
    poscar_atomtypes = file.readline().split()
    poscar_atomnumbers = file.readline().split()
    for i in range(len(poscar_atomtypes)):
        for j in range(int(poscar_atomnumbers[i])):
            c.atom_type.append(poscar_atomtypes[i])
    c.natom = len(c.atom_type)
    coord_type = file.readline().split()[0]
    if coord_type.lower() == 'direct':
        for i in range(c.natom):
            data = file.readline().split()
            c.atom_fract.append([float(data[0]), float(data[1]), float(data[2])])
    elif coord_type.lower() == 'cartesian':
        for i in range(c.natom):
            data = file.readline().split()
            c.atom_xyz.append([float(data[0]), float(data[1]), float(data[2])])
    return c


def parse_pwo(file):
    ''' Parse Quantum Espresso's .pwo and .pwi files and return Crys object '''
    c = Crys()
    # Parse cell:
    # search for the last time the cell/coord are printed and jump to that
    # line (no need to be converged). ONLY if they are not found, it reads the
    # initial input cell
    with file as myFile:
        for num, line in enumerate(myFile, 1):
            if 'CELL_PARAMETERS' in line:
                cell_line = num
    file.seek(0)
    if 'cell_line' in locals():  #read cell in vc-relax calculation
        for i in range(0, cell_line):
            skip = file.readline()  #title line
        for i in range(3):
            data = file.readline().split()
            for j in range(3):
                c.matrix[i][j] = float(data[j])
    else:  #read cell in scf or relax calculation
        while True:
            data = file.readline().split()
            if len(data) > 0 and (data[0] == "celldm(1)="):
                celldm1 = float(data[1]) / ANGS2BOHR
                skip = file.readline().split()
                skip = file.readline().split()
                skip = file.readline().split()
                for i in range(3):
                    data = file.readline().split()
                    for j in range(3):
                        c.matrix[i][j] = float(data[3 + j]) * celldm1
                break
    file.seek(0)
    # Parse atomic coordinates
    with file as myFile:
        for num, line in enumerate(myFile, 1):
            if 'ATOMIC_POSITIONS' in line:
                atomicpositions_line = num
                if line.split()[1] == 'angstrom' \
                 or line.split()[1] == '(angstrom)':
                    readfractional = False
                elif line.split()[1] == 'crystal' \
                 or line.split()[1] == '(crystal)':
                    readfractional = True
    file.seek(0)
    if 'atomicpositions_line' in locals():  #read atomic in vc-relax and relax.
        for i in range(0, atomicpositions_line):
            skip = file.readline()
        i = 0
        while True:
            data = file.readline().split()
            if len(data) < 4:  #if the coordinates are finished, break
                break
            else:
                c.atom_type.append(data[0])
                if readfractional:
                    c.atom_fract.append([float(data[1]), float(data[2]), float(data[3])])
                else:
                    c.atom_xyz.append([float(data[1]), float(data[2]), float(data[3])])
    else:  #read atomic in scf calculation
        while True:
            data = file.readline().split()
            if len(data) > 0 and (data[0] == "celldm(1)="):
                celldm1 = float(data[1]) / ANGS2BOHR
            if len(data) > 3 and (data[3] == "positions"):
                while True:
                    data = file.readline().split()
                    if len(data) < 10:  #if the file is finished stop
                        break
                    else:
                        c.atom_type.append(data[1])
                        c.atom_xyz.append([float(data[6]), float(data[7]), float(data[8])])
                break
    return c


def parse_xyz(file):
    ''' Parse .xyz file and return Crys object '''
    c = Crys()
    c.natom = int(file.readline().split()[0])
    data = file.readline().split()
    if len(data) >= 7 and data[0] == 'CELL:':
        c.length = [float(data[1]), float(data[2]), float(data[3])]
        c.angle_deg = [float(data[4]), float(data[5]), float(data[6])]
    elif len(data) >= 10 and data[0] == 'cell:':
        c.matrix[0] = [float(data[1]), float(data[2]), float(data[3])]
        c.matrix[1] = [float(data[4]), float(data[5]), float(data[6])]
        c.matrix[2] = [float(data[7]), float(data[8]), float(data[9])]
    elif len(data) >= 23 and data[0] == 'jmolscript:':  #Chargemol
        c.matrix[0] = [float(data[10]), float(data[11]), float(data[12])]
        c.matrix[1] = [float(data[15]), float(data[16]), float(data[17])]
        c.matrix[2] = [float(data[20]), float(data[21]), float(data[22])]
    elif len(data) >= 9 and data[0].split("\"")[0] == "Lattice=":  #ASE
        c.matrix[0] = [float(data[0].split("\"")[1]), float(data[1]), float(data[2])]
        c.matrix[1] = [float(data[3]), float(data[4]), float(data[5])]
        c.matrix[2] = [float(data[6]), float(data[7]), float(data[8].split('\"')[0])]
    else:
        sys.exit('WARNING: xyz file with no cell provided! EXIT.')
    for i in range(c.natom):
        data = file.readline().split()
        c.atom_type.append(data[0])
        c.atom_xyz.append([float(data[1]), float(data[2]), float(data[3])])
        if len(data) == 5:  # if there is an extra column read it as charge
            c.atom_charge.append(float(data[4]))
    return c


def parse_xyz_tm3(file):
    ''' Parse .xyz file (tailor made #3) and return Crys object '''
    c = Crys()
    junk = file.readline()  #the first line has the filename
    data = file.readline().split()
    c.length = [float(data[1]), float(data[2]), float(data[3])]
    ac.angle_deg = [float(data[4]), float(data[5]), float(data[6])]
    c.natom = int(file.readline().split()[0])
    for i in range(c.natom):
        data = file.readline().split()
        c.atom_type.append(data[0])
        c.atom_fract.append([float(data[1]), float(data[2]), float(data[3])])
        if len(data) > 4:  # if not, the Qeq method crashed: no charge
            charge.append(float(data[5]))
    return c
