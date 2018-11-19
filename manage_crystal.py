#!/usr/bin/env python
"""
Python program to read coordinates from a file and handle them. Daniele Ongari 7/11/16

inputfilename          name of the input filename.inpformat
inputformat	       input format

natoms                 number of atoms
ntypes                 number of different atom types
ABC[012]               A B C in angstrom
abc[012]               alpha,beta,gamma in radians
cell(3x3)              unit cell matrix in angstrom

  atom[index]          atomic name
    an[index]          atomic number
   xyz[index][x,y,z]   cartesian coordinates in angstrom
 fract[index][x,y,z]   fractional coordinates
charge[index]          partial charge

atom_count[an]         number of atoms for atomic number an
atomic_symbol[an]      symbol of atom for atomic number an


"""
from __future__ import print_function  #python3 like print()
from __future__ import absolute_import
import string, sys
import numpy
import math
import subprocess
import argparse
from argparse import RawTextHelpFormatter  #needed to go next line in the help text
import os
import re  #re.split('(\d+)',"O23") = ['O', '23', '']
from six.moves import range

parser = argparse.ArgumentParser(
    description=
    'Program to read, extract info and convert crystal files (by Daniele Ongari)',
    formatter_class=RawTextHelpFormatter)

parser.add_argument(
    "inputfile",
    type=str,
    help="path to the input file to read\n" +
    "IMPLEMENTED: xyz(w/CELL),pdb,cssr,pwi,pwo,cif,xsf,axsf,subsys(CP2K),\n" +
    "             restart(CP2K),inp(CP2K),cube,POSCAR(VASP) \n" +
    "             [NEXT: gaussian, dcd+atoms,POSCAR(VASP)]")

parser.add_argument(
    "-o",
    "--output",
    action="store",
    type=str,
    dest="output",
    default=None,
    help="Output filename.extension or just the extension\n" +
    "IMPLEMENTED: cif,pdb,cssr,xyz(w/CELL),pwi,subsys(CP2K),axsf,geo(GULP)")

parser.add_argument(
    "-silent",
    action="store_true",
    dest="silent",
    default=False,
    help="No output info on the screen")

parser.add_argument(
    "-show",
    action="store_true",
    dest="show",
    default=False,
    help="Show all the info\n" + "[skip -silent]")

parser.add_argument(
    "-showonly",
    action="store",
    type=str,
    dest="showonly",
    default=None,
    help="Show only the required info:\n" + "cell, CELL, xyz, fract, charge\n"
    + "[skip -silent]")

parser.add_argument(
    "-cupw",
    action="store_true",
    dest="cupw",
    default=False,
    help="Look for a Copper PaddleWheel")

parser.add_argument(
    "-void",
    action="store_true",
    dest="void",
    default=False,
    help="Compute void geometrically [NOT WORKING]")

parser.add_argument(
    "-ovlp",
    action="store_true",
    dest="ovlp",
    default=False,
    help="Look for an overlap and modify the file [WORK IN PROGRESS]")

parser.add_argument(
    "-pseudopw",
    action="store",
    type=str,
    dest="pseudopw",
    default="pbe",
    help="Pseudo for the .pwi output")

parser.add_argument(
    "-bscp2k",
    action="store",
    type=str,
    dest="bscp2k",
    default="DZVP-MOLOPT-SR-GTH",
    help="Gaussian Basis Set for CP2K")

parser.add_argument(
    "-resp",
    action="store",
    type=str,
    dest="resp",
    default=None,
    help="Read the charges from a cp2k RESP file\n" +
    "(also checking if the atoms are the same)\n" +
    "BC1: it read the first set of charges\n" +
    "BC2: Also a cp2k output file with charges is fine!\n")

parser.add_argument(
    "-readcharge",
    action="store",
    type=str,
    dest="readcharge",
    default=None,
    help="Read the charges from a simple list")

parser.add_argument(
    "-readrepeatcharge",
    action="store",
    type=str,
    dest="readrepeatcharge",
    default=None,
    help="Read the charges from REPEAT output of QE")

parser.add_argument(
    "-x",
    action="store",
    type=int,
    dest="multipl_x",
    default=1,
    help="Extend in the x direction, by the specified times")

parser.add_argument(
    "-y",
    action="store",
    type=int,
    dest="multipl_y",
    default=1,
    help="Extend in the y direction, by the specified times")

parser.add_argument(
    "-z",
    action="store",
    type=int,
    dest="multipl_z",
    default=1,
    help="Extend in the z direction, by the specified times")

parser.add_argument(
    "-cutoff",
    action="store",
    type=float,
    dest="cutoff",
    default=None,
    help="Automatically extend the UC so that the cutoff is respected\n" +
    "(TIP: use -cutoff 0 to just know the perpendicular widths!)")

parser.add_argument(
    "-chargenull",
    action="store_true",
    dest="chargenull",
    default=False,
    help="Delete the charge of the atoms")

parser.add_argument(
    "-printatoms",
    action="store_true",
    dest="printatoms",
    default=False,
    help="Print all atoms types\n" + "[skip -silent]")

parser.add_argument(
    "-printatoms_noHCO",
    action="store_true",
    dest="printatoms_noHCO",
    default=False,
    help="Print all atoms types exc. H,C,O\n" + "[skip -silent]")

parser.add_argument(
    "-transl",
    action="store",
    type=float,
    nargs='*',
    dest="transl",
    default=None,
    help="x y z translation in Angs")

parser.add_argument(
    "-mol",
    action="store_true",
    dest="mol",
    default=False,
    help="Considers a molecule for xyz: no cell!"
)  #to ad later, now putting a 50x50x50 cell is fine!

parser.add_argument(
    "-randomize",
    action="store",
    type=float,
    dest="randomize",
    default=None,
    help="Randomize the geometry by a gaussian\n" +
    "with the specified delta (angs)")

################################################################################################### START stuff for Qeq project-22
parser.add_argument(
    "-chkmetalcharge",
    action="store_true",
    dest="chkmetalcharge",
    default=False,
    help="Check if the charge on a metal (see list) is neg.\n" +
    "[skip -silent]")

parser.add_argument(
    "-chkcharge",
    action="store_true",
    dest="chkcharge",
    default=False,
    help="Check if all the charges are zero.\n" + "[skip -silent]")

parser.add_argument(
    "-chkdef2",
    action="store_true",
    dest="chkdef2",
    default=False,
    help="Check if there is a non def2 BS atom (H-La, Hf-Rn).\n" +
    "[skip -silent]")

parser.add_argument(
    "-chkmepo",
    action="store_true",
    dest="chkmepo",
    default=False,
    help="Check if there is a non MEPO atom (H,V,Cu,Zn,C,N,O,F,Cl,Br,I).\n" +
    "[skip -silent]")

parser.add_argument(
    "-avgcharges",
    action="store_true",
    dest="avgcharges",
    default=False,
    help="Use average charges from DDEC")

parser.add_argument(
    "-normalizecharges",
    action="store_true",
    dest="normalizecharges",
    default=False,
    help="Normalize the charges to have a null total charge.")

################################################################################################### END stuff for Qeq project-22

parser.add_argument(
    "-tm1",
    action="store_true",
    dest="tailormade1",
    default=False,
    help="Tailor-made 1: read  .cif CoRE MOF w/DDEC charges")

parser.add_argument(
    "-tm2",
    action="store_true",
    dest="tailormade2",
    default=False,
    help="Tailor-made 2: print .cif for EQeq")

parser.add_argument(
    "-tm3",
    action="store_true",
    dest="tailormade3",
    default=False,
    help="Tailor-made 3: read  .xyz for B.Wells Qeq")

parser.add_argument(
    "-tm4",
    action="store_true",
    dest="tailormade4",
    default=False,
    help="Tailor-made 4: print .xyz for B.Wells Qeq")

parser.add_argument(
    "-tm5",
    action="store_true",
    dest="tailormade5",
    default=False,
    help="Tailor-made 5: print .xyz for B.Wells Qeq w/zero FC")

parser.add_argument(
    "-tm6",
    action="store_true",
    dest="tailormade6",
    default=False,
    help="Tailor-made 6: read GULP's cif")

args = parser.parse_args()

############################################################################# UTILITIES and STANDARD INFOS about the periodic table: atomic_symbol/name/vdw/mass
from atomic_data import *  #import all the data stored in the file atom_data.py

atom_count = [0] * 119  #anything assigned to 0, H_index=1, He_index=2, ...

ANGS2BOHR = 1.88973
AVOGCONST = 6.022E+23


def is_number(s):  #checks if a string is a number or not
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass


############################################################################## INPUT
if not args.silent: print()

#reading input file: name and format (notice that if there is a path it becomes part of the name, to have the output in the same place)
if not os.path.isfile(args.inputfile):
    sys.exit("ERROR: The file %s doesn't exist!" % args.inputfile)
inputfilename = os.path.splitext(args.inputfile)[0]
inputformat = os.path.splitext(args.inputfile)[1][1:]  
file = open(inputfilename + "." + inputformat, 'r')
"""/
if inputformat=='dcd':
	from pwtools import dcd
	cc,co = dcd.read_dcd_data('./'+inputfilename+'.'+inputformat, convang=False) #convang=True usefull to read LAMMPS data, because it converts angles from cosine to degree
	ABC=[ cc[-1][0], cc[-1][1], cc[-1][2] ]
	abc=[ math.acos(cc[-1][3]), math.acos(cc[-1][4]), math.acos(cc[-1][5]) ]
	cell=numpy.matrix([[                 ABC[0],                    0.0,                                                                           0.0],
			   [ABC[1]*math.cos(abc[2]),ABC[1]*math.sin(abc[2]),                                                                           0.0],
			   [ABC[2]*math.cos(abc[1]),ABC[2]*math.cos(abc[0]),math.sqrt(ABC[2]**2-(ABC[2]*math.cos(abc[1]))**2-(ABC[2]*math.cos(abc[0]))**2)]])
	xyz=co[-1]
	file.close()
	user_input = raw_input("Give me a file with the atoms: ")
	file = open(user_input,'r')
	natoms=int(file.readlines()[0].split()[0])
	atom=[]
	an=[]
	#xyz=[]
	charge=[]
	file = open(user_input,'r')
	for i in range(0,natoms):
		line = file.readlines()[i+2]
		data = line.split()
		atom.append(data[0])
		#an.append(atomic_symbol.index(atom[i]))
                #atom_count[an[i]]+=1
		an.append(1)
                atom_count[an[i]]+=1
		#xyz.append([float(data[1]), float(data[2]), float(data[3])])
		charge.append(0.000)       """

if inputformat == 'pdb':
    while True:
        line = file.readline()
        if line.split()[0] == 'CRYST1' or line.split()[0] == 'COMPD':
            break
    ABC = [float(line[0o6:15]), float(line[15:24]), float(line[24:33])]
    abc = [
        math.radians(float(line[33:40])),
        math.radians(float(line[40:47])),
        math.radians(float(line[47:54]))
    ]
    #read atom[index]
    atom = []
    an = []
    xyz = []
    i = 0
    while True:
        line = file.readline()
        data = line.split()
        if len(data) == 0:  #if the file is finished stop
            break
        elif data[0] == 'END' or data[0] == 'ENDMDL':  #if the file says "END"
            break
        elif data[0] != "ATOM" and data[0] != "HETATM":  #avoid other stuff
            donothing = True
        else:
            atom.append(data[2])
            an.append(atomic_symbol.index(atom[i]))
            atom_count[an[i]] += 1
            xyz.append(
                [float(line[30:38]),
                 float(line[38:46]),
                 float(line[46:54])])
            i = i + 1
        natoms = i

if inputformat == 'POSCAR':
    junk_title = file.readline()
    junk_symm = file.readline()
    a_vect = file.readline().split()
    b_vect = file.readline().split()
    c_vect = file.readline().split()
    cell = numpy.matrix(
        [[float(a_vect[0]),
          float(a_vect[1]),
          float(a_vect[2])],
         [float(b_vect[0]),
          float(b_vect[1]),
          float(b_vect[2])],
         [float(c_vect[0]),
          float(c_vect[1]),
          float(c_vect[2])]])
    poscar_atomtypes = file.readline().split()
    poscar_atomnumbers = file.readline().split()
    atom = []
    an = []
    for i in range(len(poscar_atomtypes)):
        for j in range(int(poscar_atomnumbers[i])):
            atom.append(poscar_atomtypes[i])
            an.append(atomic_symbol.index(poscar_atomtypes[i]))
            atom_count[atomic_symbol.index(poscar_atomtypes[i])] += 1
    natoms = len(atom)
    coord_type = file.readline().split()[0]
    if coord_type == 'Direct' or coord_type == 'direct':
        fract = []
        for i in range(natoms):
            coord = file.readline().split()
            fract.append([float(coord[0]), float(coord[1]), float(coord[2])])
    elif coord_type == 'Cartesian' or coord_type == 'cartesian':
        xyz = []
        for i in range(natoms):
            coord = file.readline().split()
            xyz.append([float(coord[0]), float(coord[1]), float(coord[2])])

if inputformat == 'cssr':  #file format description: http://www.chem.cmu.edu/courses/09-560/docs/msi/modenv/D_Files.html#944777
    celltemp = file.readline().split()
    ABC = [float(celltemp[0]), float(celltemp[1]), float(celltemp[2])]
    celltemp = file.readline().split()
    abc = [
        math.radians(float(celltemp[0])),
        math.radians(float(celltemp[1])),
        math.radians(float(celltemp[2]))
    ]
    data = file.readline().split()
    natoms = int(data[0])
    junk = file.readline()
    atom = []
    an = []
    fract = []
    charge = []
    for i in range(natoms):
        data = file.readline().split(
        )  #Atom serial number, atom name, x, y, z coordinates, bonding connectivities (max 8), charge
        atom.append(data[1])
        an.append(atomic_symbol.index(atom[i]))
        atom_count[an[i]] += 1
        fract.append([float(data[2]), float(data[3]), float(data[4])])
        if len(data) == 14:
            charge.append(float(data[13]))
        else:
            charge.append(0.0)

if inputformat == 'xyz':
    if not args.tailormade3:
        #read number of atoms
        line = file.readline()
        natoms = int(line.split()[0])

        #read cell in my way of writing it as a comment of xyz
        line = file.readline()

        if len(line) == 0 or (line.split()[0] != 'CELL:'
                              and line.split()[0] != 'cell:'
                              and line.split()[0] != 'jmolscript:'
                              ):  #set a 50x50x50 cell if CELL is not specified
            if not args.silent:
                print(len(line))
                print("WARNING: no CELL properly specified... using 50 50 50")
                ABC = [50., 50., 50.]
                abc = [math.radians(90.), math.radians(90.), math.radians(90.)]

        else:
            celltemp = line.split()

            if celltemp[0] == 'CELL:':  #ABC alpha/beta/gamma
                ABC = [
                    float(celltemp[1]),
                    float(celltemp[2]),
                    float(celltemp[3])
                ]
                abc = [
                    math.radians(float(celltemp[4])),
                    math.radians(float(celltemp[5])),
                    math.radians(float(celltemp[6]))
                ]

            elif celltemp[0] == 'cell:':  #cellmatrix
                cell = numpy.matrix([[
                    float(celltemp[1]),
                    float(celltemp[2]),
                    float(celltemp[3])
                ], [
                    float(celltemp[4]),
                    float(celltemp[5]),
                    float(celltemp[6])
                ], [
                    float(celltemp[7]),
                    float(celltemp[8]),
                    float(celltemp[9])
                ]])
            elif celltemp[0] == 'jmolscript:':  #working for ddec output
                cell = numpy.matrix([[
                    float(celltemp[10]),
                    float(celltemp[11]),
                    float(celltemp[12])
                ],
                                     [
                                         float(celltemp[15]),
                                         float(celltemp[16]),
                                         float(celltemp[17])
                                     ],
                                     [
                                         float(celltemp[20]),
                                         float(celltemp[21]),
                                         float(celltemp[22])
                                     ]])

        #read atom[index]
        atom = []
        an = []
        xyz = []
        charge = []
        for i in range(0, natoms):
            line = file.readline()
            data = line.split()
            atom.append(data[0])
            an.append(atomic_symbol.index(atom[i]))
            atom_count[an[i]] += 1
            xyz.append([float(data[1]), float(data[2]), float(data[3])])
            if len(data) == 5:
                charge.append(float(
                    data[4]))  #if there is an extra column is the charge!
        if len(charge) == 0: del charge

    if args.tailormade3:
        junk = file.readline()
        celltemp = file.readline().split()
        ABC = [float(celltemp[1]), float(celltemp[2]), float(celltemp[3])]
        abc = [
            math.radians(float(celltemp[4])),
            math.radians(float(celltemp[5])),
            math.radians(float(celltemp[6]))
        ]
        natoms = int(file.readline().split()[0])
        atom = []
        an = []
        fract = []
        charge = []
        for i in range(0, natoms):
            line = file.readline()
            data = line.split()
            atom.append(data[0])
            an.append(atomic_symbol.index(atom[i]))
            atom_count[an[i]] += 1
            fract.append([float(data[1]), float(data[2]), float(data[3])])
            if len(data) > 4:
                charge.append(float(data[5]))
            else:
                charge.append(float(0))

if (inputformat == 'pwo') or (inputformat == 'pwi'):
    # search for the last time the cell/coord are printed and jump to that
    # line (no need to be converged). ONLY if they are not found, it reads the
    # initial input cell
    with file as myFile:
        for num, line in enumerate(myFile, 1):
            if 'CELL_PARAMETERS' in line:
                cell_line = num
    file.close()

    file = open(args.inputfile, 'r')
    if 'cell_line' in locals():  #read cell in vc-relax calculation
        for i in range(0, cell_line):
            skip = file.readline()  #title line

        celltempA = file.readline().split()
        celltempB = file.readline().split()
        celltempC = file.readline().split()
        cell = numpy.matrix(
            [[float(celltempA[0]),
              float(celltempA[1]),
              float(celltempA[2])],
             [float(celltempB[0]),
              float(celltempB[1]),
              float(celltempB[2])],
             [float(celltempC[0]),
              float(celltempC[1]),
              float(celltempC[2])]])

    else:  #read cell in scf or relax calculation
        while True:
            data = file.readline().split()
            if len(data) > 0 and (data[0] == "celldm(1)="):
                celldm1 = float(data[1]) / ANGS2BOHR
                skip = file.readline().split()
                skip = file.readline().split()
                skip = file.readline().split()
                celltempA = file.readline().split()
                celltempB = file.readline().split()
                celltempC = file.readline().split()
                cell = numpy.matrix([[
                    float(celltempA[3]),
                    float(celltempA[4]),
                    float(celltempA[5])
                ],
                                     [
                                         float(celltempB[3]),
                                         float(celltempB[4]),
                                         float(celltempB[5])
                                     ],
                                     [
                                         float(celltempC[3]),
                                         float(celltempC[4]),
                                         float(celltempC[5])
                                     ]])
                cell *= celldm1
                break

    file.close()

    #atomic

    file = open(args.inputfile, 'r')
    with file as myFile:
        for num, line in enumerate(myFile, 1):
            if 'ATOMIC_POSITIONS' in line:
                atomic_line = num
                if line.split()[1] == 'angstrom' or line.split(
                )[1] == '(angstrom)':
                    readfractional = False
                elif line.split()[1] == 'crystal' or line.split(
                )[1] == '(crystal)':
                    readfractional = True
    file.close()

    file = open(args.inputfile, 'r')
    if 'atomic_line' in locals(
    ):  #read atomic in vc-relax and relax calculation
        atom = []
        an = []

        if readfractional: fract = []
        else: xyz = []

        for i in range(0, atomic_line):
            skip = file.readline()
        i = 0
        while True:
            data = file.readline().split()
            if len(data) < 4:  #if the file is finished stop
                break
            else:
                atom.append(
                    re.split('(\d+)', data[0])
                    [0])  #takes only the atomtype from a label like "Cu34"
                an.append(atomic_symbol.index(atom[i]))
                atom_count[an[i]] += 1
                if readfractional:
                    fract.append(
                        [float(data[1]),
                         float(data[2]),
                         float(data[3])])
                else:
                    xyz.append(
                        [float(data[1]),
                         float(data[2]),
                         float(data[3])])
                i = i + 1
        natoms = i

    else:  #read atomic in scf calculation
        atom = []
        an = []
        xyz = []
        while True:
            data = file.readline().split()
            if len(data) > 0 and (data[0] == "celldm(1)="):
                celldm1 = float(data[1]) / ANGS2BOHR
            if len(data) > 3 and (data[3] == "positions"):
                i = 0
                while True:
                    data = file.readline().split()
                    if len(data) < 10:  #if the file is finished stop
                        break
                    else:
                        atom.append(data[1])
                        an.append(atomic_symbol.index(atom[i]))
                        atom_count[an[i]] += 1
                        xyz.append([
                            float(data[6]) * celldm1,
                            float(data[7]) * celldm1,
                            float(data[8]) * celldm1
                        ])
                        i = i + 1
                natoms = i
                break

if inputformat == 'cif':
    if not (args.tailormade1 or args.tailormade6):
        if not args.silent:
            print("**** BE CAREFUL: cif reading is tested only for the format")
        if not args.silent:
            print("****             of this output (P1,label+symbol+fract)")
        if not args.silent:
            print("****             (also CoRE MOF cif are readable!)")
        if not args.silent:
            print("****             (ok for VESTA, AVOGADRO, ZEO++)")
        if not args.silent:
            print("**** Use obabel to have: label, type, x, y, z, occupancy")
        if not args.silent: print()
        ABC = [0] * 3
        abc = [0] * 3
        VESTA_CIF = False
        AVOGADRO_CIF = False
        ZEOPP_CIF = False
        while True:
            line = file.readline()
            data = line.split()
            if len(data) > 0 and (data[0] == "data_VESTA_phase_1"):
                VESTA_CIF = True
            if len(data) > 0 and (
                    line ==
                    "# CIF file generated by openbabel 2.3.2, see http://openbabel.sf.net\n"
            ):
                AVOGADRO_CIF = True
            if len(data) > 0 and (line == "# CIF file created by Zeo++\n"):
                ZEOPP_CIF = True
            if len(data) > 0 and (data[0] == "_cell_length_a"):
                ABC[0] = float(data[1])
            if len(data) > 0 and (data[0] == "_cell_length_b"):
                ABC[1] = float(data[1])
            if len(data) > 0 and (data[0] == "_cell_length_c"):
                ABC[2] = float(data[1])
            if len(data) > 0 and (data[0] == "_cell_angle_alpha"):
                abc[0] = math.radians(float(data[1]))
            if len(data) > 0 and (data[0] == "_cell_angle_beta"):
                abc[1] = math.radians(float(data[1]))
            if len(data) > 0 and (data[0] == "_cell_angle_gamma"):
                abc[2] = math.radians(float(data[1]))
            if len(data) > 0 and (data[0] == "_atom_site_fract_x"
                                  or data[0] == "_atom_site_Cartn_x"):
                atom = []
                an = []
                if (AVOGADRO_CIF): xyz = []
                else: fract = []
                charge = []
                i = 0
                while True:
                    data = file.readline().split()

                    if len(data) == 0: break  #end of file
                    if data[0] == "loop_": break  #flag loop_ (CoRE mof)

                    if (data[0][0] != "_"):
                        if (VESTA_CIF):
                            atom.append(re.split('(\d+)', data[7])[0])
                        elif (AVOGADRO_CIF):
                            atom.append(re.split('(\d+)', data[0])[0])
                        elif (ZEOPP_CIF):
                            atom.append(data[1])
                        else:
                            atom.append(re.split('(\d+)', data[0])[0])
                        an.append(atomic_symbol.index(atom[i]))
                        atom_count[an[i]] += 1
                        if (AVOGADRO_CIF):
                            xyz.append([
                                float(data[2]),
                                float(data[3]),
                                float(data[4])
                            ])
                        else:
                            fract.append([
                                float(data[2]),
                                float(data[3]),
                                float(data[4])
                            ])
                        if len(data) > 5 and is_number(data[5]):
                            charge.append(
                                float(data[5])
                            )  #BE CAREFUL: it could read other numbers!!!
                        else:
                            charge.append(0.0)
                        i += 1
                natoms = i
                break

    else:  #Tailormade1: DDEC CoRE MOF (only site_label and no site_type_symbol)
        if args.tailormade1:
            if not args.silent:
                print("**** READING .cif using tailor-made1 settings")
            if not args.silent:
                print("****         = CoRE MOF DDEC w/charges       ")
        if args.tailormade6:
            if not args.silent:
                print("**** READING .cif using tailor-made6 settings")
            if not args.silent:
                print("****         = GULP's cif w/charges          ")
            ABC = [0] * 3
            abc = [0] * 3

            while True:
                line = file.readline()
                data = line.split()
                if len(data) > 0 and (data[0] == "_cell_length_a"):
                    ABC[0] = float(data[1])
                if len(data) > 0 and (data[0] == "_cell_length_b"):
                    ABC[1] = float(data[1])
                if len(data) > 0 and (data[0] == "_cell_length_c"):
                    ABC[2] = float(data[1])
                if len(data) > 0 and (data[0] == "_cell_angle_alpha"):
                    abc[0] = math.radians(float(data[1]))
                if len(data) > 0 and (data[0] == "_cell_angle_beta"):
                    abc[1] = math.radians(float(data[1]))
                if len(data) > 0 and (data[0] == "_cell_angle_gamma"):
                    abc[2] = math.radians(float(data[1]))
                if len(data) > 0 and (data[0] == "_atom_site_fract_x"):
                    atom = []
                    an = []
                    fract = []
                    charge = []
                    i = 0
                    while True:
                        data = file.readline().split()
                        if len(data) == 0: break  #end of file

                        if (data[0][0] != "_"):
                            atom.append(
                                re.split('(\d+)', data[0])[0]
                            )  #takes only the atomtype from a label like "Cu34"
                            an.append(atomic_symbol.index(atom[i]))
                            atom_count[an[i]] += 1
                            if args.tailormade1:
                                fract.append([
                                    float(data[1]),
                                    float(data[2]),
                                    float(data[3])
                                ])
                            if args.tailormade6:
                                fract.append([
                                    float(data[3]),
                                    float(data[4]),
                                    float(data[5])
                                ])
                            if args.tailormade1: charge.append(float(data[4]))
                            if args.tailormade6: charge.append(float(data[6]))
                            i += 1
                    natoms = i
                    break

if inputformat == 'xsf' or inputformat == 'axsf':
    while True:
        line = file.readline()
        if line.split()[0] == 'PRIMVEC':
            break
    celltempA = file.readline().split()
    celltempB = file.readline().split()
    celltempC = file.readline().split()
    cell = numpy.matrix(
        [[float(celltempA[0]),
          float(celltempA[1]),
          float(celltempA[2])],
         [float(celltempB[0]),
          float(celltempB[1]),
          float(celltempB[2])],
         [float(celltempC[0]),
          float(celltempC[1]),
          float(celltempC[2])]])
    while True:
        line = file.readline()
        if line.split()[0] == 'PRIMCOORD':
            break
    line = file.readline()
    natoms = int(line.split()[0])

    atom = []
    an = []
    xyz = []
    for i in range(0, natoms):
        data = file.readline().split()
        if is_number(
                data[0]
        ):  #In .xsf is it specified the atom type number while in my .axsf the atom's name.
            an.append(int(data[0]))
            atom.append(atomic_symbol[an[i]])
        else:
            atom.append(data[0])
            an.append(atomic_symbol.index(atom[i]))
        atom_count[an[i]] += 1
        xyz.append([float(data[1]), float(data[2]), float(data[3])])

if inputformat == 'subsys' or inputformat == 'inp' or inputformat == 'restart':  #CP2K files
    while True:
        data = file.readline().split()
        if len(data) > 0 and (data[0] == "A"): celltempA = data
        if len(data) > 0 and (data[0] == "B"): celltempB = data
        if len(data) > 0 and (data[0] == "C"): celltempC = data
        if len(data) > 0 and (data[0] == "&COORD"): break

    if (celltempA[1] == '[angstrom]' or celltempA[1] == '[ANGSTROM]'):
        cell = numpy.matrix(
            [[float(celltempA[2]),
              float(celltempA[3]),
              float(celltempA[4])],
             [float(celltempB[2]),
              float(celltempB[3]),
              float(celltempB[4])],
             [float(celltempC[2]),
              float(celltempC[3]),
              float(celltempC[4])]])
    else:  #no units specified (default=angstrom)
        cell = numpy.matrix(
            [[float(celltempA[1]),
              float(celltempA[2]),
              float(celltempA[3])],
             [float(celltempB[1]),
              float(celltempB[2]),
              float(celltempB[3])],
             [float(celltempC[1]),
              float(celltempC[2]),
              float(celltempC[3])]])
    atom = []
    an = []
    xyz = []
    scaled_coord = False  #default (*SCALED*)
    i = 0
    while True:
        data = file.readline().split()
        if len(data) == 0: donothing = True
        elif data[0] == "SCALED" and (data[1] == 'T' or data[1] == 'TRUE'
                                      or data[1] == '.TRUE.'):
            scaled_coord = True  #Can be before or after the coordinates (*SCALED*)
        elif data[0] == "SCALED" and (data[1] == 'F' or data[1] == 'FALSE'
                                      or data[1] == '.FALSE.'):
            donothing = True  #default (*SCALED*)
        elif data[0] == "&END":
            natoms = i
            break
        else:
            atom.append(re.split(
                '(\d+)',
                data[0])[0])  #takes only the atomtype from a label like "Cu34"
            an.append(atomic_symbol.index(atom[i]))
            atom_count[an[i]] += 1
            xyz.append([float(data[1]),
                        float(data[2]),
                        float(data[3])
                        ])  #They will be scaled later if necessary (*SCALED*)
            i += 1
    if scaled_coord:  #Using scaled coordinates (*SCALED*)
        fract = xyz
        del xyz

if inputformat == 'cube':
    junk = file.readline()  #header1
    junk = file.readline()  #header2

    line = file.readline()

    natoms = int(line.split()[0])
    x_orig = float(line.split()[1])
    y_orig = float(line.split()[2])
    z_orig = float(line.split()[3])

    celltempA = file.readline().split()
    celltempB = file.readline().split()
    celltempC = file.readline().split()

    cell = numpy.matrix(
        [[
            float(celltempA[0]) * float(celltempA[1]) / ANGS2BOHR,
            float(celltempA[0]) * float(celltempA[2]) / ANGS2BOHR,
            float(celltempA[0]) * float(celltempA[3]) / ANGS2BOHR
        ],
         [
             float(celltempB[0]) * float(celltempB[1]) / ANGS2BOHR,
             float(celltempB[0]) * float(celltempB[2]) / ANGS2BOHR,
             float(celltempB[0]) * float(celltempB[3]) / ANGS2BOHR
         ],
         [
             float(celltempC[0]) * float(celltempC[1]) / ANGS2BOHR,
             float(celltempC[0]) * float(celltempC[2]) / ANGS2BOHR,
             float(celltempC[0]) * float(celltempC[3]) / ANGS2BOHR
         ]])
    atom = []
    an = []
    xyz = []
    for i in range(0, natoms):
        data = file.readline().split()
        an.append(int(data[0]))
        atom.append(atomic_symbol[an[i]])
        atom_count[an[i]] += 1
        xyz.append([
            float(data[2]) / ANGS2BOHR,
            float(data[3]) / ANGS2BOHR,
            float(data[4]) / ANGS2BOHR
        ])

file.close()

if not args.resp == None:
    if 'charge' in locals():
        if not args.silent:
            print(" ... THERE WERE ALREADY CHARGES BUT I'M OVERWRITING THEM!")
    charge = [0] * natoms
    with open(args.resp) as openfileobject:
        i = 0
        for line in openfileobject:
            data = line.split()
            if not len(data) < 4 and data[0] == 'RESP' and data[2] == atom[
                    i]:  #Removed condition {and int(data[1])==(i+1)} so that i can assemble the RESP file without changing the indexes !USE WITH CARE
                charge[i] = float(data[3])
                i = i + 1
            if i == natoms: break

if not args.readcharge == None:
    if 'charge' in locals():
        if not args.silent:
            print(" ... THERE WERE ALREADY CHARGES BUT I'M OVERWRITING THEM!")
    charge = [0] * natoms
    with open(args.readcharge) as openfileobject:
        i = 0
        for line in openfileobject:
            data = line.split()
            charge[i] = float(data[0])
            i = i + 1

if not args.readrepeatcharge == None:
    if 'charge' in locals():
        if not args.silent:
            print(" ... THERE WERE ALREADY CHARGES BUT I'M OVERWRITING THEM!")
    charge = [0] * natoms
    with open(args.readrepeatcharge) as openfileobject:
        if not args.silent:
            print("*** CHARGES from QE>REPEAT.out: multiplying by -0.5")
        countline = 0
        for line in openfileobject:
            data = line.split()
            if (countline - 17) >= 0 and (
                    countline - 17
            ) < natoms:  #Notice: there are 17 info lines in the header of REPEAT.out
                charge[countline - 17] = float(data[6]) * (-0.5)
            countline = countline + 1

    if not args.silent: print()
    if not args.silent:
        print(" ... %d atomic charges taken from %s" % (i, args.resp))
    if not i == natoms:
        print(
            "WARNING: the number of charges and atoms are different! Chech the order of your atoms!"
        )

############################################################################# DO SOMETHING
#check if xyz are really cartesian (angstrom) and if fract are really fractional coordinates.

if 'cell' in locals():  #make uc ABC+abc if it was read in cell
    if not args.silent: print()
    if not args.silent: print(" ... converting cell (matrix) to CELL (ABCabc)")
    #cell= numpy.tril(cell)
    ABC = [0] * 3
    abc = [0] * 3
    ABC[0] = math.sqrt(
        cell.item((0, 0)) * cell.item((0, 0)) + cell.item((0, 1)) * cell.item(
            (0, 1)) + cell.item((0, 2)) * cell.item((0, 2)))
    ABC[1] = math.sqrt(
        cell.item((1, 0)) * cell.item((1, 0)) + cell.item((1, 1)) * cell.item(
            (1, 1)) + cell.item((1, 2)) * cell.item((1, 2)))
    ABC[2] = math.sqrt(
        cell.item((2, 0)) * cell.item((2, 0)) + cell.item((2, 1)) * cell.item(
            (2, 1)) + cell.item((2, 2)) * cell.item((2, 2)))
    abc[0] = math.acos((cell.item((1, 0)) * cell.item((2, 0)) + cell.item(
        (1, 1)) * cell.item((2, 1)) + cell.item((1, 2)) * cell.item(
            (2, 2))) / ABC[1] / ABC[2])  #alpha=B^C
    abc[1] = math.acos((cell.item((0, 0)) * cell.item((2, 0)) + cell.item(
        (0, 1)) * cell.item((2, 1)) + cell.item((0, 2)) * cell.item(
            (2, 2))) / ABC[0] / ABC[2])  #beta=A^C
    abc[2] = math.acos((cell.item((0, 0)) * cell.item((1, 0)) + cell.item(
        (0, 1)) * cell.item((1, 1)) + cell.item((0, 2)) * cell.item(
            (1, 2))) / ABC[0] / ABC[1])  #gamma=A^B

elif 'ABC' in locals(
):  #make uc matrix if it was read in ABC+abc. Copied from Raspa>framework.c>UnitCellBox
    if not args.silent: print()
    if not args.silent:
        print(" ... converting CELL (ABCabc) to cell (matrix) ")
    tempd = (math.cos(abc[0]) -
             math.cos(abc[2]) * math.cos(abc[1])) / math.sin(abc[2])
    cell = numpy.matrix(
        [[ABC[0], 0.0, 0.0],
         [ABC[1] * math.cos(abc[2]), ABC[1] * math.sin(abc[2]), 0.0],
         [
             ABC[2] * math.cos(abc[1]), ABC[2] * tempd,
             ABC[2] * math.sqrt(1 - (math.cos(abc[1]))**2 - (tempd)**2)
         ]])

from numpy.linalg import inv
invcell = inv(cell)

if 'fract' in locals():  #convert in cartesian
    if not args.silent: print()
    if not args.silent:
        print(" ... converting fractional coordinates in cartesian")
    xyz = []
    for i in range(0, natoms):
        x = fract[i][0] * cell.item((0, 0)) + fract[i][1] * cell.item(
            (1, 0)) + fract[i][2] * cell.item((2, 0))
        y = fract[i][0] * cell.item((0, 1)) + fract[i][1] * cell.item(
            (1, 1)) + fract[i][2] * cell.item((2, 1))
        z = fract[i][0] * cell.item((0, 2)) + fract[i][1] * cell.item(
            (1, 2)) + fract[i][2] * cell.item((2, 2))
        xyz.append([x, y, z])
elif 'xyz' in locals():  #convert in fractionals
    if not args.silent: print()
    if not args.silent:
        print(" ... converting cartesian coordinates to fractional")
    fract = []
    for i in range(0, natoms):
        x = xyz[i][0] * invcell.item((0, 0)) + xyz[i][1] * invcell.item(
            (1, 0)) + xyz[i][2] * invcell.item((2, 0))
        y = xyz[i][0] * invcell.item((0, 1)) + xyz[i][1] * invcell.item(
            (1, 1)) + xyz[i][2] * invcell.item((2, 1))
        z = xyz[i][0] * invcell.item((0, 2)) + xyz[i][1] * invcell.item(
            (1, 2)) + xyz[i][2] * invcell.item((2, 2))
        fract.append([x, y, z])

if not 'charge' in locals():
    if not args.silent: print()
    if not args.silent:
        print(" ... no atomic charge found: 0 charge for each atom")
    charge = [0] * natoms

if args.chargenull:
    if not args.silent: print("*** chargenull: DELETING ALL THE CHARGES! ***")
    charge = [0] * natoms

if args.avgcharges:
    if not args.silent:
        print("*** avgcharges: Taking charges from atomic_ddecavgcharges ***")
    for i in range(0, natoms):
        charge[i] = atomic_ddecavgcharges[an[i]]

if args.normalizecharges:
    if not args.silent: print("")
    if not args.silent: print("*** NORMALIZING CHARGES ***")
    pos_charge = 0
    neg_charge = 0
    for i in range(0, natoms):
        if charge[i] > 0:
            pos_charge += charge[i]
        else:
            neg_charge += charge[i]

    tot_charge = pos_charge + neg_charge
    tot_abs = pos_charge - neg_charge
    pos_fract = pos_charge / tot_abs

    if not args.silent: print("total charge: %f" % tot_charge)
    if not args.silent: print("positive charges: %f" % pos_charge)
    if not args.silent: print("negative charges: %f" % neg_charge)
    if not args.silent: print("total absolute ch.: %f" % tot_abs)

    for i in range(0, natoms):
        if charge[i] > 0:
            charge[i] = charge[
                i] - tot_charge * pos_fract * charge[i] / pos_charge
        else:
            charge[i] = charge[i] - tot_charge * (
                1 - pos_fract) * charge[i] / neg_charge

############################################################################ APPLY TRANSLATION / RANDOMIZE
if args.transl != None:

    if not args.silent: print()
    if not args.silent:
        print("*** TRANSLATING coordinates by %f %f %f Angs" %
              (args.transl[0], args.transl[1], args.transl[2]))

    xyz_transl = []
    for i in range(0, natoms):
        x = xyz[i][0] + args.transl[0]
        y = xyz[i][1] + args.transl[1]
        z = xyz[i][2] + args.transl[2]
        xyz_transl.append([x, y, z])
    xyz = xyz_transl  #overwrite old xyz

    fract = []  #clear old fract
    for i in range(0, natoms):
        x = xyz[i][0] * invcell.item((0, 0)) + xyz[i][1] * invcell.item(
            (1, 0)) + xyz[i][2] * invcell.item((2, 0))
        y = xyz[i][1] * invcell.item((0, 1)) + xyz[i][1] * invcell.item(
            (1, 1)) + xyz[i][2] * invcell.item((2, 1))
        z = xyz[i][2] * invcell.item((0, 2)) + xyz[i][1] * invcell.item(
            (1, 2)) + xyz[i][2] * invcell.item((2, 2))
        fract.append([x, y, z])

if args.randomize != None:

    if not args.silent: print()
    if not args.silent:
        print(
            "*** RANDOMIZING XYZ coordinates by a normal distrib with delta=%f Angs"
            % args.randomize)

    xyz_random = []
    for i in range(0, natoms):
        x = xyz[i][0] + numpy.random.normal(0, args.randomize, 1)
        y = xyz[i][1] + numpy.random.normal(0, args.randomize, 1)
        z = xyz[i][2] + numpy.random.normal(0, args.randomize, 1)
        xyz_random.append([x, y, z])
    xyz = xyz_random  #overwrite old xyz

    fract = []  #clear old fract
    for i in range(0, natoms):
        x = float(
            xyz[i][0] * invcell.item((0, 0)) + xyz[i][1] * invcell.item(
                (1, 0)) + xyz[i][2] * invcell.item((2, 0))
        )  #no clear why float is needed here: otherwise it is considered a str
        y = float(xyz[i][1] * invcell.item((0, 1)) + xyz[i][1] * invcell.item(
            (1, 1)) + xyz[i][2] * invcell.item((2, 1)))
        z = float(xyz[i][2] * invcell.item((0, 2)) + xyz[i][1] * invcell.item(
            (1, 2)) + xyz[i][2] * invcell.item((2, 2)))
        fract.append([x, y, z])

############################################################################# CUTOFF TEST
if not args.cutoff == None:  #copied from raspa/framework.c/CellProperties(line:6184)

    if args.multipl_x > 1 or args.multipl_y > 1 or args.multipl_z > 1:
        sys.exit("Why did you ask for both -cutoff and -x -y -z ????")

    ax = cell.item((0, 0))
    ay = cell.item((0, 1))
    az = cell.item((0, 2))
    bx = cell.item((1, 0))
    by = cell.item((1, 1))
    bz = cell.item((1, 2))
    cx = cell.item((2, 0))
    cy = cell.item((2, 1))
    cz = cell.item((2, 2))

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
    VOLUME = math.fabs(ax * bxc1 + ay * bxc2 + az * bxc3)

    # calculate cell perpendicular widths
    perp_x = VOLUME / math.sqrt(bxc1**2 + bxc2**2 + bxc3**2)
    perp_y = VOLUME / math.sqrt(cxa1**2 + cxa2**2 + cxa3**2)
    perp_z = VOLUME / math.sqrt(axb1**2 + axb2**2 + axb3**2)

    print()
    print("CUTOFF_TEST | Cutoff: %.1f" % (args.cutoff))
    print("CUTOFF_TEST | Cell perpendicular widths: %.3f %.3f %.3f" %
          (perp_x, perp_y, perp_z))

    # compute how big the cell must be
    args.multipl_x = int(math.ceil(2 * args.cutoff / perp_x))
    args.multipl_y = int(math.ceil(2 * args.cutoff / perp_y))
    args.multipl_z = int(math.ceil(2 * args.cutoff / perp_z))

    if args.multipl_x > 1 or args.multipl_y > 1 or args.multipl_z > 1:
        print("CUTOFF_TEST | Expansion_done: %d %d %d for %s" %
              (args.multipl_x, args.multipl_y, args.multipl_z, args.inputfile))
    else:
        print("CUTOFF_TEST | Expansion_unnecesary: 1 1 1 for %s" %
              (args.inputfile))

############################################################################# CELL EXTENSION
if args.multipl_x > 1 or args.multipl_y > 1 or args.multipl_z > 1:

    if args.multipl_x > 1:
        for i in range(1, args.multipl_x):
            for j in range(0, natoms):
                atom.append(atom[j])
                an.append(atomic_symbol.index(atom[j]))
                atom_count[an[j]] += 1
                charge.append(charge[j])
                xyz.append([
                    xyz[j][0] + i * cell[0, 0], xyz[j][1] + i * cell[0, 1],
                    xyz[j][2] + i * cell[0, 2]
                ])
        ABC[0] *= args.multipl_x
        cell[0, 0] *= args.multipl_x  #nonzero
        cell[0, 1] *= args.multipl_x  #zero
        cell[0, 2] *= args.multipl_x  #zero
        natoms *= args.multipl_x

    if args.multipl_y > 1:
        for i in range(1, args.multipl_y):
            for j in range(0, natoms):
                atom.append(atom[j])
                an.append(atomic_symbol.index(atom[j]))
                atom_count[an[j]] += 1
                charge.append(charge[j])
                xyz.append([
                    xyz[j][0] + i * cell[1, 0], xyz[j][1] + i * cell[1, 1],
                    xyz[j][2] + i * cell[1, 2]
                ])
        ABC[1] *= args.multipl_y
        cell[1, 0] *= args.multipl_y  #nonzero
        cell[1, 1] *= args.multipl_y  #nonzero
        cell[1, 2] *= args.multipl_y  #zero
        natoms *= args.multipl_y

    if args.multipl_z > 1:
        for i in range(1, args.multipl_z):
            for j in range(0, natoms):
                atom.append(atom[j])
                an.append(atomic_symbol.index(atom[j]))
                atom_count[an[j]] += 1
                charge.append(charge[j])
                xyz.append([
                    xyz[j][0] + i * cell[2, 0], xyz[j][1] + i * cell[2, 1],
                    xyz[j][2] + i * cell[2, 2]
                ])
        ABC[2] *= args.multipl_z
        cell[2, 0] *= args.multipl_z  #nonzero
        cell[2, 1] *= args.multipl_z  #nonzero
        cell[2, 2] *= args.multipl_z  #nonzero
        natoms *= args.multipl_z

    fract = []
    invcell = inv(cell)
    for i in range(0, natoms):
        x = xyz[i][0] * invcell[0, 0] + xyz[i][1] * invcell[1, 0] + xyz[i][
            2] * invcell[2, 0]
        y = xyz[i][1] * invcell[0, 1] + xyz[i][1] * invcell[1, 1] + xyz[i][
            2] * invcell[2, 1]
        z = xyz[i][2] * invcell[0, 2] + xyz[i][1] * invcell[1, 2] + xyz[i][
            2] * invcell[2, 2]
        fract.append([x, y, z])

########################################################################################### COMPUTE INFO
# count atoms
if not args.silent: print()
ntypes = 0
for i in range(0, len(atom_count)):
    if atom_count[i] != 0:
        ntypes += 1
        if not args.silent:
            print('{0:>5} {1:3} atoms'.format(atom_count[i], atomic_symbol[i]))

if not args.silent: print(" ---- --- ----- ")
if not args.silent: print('{0:>5} {1:3} atoms'.format(natoms, 'tot'))

#compute volume (http://www.fxsolver.com/browse/formulas/Triclinic+crystal+system+(Unit+cell's+volume))
volume = ABC[0] * ABC[1] * ABC[2] * math.sqrt(
    1 - (math.cos(abc[0]))**2 - (math.cos(abc[1]))**2 - (math.cos(abc[2]))**2 +
    2 * math.cos(abc[0]) * math.cos(abc[1]) * math.cos(abc[2]))
if not args.silent: print()
if not args.silent: print("Volume: %.3f (Angtrom^3/u.c.)" % volume)

#compute density
weight = 0
for i in range(0, len(atom_count)):
    if atom_count[i] != 0:
        weight += atom_count[i] * atomic_mass[i]  #g/mol_uc
rho = weight / volume / AVOGCONST * 1E+10**3 / 1000  #Kg/m3
if not args.silent:
    print(
        "Density: %.5f (kg/m3), %.5f (g/cm3), %.5f (g/molUC)" %
        (rho, rho / 1000, weight), )

#compute conversion to mol/kg
molkg = 1000 / weight  #mol/g
if not args.silent:
    print("Conversion: 1 molec./u.c. = %.5f (mol/kg)" % (molkg))

#check the NET_CHARGE
if not args.silent: print()
#if not args.silent: print("Net charge: %.10f" %sum(charge))
if sum(charge) == 0 and max(charge) < 0.001:
    if not args.silent: print("NET_CHARGE: all the charges are zero.")
elif sum(charge) > -0.001 and sum(charge) < +0.001:
    if not args.silent: print("NET_CHARGE: negligible (|sum|<0.001).")
    #charge[0]=charge[0]-sum(charge)
    #print("*** Negligible error due to the rounding: subtracted from the first atom!")
    #print("*** Now the net charge is %.10f" %sum(charge))
else:
    if not args.silent:
        print("NET_CHARGE: nonzero (%.3f). ***WARNING***" % sum(charge))

#check negative charge on metals [skip -silent]
if args.chkmetalcharge:
    metal_list = [
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Hf",
        "Ta", "W"
    ]
    found_met = False
    found_met_neg = False
    found_met_nonzero = False
    found_met_notnumber = False
    for i in range(0, natoms):
        if (atom[i] in metal_list):
            found_met = True
            if math.isnan(charge[i]):
                print("CHK_METAL_CHARGE: not_number >>> %s=%s" % (atom[i],
                                                                  charge[i]))
                found_met_notnumber = True
            elif charge[i] < 0:
                print("CHK_METAL_CHARGE: found_negative >>> %s=%.3f" %
                      (atom[i], charge[i]))
                found_met_neg = True
            elif charge[i] > 0:
                found_met_nonzero = True
        if found_met_notnumber or found_met_neg: break
    if not found_met: print("CHK_METAL_CHARGE: no_metals")
    if found_met and not found_met_neg and not found_met_notnumber and not found_met_nonzero:
        print("CHK_METAL_CHARGE: all_zero")
    if found_met and not found_met_neg and not found_met_notnumber and found_met_nonzero:
        print("CHK_METAL_CHARGE: ok_positive")

#check if the charges are assigned / all zero / nan [skip -silent]
if args.chkcharge:
    found_nonzero = False
    found_notnumber = False
    found_weird = False
    for i in range(0, natoms):
        if math.isnan(charge[i]):
            print("CHK_CHARGE: not_number >>> %s=%s" % (atom[i], charge[i]))
            found_notnumber = True
            break
        elif charge[i] > 3 or charge[i] < -3:
            found_weird = True
            print("CHK_CHARGE: weird_charges")
            break
        elif charge[i] != 0:
            found_nonzero = True
            print("CHK_CHARGE: charged_framework")
            break
    if not found_notnumber and not found_weird and not found_nonzero:
        print("CHK_CHARGE: all_zero")

#check non-def2 atoms [skip -silent]
if args.chkdef2:
    found_nondef2 = False
    for i in range(0, natoms):
        if 58 <= an[i] <= 71 or an[i] >= 87:
            found_nondef2 = True
            print("CHK_def2: found %s" % atom[i])
            break
    if not found_nondef2: print("CHK_def2: ok")

#check non-mepo atoms [skip -silent]
if args.chkmepo:
    mepo_list = ["H", "V", "Cu", "Zn", "C", "N", "O", "F", "Cl", "Br", "I"]
    found_nonmepo = False
    for i in range(0, natoms):
        if not (atom[i] in mepo_list):
            found_nonmepo = True
            print("CHK_mepo: found %s" % atom[i])
            break
    if not found_nonmepo: print("CHK_mepo: ok")

#number of electrons
nelectrons = 0
for i in range(0, len(atom_count)):
    if atom_count[i] != 0:
        nelectrons += atom_count[i] * i  #nuber_of_atoms_with_AN*AN
if not args.silent: print("Tot. electrons: %d" % nelectrons)

#print atoms on one line for info
if args.printatoms:
    for i in range(0, len(atom_count)):
        if atom_count[i] != 0:
            print(atomic_symbol[i], end='_')
    print("")
if args.printatoms_noHCO:
    for i in range(0, len(atom_count)):
        if atom_count[i] != 0 and atomic_symbol[i]!="H" \
                                     and atomic_symbol[i]!="C" \
                                     and atomic_symbol[i]!="O":
            print(atomic_symbol[i], end='_')
    print("")
################################################################################################# OUTPUT OPERATIONS
if args.output == None:  #CHECK IF AN OUTPUT IS DEFINED
    outputfile = 'NOTHING'
else:
    if len(args.output.split(".")) > 1:  #output defined as name.format
        outputfilename = os.path.splitext(args.output)[0]
        outputformat = os.path.splitext(args.output)[1][1:]
        outputfile = outputfilename + "." + outputformat
    else:  #output defined as format
        outputfilename = inputfilename
        outputformat = args.output
        outputfile = outputfilename + "." + outputformat

if not args.silent: print()
if not args.silent:
    print("***************************************************")
if not args.silent:
    print("  Converting %s to %s" % (inputfilename + "." + inputformat,
                                     outputfile))
if not args.silent:
    print("***************************************************")
if not args.silent: print()

#show and showonly
if args.show:
    print(
        "cell ---------------------------------------------------------------")
if args.show or args.showonly == "cell":
    print("     %10.5f %10.5f %10.5f" % (cell.item((0, 0)), cell.item(
        (0, 1)), cell.item((0, 2))))
if args.show or args.showonly == "cell":
    print("     %10.5f %10.5f %10.5f" % (cell.item((1, 0)), cell.item(
        (1, 1)), cell.item((1, 2))))
if args.show or args.showonly == "cell":
    print("     %10.5f %10.5f %10.5f" % (cell.item((2, 0)), cell.item(
        (2, 1)), cell.item((2, 2))))
if args.show:
    print(
        "CELL (ABC, abc) ----------------------------------------------------")
if args.show or args.showonly == "CELL":
    print(" %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  " %
          (ABC[0], ABC[1], ABC[2], math.degrees(abc[0]), math.degrees(abc[1]),
           math.degrees(abc[2])))
if args.show:
    print(
        "xyz ----------------------------------------------------------------")
if args.show or args.showonly == "xyz":
    for i in range(0, natoms):
        print("%3s %10.5f %10.5f %10.5f " % (atom[i], xyz[i][0], xyz[i][1],
                                             xyz[i][2]))
if args.show:
    print(
        "fract --------------------------------------------------------------")
if args.show or args.showonly == "fract":
    for i in range(0, natoms):
        print("%3s %10.5f %10.5f %10.5f " % (atom[i], fract[i][0], fract[i][1],
                                             fract[i][2]))
if args.show:
    print(
        "charges --------------------------------------------------------------"
    )
if args.show or args.showonly == "charge":
    for i in range(0, natoms):
        print("%3s %10.5f " % (atom[i], charge[i]))
if args.show or args.showonly != None: sys.exit()

#computing void matematically: not working because of three spheres overlapping
if args.void:
    volsphere = 0
    volumeocc = 0
    d = [0] * 3
    s = [0] * 3
    t = [0] * 3
    for i in range(0, natoms):
        Ri = atomic_vdw_UFF[an[i]]
        volsphere += (4. / 3.) * math.pi * (Ri**3)
        volumeocc += (4. / 3.) * math.pi * (Ri**3)
        for j in range((i + 1), natoms):
            Rj = atomic_vdw_UFF[an[j]]

            d[0] = xyz[i][0] - xyz[j][0]
            d[1] = xyz[i][1] - xyz[j][1]
            d[2] = xyz[i][2] - xyz[j][2]

            s[0] = invcell.item((0, 0)) * d[0] + invcell.item(
                (1, 0)) * d[1] + invcell.item((2, 0)) * d[2]
            s[1] = invcell.item((0, 1)) * d[0] + invcell.item(
                (1, 1)) * d[1] + invcell.item((2, 1)) * d[2]
            s[2] = invcell.item((0, 2)) * d[0] + invcell.item(
                (1, 2)) * d[1] + invcell.item((2, 2)) * d[2]

            t[0] = s[0] - int(round(s[0]))
            t[1] = s[1] - int(round(s[1]))
            t[2] = s[2] - int(round(s[2]))

            d[0] = cell.item((0, 0)) * t[0] + cell.item(
                (1, 0)) * t[1] + cell.item((2, 0)) * t[2]
            d[1] = cell.item((0, 1)) * t[0] + cell.item(
                (1, 1)) * t[1] + cell.item((2, 1)) * t[2]
            d[2] = cell.item((0, 2)) * t[0] + cell.item(
                (1, 2)) * t[1] + cell.item((2, 2)) * t[2]

            mindist = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2)

            if (Ri + Rj < mindist): V_ovlp = 0
            else:
                V_ovlp = math.pi * (Ri + Rj - mindist)**2 * (
                    mindist**2 + 2 * mindist * Rj - 3 * Rj**2 + 2 * mindist *
                    Ri + 6 * Rj * Ri - 3 * Ri**2) / (12. * mindist)
            volumeocc -= V_ovlp  #works only if there are not three spheres with a common overlapping

#print("Volume without atom spheres:  %.3f (Angtrom^3/u.c.)" %(volume-volumeocc))
    print("Void fraction (cons. ovlp):       %.6f [= %.3f non-void (A^3)]" %
          (1 - volumeocc / volume, volumeocc))
    print("Void fraction (negl. ovlp):       %.6f [= %.3f non-void (A^3)]" %
          (1 - volsphere / volume, volsphere))
    print()
    #sys.exit("YOU JUST ASKED for VOID: no external files printed!")
    sys.exit()

#This function checks if there are copper paddlewheels (= a copper atom with a close Cu and 4 close O)
if args.cupw:
    ncupw_act = 0
    ncupw_sol = 0
    ncupw_wrd = 0
    for i in range(0, natoms):
        closeCu = 0
        closeO = 0
        closeX = 0

        if (an[i] == 29):
            for j in range(0, natoms):
                if (j != i):
                    d = [0] * 3
                    s = [0] * 3
                    t = [0] * 3

                    #d[0]=xyz[i][0]-xyz[j][0]
                    #d[1]=xyz[i][1]-xyz[j][1]
                    #d[2]=xyz[i][2]-xyz[j][2]
                    #s[0]=invcell.item((0,0))*d[0]+invcell.item((1,0))*d[1]+invcell.item((2,0))*d[2]
                    #s[1]=invcell.item((0,1))*d[0]+invcell.item((1,1))*d[1]+invcell.item((2,1))*d[2]
                    #s[2]=invcell.item((0,2))*d[0]+invcell.item((1,2))*d[1]+invcell.item((2,2))*d[2]

                    s[0] = fract[i][0] - fract[j][0]
                    s[1] = fract[i][1] - fract[j][1]
                    s[2] = fract[i][2] - fract[j][2]

                    t[0] = s[0] - int(round(s[0]))
                    t[1] = s[1] - int(round(s[1]))
                    t[2] = s[2] - int(round(s[2]))

                    d[0] = cell.item((0, 0)) * t[0] + cell.item(
                        (1, 0)) * t[1] + cell.item((2, 0)) * t[2]
                    d[1] = cell.item((0, 1)) * t[0] + cell.item(
                        (1, 1)) * t[1] + cell.item((2, 1)) * t[2]
                    d[2] = cell.item((0, 2)) * t[0] + cell.item(
                        (1, 2)) * t[1] + cell.item((2, 2)) * t[2]

                    mindist = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2)

                    if (an[j] == 29) and (1.8 < mindist < 2.8):
                        closeCu += 1
                        #print("%d %d %f     %d %d" %(i,j,mindist,closeO,closeX)
                    if (an[j] == 8) and (1.5 < mindist < 2.5): closeO += 1
                    if (an[j] != 29) and (an[j] != 8) and (1.8 < mindist <
                                                           2.5):
                        closeX += 1
                        #print("%d %d %s  %f------%f %f %f " %(i,j,an[j],mindist,s[0],s[1],s[2]); print(fract[i]; print(fract[j]

            if (closeCu == 1) and (closeO == 4):
                ncupw_act += 1
                #elif (closeCu==1) and (closeO==4) and (closeX>0):  ncupw_sol+=1
            elif (closeCu == 1) and (closeO == 5):
                ncupw_sol += 1
            elif (closeCu < 0):
                ncupw_wrd += 1

    print("Cu-paddlewheel activated found:    %d" % ncupw_act)
    print(
        "Cu-paddlewheel solvated  found:    %d   (only Oxygen atoms considered)"
        % ncupw_sol)
    print("Cu-paddlewheel WEIRD     found:    %d" % ncupw_wrd)
    print()

    #print that I found it
    if True:
        ofile = open("000_cupw_found.txt", 'a')
        if (ncupw_act > 0) or (ncupw_act > 0) or (ncupw_wrd < 0):
            print(args.inputfile, file=ofile)
            print(
                "Cu-paddlewheel activated found:    %d" % ncupw_act,
                file=ofile)
            print(
                "Cu-paddlewheel solvated  found:    %d   (only Oxygen atoms considered)"
                % ncupw_sol,
                file=ofile)
            print(
                "Cu-paddlewheel WEIRD     found:    %d" % ncupw_wrd,
                file=ofile)
            print(" ", file=ofile)

    sys.exit()

#This function checks if two atoms overlap because of a bad PBC wrap
if args.ovlp:
    jlist = []
    for i in range(0, natoms):
        for j in range(i + 1, natoms):
            d = [0] * 3
            s = [0] * 3
            t = [0] * 3

            s[0] = fract[i][0] - fract[j][0]
            s[1] = fract[i][1] - fract[j][1]
            s[2] = fract[i][2] - fract[j][2]

            t[0] = s[0] - int(round(s[0]))
            t[1] = s[1] - int(round(s[1]))
            t[2] = s[2] - int(round(s[2]))

            d[0] = cell.item((0, 0)) * t[0] + cell.item(
                (1, 0)) * t[1] + cell.item((2, 0)) * t[2]
            d[1] = cell.item((0, 1)) * t[0] + cell.item(
                (1, 1)) * t[1] + cell.item((2, 1)) * t[2]
            d[2] = cell.item((0, 2)) * t[0] + cell.item(
                (1, 2)) * t[1] + cell.item((2, 2)) * t[2]

            mindist = math.sqrt(d[0]**2 + d[1]**2 + d[2]**2)

            if (mindist < 0.2):
                print("Overlap found between:")
                print("%3s %9.5f %9.5f %9.5f " % (atom[i], xyz[i][0],
                                                  xyz[i][1], xyz[i][2]))
                print("%3s %9.5f %9.5f %9.5f " % (atom[j], xyz[j][0],
                                                  xyz[j][1], xyz[j][2]))
                if atom[i] == atom[
                        j]:  #the two atoms are the same and they are overlapping
                    jlist.append(j)
                else:  #something weird is happening, two different atoms are overlapping
                    print(
                        "!!!!!! CHECK THE CRYSTAL, DIFFERENT ATOMS OVERLAPPING !!!!!!"
                    )
                    sys.exit(
                        "!!!!!! CHECK THE CRYSTAL, DIFFERENT ATOMS OVERLAPPING !!!!!!"
                    )

    if len(jlist) > 0:  #correct overlaps
        natoms = natoms - len(jlist)
        atom = [i for j, i in enumerate(atom) if j not in jlist]
        an = [i for j, i in enumerate(an) if j not in jlist]
        xyz = [i for j, i in enumerate(xyz) if j not in jlist]
        fract = [i for j, i in enumerate(fract) if j not in jlist]
        charge = [i for j, i in enumerate(charge) if j not in jlist]
        print("OVERLAPS FOUND: %d" % len(jlist))
        #continue and overwrite the file
        outputfile = args.inputfile
        outputformat = args.inputfile.split(".")[-1]
    else:
        print("OVERLAPS FOUND: %d" % len(jlist))
        sys.exit()

############################################################################## Write converted files
if args.output != None:

    ofile = open(outputfile, 'w+')

    if outputformat == "cif":
        if not args.tailormade2:
            print("data_crystal", file=ofile)
            print(" ", file=ofile)
            print("_cell_length_a    %.5f" % ABC[0], file=ofile)
            print("_cell_length_b    %.5f" % ABC[1], file=ofile)
            print("_cell_length_c    %.5f" % ABC[2], file=ofile)
            print("_cell_angle_alpha %.5f" % math.degrees(abc[0]), file=ofile)
            print("_cell_angle_beta  %.5f" % math.degrees(abc[1]), file=ofile)
            print("_cell_angle_gamma %.5f" % math.degrees(abc[2]), file=ofile)
            print(" ", file=ofile)
            print("_symmetry_space_group_name_Hall 'P 1'", file=ofile)
            print("_symmetry_space_group_name_H-M  'P 1'", file=ofile)
            print(" ", file=ofile)
            print("loop_", file=ofile)
            print("_symmetry_equiv_pos_as_xyz", file=ofile)
            print(" 'x,y,z' ", file=ofile)
            print(" ", file=ofile)
            print("loop_", file=ofile)
            print("_atom_site_label", file=ofile)
            print("_atom_site_type_symbol", file=ofile)
            print("_atom_site_fract_x", file=ofile)
            print("_atom_site_fract_y", file=ofile)
            print("_atom_site_fract_z", file=ofile)
            print("_atom_site_charge", file=ofile)
            for i in range(0, natoms):
                label = atom[
                    i]  #removed: label=atom[i]+"_"+str(i+1) because the number makes Raspa extremely verbose
                print(
                    ('{0:10} {1:5} {2:>9.5f} {3:>9.5f} {4:>9.5f} {5:>14.10f}'.
                     format(label, atom[i], fract[i][0], fract[i][1],
                            fract[i][2], charge[i])),
                    file=ofile)
        #Better to keep a lot of decimals to avoid small net charged

        if args.tailormade2:

            print("****PRINTING .CIF TAILOR-MADE2 FOR EQeq***", file=ofile)

            print("data_crystal", file=ofile)
            print("loop_", file=ofile)
            print("_symmetry_equiv_pos_as_xyz", file=ofile)
            print(" 'x,y,z' ", file=ofile)
            print("loop_", file=ofile)
            print("_cell_length_a    %.5f" % ABC[0], file=ofile)
            print("_cell_length_b    %.5f" % ABC[1], file=ofile)
            print("_cell_length_c    %.5f" % ABC[2], file=ofile)
            print("_cell_angle_alpha %.5f" % math.degrees(abc[0]), file=ofile)
            print("_cell_angle_beta  %.5f" % math.degrees(abc[1]), file=ofile)
            print("_cell_angle_gamma %.5f" % math.degrees(abc[2]), file=ofile)
            print("_symmetry_space_group_name_Hall 'P 1'", file=ofile)
            print("_symmetry_space_group_name_H-M  'P 1'", file=ofile)
            print("_atom_site_label", file=ofile)
            print("_atom_site_type_symbol", file=ofile)
            print("_atom_site_fract_x", file=ofile)
            print("_atom_site_fract_y", file=ofile)
            print("_atom_site_fract_z", file=ofile)
            for i in range(0, natoms):
                label = atom[i]  #removed: label=atom[i]+"_"+str(i+1)
                print(('{0:10} {1:5} {2:>9.5f} {3:>9.5f} {4:>9.5f}'.format(
                    label, atom[i], fract[i][0], fract[i][1], fract[i][2])),
                      file=ofile)
            print("_loop", file=ofile)

    if outputformat == "pdb":
        print((
            'CRYST1{0:>9.3f}{1:>9.3f}{2:>9.3f}{3:>7.2f}{4:>7.2f}{5:>7.2f} P 1           1'
            .format(ABC[0], ABC[1], ABC[2], math.degrees(abc[0]),
                    math.degrees(abc[1]), math.degrees(abc[2]))),
              file=ofile)
        for i in range(0, natoms):
            print(
                "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.5f%8.5f%8.5f%6.2f%6.2f          %2s%2s"
                % ("ATOM", i + 1, atom[i], "", "XXX", "X", 1, "", fract[i][0],
                   fract[i][1], fract[i][2], 1.00, 0.00, atom[i], ""),
                file=ofile)

    if outputformat == "cssr":
        print(
            "                               %.3f  %.3f  %.3f" %
            (ABC[0], ABC[1], ABC[2]),
            file=ofile)
        print(
            "                %.3f   %.3f   %.3f   SPGR =  1 P 1         OPT = 1"
            % (math.degrees(abc[0]), math.degrees(abc[1]), math.degrees(
                abc[2])),
            file=ofile)
        print("%d   0" % (natoms), file=ofile)
        print("0 %s       : %s" % (inputfilename, inputfilename), file=ofile)
        for i in range(0, natoms):
            print(
                "%4d %3s %8.5f %8.5f %8.5f    0  0  0  0  0  0  0  0  %7.5f" %
                (i + 1, atom[i], fract[i][0], fract[i][1], fract[i][2],
                 charge[i]),
                file=ofile)

    if outputformat == "xyz":
        if not (args.tailormade4 or args.tailormade5):
            print("%d" % (natoms), file=ofile)
            print(
                "CELL:  %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
                (ABC[0], ABC[1], ABC[2], math.degrees(abc[0]),
                 math.degrees(abc[1]), math.degrees(abc[2])),
                file=ofile)
            for i in range(0, natoms):
                print(
                    "%3s %9.5f %9.5f %9.5f " % (atom[i], xyz[i][0], xyz[i][1],
                                                xyz[i][2]),
                    file=ofile)
        if args.tailormade4:
            print(
                "****PRINTING .xyz TAILOR-MADE4 FOR Qeq program by B.Wells***",
                file=ofile)
            print(
                "      FRAC       %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
                (ABC[0], ABC[1], ABC[2], math.degrees(abc[0]),
                 math.degrees(abc[1]), math.degrees(abc[2])),
                file=ofile)
            print("%d" % (natoms), file=ofile)
            for i in range(0, natoms):
                print(
                    "%3s %9.5f %9.5f %9.5f " % (atom[i], fract[i][0],
                                                fract[i][1], fract[i][2]),
                    file=ofile)
        if args.tailormade5:
            print(
                "****PRINTING .xyz TAILOR-MADE5 FOR Qeq program by B.Wells (with FC=0)***",
                file=ofile)
            print(
                "      FRAC       %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
                (ABC[0], ABC[1], ABC[2], math.degrees(abc[0]),
                 math.degrees(abc[1]), math.degrees(abc[2])),
                file=ofile)
            print("%d" % (natoms), file=ofile)
            for i in range(0, natoms):
                print(
                    "%3s %9.5f %9.5f %9.5f   xx   0.000   0.000" %
                    (atom[i], fract[i][0], fract[i][1], fract[i][2]),
                    file=ofile)

    if outputformat == "pwi":
        if not args.silent:
            print("QE input .pwi using the pseudo: %s" % (args.pseudopw))
        if not args.silent: print()
        print(" &CONTROL ", file=ofile)
        print("    calculation = 'vc-relax' ", file=ofile)
        print("    verbosity   = 'high' ", file=ofile)
        print("    !restart_mode= 'restart' ", file=ofile)
        print("    wf_collect  = .true. ", file=ofile)
        print("    outdir      = './' ", file=ofile)
        print("    prefix      = 'pwscf' ", file=ofile)
        print(
            "    pseudo_dir  = '/scratch/ongari/0_LIBRARIES/2_espresso/%s.1.0.0' "
            % (args.pseudopw),
            file=ofile)
        print("      !nstep        = 50", file=ofile)
        print(
            "      !etot_conv_thr= 1.0D-4", file=ofile
        )  # Note that etot_conv_thr is extensive: it can be hard to converge for big systems!
        print(
            "      !forc_conv_thr= 1.0D-3", file=ofile
        )  # Note that forc_conv_thr is intensive: Ok if max_forc < forc_conv_thr.
        print("    !max_seconds  = 3500", file=ofile)
        print("    !disk_io      = 'none'", file=ofile)
        print(" / ", file=ofile)
        print(" &SYSTEM ", file=ofile)
        print("    ibrav = 0 ", file=ofile)
        print("    nat   = %d " % (natoms), file=ofile)
        print("    ntyp  = %d " % (ntypes), file=ofile)
        print("    !nosym  = .true. ", file=ofile)
        print("      ecutwfc = 70 ", file=ofile)
        print("      ecutrho = 350 ", file=ofile)
        print("    occupations = 'smearing' ", file=ofile)
        print("    smearing    = 'gaussian' ", file=ofile)
        print("    degauss     = 0.02 ", file=ofile)
        print("      tot_charge                 = 0.0 ", file=ofile)
        print("      nspin                      = 1 ", file=ofile)
        print("      tot_magnetization          = -1 ", file=ofile)
        print("      !starting_magnetization(i) = +0.1 ", file=ofile)
        print("    !vdw_corr    = 'grimme-d2' ", file=ofile)
        print(
            "    !input_dft   = 'vdw-df2-b86r' ", file=ofile
        )  # REMEMBER TO: generate_vdW_kernel_table.x     (vc-relax not implemented for nspin>1)
        print(
            "    !input_dft   = 'rvv10' ", file=ofile
        )  # REMEMBER TO: generate_rVV10_kernel_table.x   (vc-relax ok)
        print(" / ", file=ofile)
        print(" &ELECTRONS ", file=ofile)
        print("    scf_must_converge = .false. ", file=ofile)
        print("    electron_maxstep  = 100 ", file=ofile)
        print("    conv_thr          = 1.0d-6 ", file=ofile)
        print("    mixing_mode       = 'local-TF' ", file=ofile)
        print("    mixing_beta       = 0.7 ", file=ofile)
        print("    diagonalization   = 'david' ", file=ofile)
        print(" / ", file=ofile)
        print(" &IONS ", file=ofile)
        print("    ion_dynamics = 'bfgs' ", file=ofile)
        print(" / ", file=ofile)
        print(" &CELL ", file=ofile)
        print("    cell_dynamics  = 'bfgs' ", file=ofile)
        print("    press          = 0.0 ", file=ofile)
        print(
            "    !cell_factor    = 2.0 ", file=ofile
        )  #default: 1.0 for all, 2.0 for vc-relax/md (Not enough space allocated for radial FFT: try restarting with a larger cell_factor.)
        print("    press_conv_thr = 0.5 ", file=ofile)
        print("    cell_dofree    = 'all' ", file=ofile)
        print(" / ", file=ofile)
        print("ATOMIC_SPECIES ", file=ofile)
        for i in range(0, len(atom_count)):
            if atom_count[i] != 0:
                print(
                    "%3s %8.3f  %s" % (atomic_symbol[i], atomic_mass[i],
                                       atomic_pseudo[args.pseudopw][i]),
                    file=ofile)  #add pseudo!
        print(" ", file=ofile)
        print("K_POINTS gamma ", file=ofile)
        print(" ", file=ofile)
        print(
            "CELL_PARAMETERS angstrom ", file=ofile
        )  #It should be very precise (http://pw_forum.pwscf.narkive.com/26uqaajr/crash-in-routine-set-sym-bl)
        print(
            "%11.8f %11.8f %11.8f" % (cell.item((0, 0)), cell.item(
                (0, 1)), cell.item((0, 2))),
            file=ofile)
        print(
            "%11.8f %11.8f %11.8f" % (cell.item((1, 0)), cell.item(
                (1, 1)), cell.item((1, 2))),
            file=ofile)
        print(
            "%11.8f %11.8f %11.8f" % (cell.item((2, 0)), cell.item(
                (2, 1)), cell.item((2, 2))),
            file=ofile)
        print(" ", file=ofile)
        print("ATOMIC_POSITIONS angstrom ", file=ofile)
        for i in range(0, natoms):
            print(
                "%3s %12.8f %12.8f %12.8f " % (atom[i], xyz[i][0], xyz[i][1],
                                               xyz[i][2]),
                file=ofile)

    if outputformat == "subsys":
        print(
            "##### Include it to the main cp2k.inp using: @INCLUDE '%s.subsys'"
            % outputfilename,
            file=ofile)
        print("  &SUBSYS", file=ofile)
        print("    &CELL", file=ofile)
        print("      PERIODIC XYZ", file=ofile)
        print("      MULTIPLE_UNIT_CELL 1 1 1", file=ofile)
        print("      SYMMETRY NONE", file=ofile)
        print(
            "      A [angstrom] %8.5f %8.5f %8.5f" % (cell.item(
                (0, 0)), cell.item((0, 1)), cell.item((0, 2))),
            file=ofile)
        print(
            "      B [angstrom] %8.5f %8.5f %8.5f" % (cell.item(
                (1, 0)), cell.item((1, 1)), cell.item((1, 2))),
            file=ofile)
        print(
            "      C [angstrom] %8.5f %8.5f %8.5f" % (cell.item(
                (2, 0)), cell.item((2, 1)), cell.item((2, 2))),
            file=ofile)
        print("    &END CELL", file=ofile)
        print(" ", file=ofile)
        print("    &COORD", file=ofile)
        print("      SCALED .FALSE.", file=ofile)
        for i in range(0, natoms):
            print(
                "%3s %9.5f %9.5f %9.5f " % (atom[i], xyz[i][0], xyz[i][1],
                                            xyz[i][2]),
                file=ofile)
        print("    &END COORD", file=ofile)
        print(" ", file=ofile)
        for i in range(0, len(atom_count)):
            if atom_count[i] != 0:
                print("    &KIND %s" % (atomic_symbol[i]), file=ofile)
                print("      BASIS_SET %s" % (args.bscp2k), file=ofile)
                print("      POTENTIAL GTH-PBE", file=ofile)
                print("    &END KIND", file=ofile)
                print("    &KIND %s_GHOST" % (atomic_symbol[i]), file=ofile)
                print("      BASIS_SET %s" % (args.bscp2k), file=ofile)
                print("      GHOST", file=ofile)
                print("    &END KIND", file=ofile)
                print(" ", file=ofile)
        print("  &END SUBSYS", file=ofile)

    if outputformat == "axsf":
        print("ANIMSTEPS 1", file=ofile)
        print("CRYSTAL", file=ofile)
        print("PRIMVEC 1", file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((0, 0)), cell.item(
                (0, 1)), cell.item((0, 2))),
            file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((1, 0)), cell.item(
                (1, 1)), cell.item((1, 2))),
            file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((2, 0)), cell.item(
                (2, 1)), cell.item((2, 2))),
            file=ofile)
        print("PRIMCOORD 1", file=ofile)
        print("%d 1" % (natoms), file=ofile)
        for i in range(0, natoms):
            print(
                "%3s %8.3f %8.3f %8.3f " % (atom[i], xyz[i][0], xyz[i][1],
                                            xyz[i][2]),
                file=ofile)

    if outputformat == "geo":
        print("Made with manage_crystal.py", file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((0, 0)), cell.item(
                (0, 1)), cell.item((0, 2))),
            file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((1, 0)), cell.item(
                (1, 1)), cell.item((1, 2))),
            file=ofile)
        print(
            "     %8.5f %8.5f %8.5f" % (cell.item((2, 0)), cell.item(
                (2, 1)), cell.item((2, 2))),
            file=ofile)
        print("%d " % (natoms), file=ofile)
        for i in range(0, natoms):
            print(
                "%d %8.3f %8.3f %8.3f %8.3f" % (an[i], xyz[i][0], xyz[i][1],
                                                xyz[i][2], charge[i]),
                file=ofile)

    ofile.close()

############################################################################## Write supplementary outputs
"""
if args.printatoms!=None:
   ofile=open(args.printatoms, 'w+')

   if not args.silent: print("*** Printing atoms to", args.printatoms)
   if not args.silent: print("*** 1st line: all atoms")
   if not args.silent: print("*** 2nd line: all atoms except for H, C, O")
   if not args.silent: print()

   for i in range(0,len(atom_count)):
	if atom_count[i] != 0:
           print(atomic_symbol[i],end='_',file=ofile)
   print("",file=ofile)

   for i in range(0,len(atom_count)):
	if atom_count[i] != 0 and atomic_symbol[i]!="H" \
                              and atomic_symbol[i]!="C" \
                              and atomic_symbol[i]!="O":
           print(atomic_symbol[i],end='_',file=ofile)
   print("",file=ofile)
   ofile.close()
"""
