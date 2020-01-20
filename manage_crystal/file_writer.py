"""Writing instructions for different formats.
Functions are sorted in alphabetic order.
"""

from manage_crystal.periodic_table import ptab_mass, ptab_qepseudo
from manage_crystal import Crys
import os
import sys


def write_to_filepath(crys, filepath, tm=None, pseudopw=None, bscp2k=None, potcp2k=None, fract=None):
    outputformat = os.path.splitext(filepath)[1][1:]
    ofile = open(filepath, 'w+')
    if outputformat == "axsf":
        write_axsf(ofile, crys)
    elif outputformat == "cif":
        write_cif(ofile, crys)
    elif outputformat == "cssr":
        write_cssr(ofile, crys)
    elif outputformat == "pdb":
        write_pdb(ofile, crys)
    elif outputformat == "poscar":
        write_poscar(ofile, crys, fract)
    elif outputformat == "pwi":
        write_pwi(ofile, crys, pseudopw)
    elif outputformat == "subsys":
        write_subsys(ofile, crys, bscp2k, potcp2k, fract)
    elif outputformat == "xyz" and tm == 0:
        write_xyz(ofile, crys)
    elif outputformat == "xyz" and tm == 4:
        write_xyz_tm4(ofile, crys)
    else:
        sys.exit("WARNING: Output file format not implemented. EXIT.")
    ofile.close()
    return


def write_axsf(ofile, c):
    ''' Write the Crys in .axsf file format'''
    # Better than .xsf fecause it can be readen (with cell) bt VMD
    print("ANIMSTEPS 1", file=ofile)
    print("CRYSTAL", file=ofile)
    print("PRIMVEC 1", file=ofile)
    for k in range(3):
        print("     %8.5f %8.5f %8.5f" % (c.matrix[k][0], c.matrix[k][1], c.matrix[k][2]), file=ofile)
    print("PRIMCOORD 1", file=ofile)
    print("%d 1" % (c.natom), file=ofile)
    for i in range(c.natom):
        print("%3s %8.3f %8.3f %8.3f " % (c.atom_element[i], c.atom_xyz[i][0], c.atom_xyz[i][1], c.atom_xyz[i][2]),
              file=ofile)
    return


def write_cif(ofile, c):
    ''' Write the Crys in .cif file format'''
    # BC: don't modify, if possible, to be coherent with git deposited files

    # if all charges are 0 do not print them
    print_charges = (max([abs(x) for x in c.atom_charge]) > 0)

    print("data_crystal", file=ofile)
    print(" ", file=ofile)
    print("_cell_length_a    %.5f" % c.length[0], file=ofile)
    print("_cell_length_b    %.5f" % c.length[1], file=ofile)
    print("_cell_length_c    %.5f" % c.length[2], file=ofile)
    print("_cell_angle_alpha %.5f" % c.angle_deg[0], file=ofile)
    print("_cell_angle_beta  %.5f" % c.angle_deg[1], file=ofile)
    print("_cell_angle_gamma %.5f" % c.angle_deg[2], file=ofile)
    print(" ", file=ofile)
    # To be versatile, P1 symmetry is specified in Hall, H-M and equiv_pos
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
    # Fractional cordinates are the newest standard for cifs
    print("_atom_site_fract_x", file=ofile)
    print("_atom_site_fract_y", file=ofile)
    print("_atom_site_fract_z", file=ofile)
    if print_charges:
        print("_atom_site_charge", file=ofile)
    for i in range(c.natom):
        print("{:10} {:5} {:>9.5f} {:>9.5f} {:>9.5f} ".format(c.atom_element[i], c.atom_element[i], c.atom_fract[i][0],
                                                              c.atom_fract[i][1], c.atom_fract[i][2]),
              end="",
              file=ofile)
        if print_charges:  # add charge with many decimals to avoid net charge
            print('{:>14.10f}'.format(c.atom_charge[i]), file=ofile)
        else:  # add \n
            print('', file=ofile)
    return


def write_cssr(ofile, c):
    ''' Write the Crys in .cssr file format'''
    print("                               %.3f  %.3f  %.3f" % (c.length[0], c.length[1], c.length[2]), file=ofile)
    print("                %.3f   %.3f   %.3f   SPGR =  1 P 1         OPT = 1" %
          (c.angle_deg[0], c.angle_deg[1], c.angle_deg[2]),
          file=ofile)
    print("%d   0" % (c.natom), file=ofile)
    print("0 %s       : %s" % ("xxxxxxx", "xxxxxxx"), file=ofile)
    for i in range(c.natom):
        print("%4d %3s %8.5f %8.5f %8.5f    0  0  0  0  0  0  0  0  %7.5f" %
              (i + 1, c.atom_element[i], c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2], c.atom_charge[i]),
              file=ofile)
    return


def write_pdb(ofile, c):
    ''' Write the Crys in .pdb file format'''
    print("CRYST1{0:>9.3f}{1:>9.3f}{2:>9.3f}{3:>7.2f}{4:>7.2f}{5:>7.2f} P 1           1".format(
        c.length[0], c.length[1], c.length[2], c.angle_deg[0], c.angle_deg[1], c.angle_deg[2]),
          file=ofile)
    for i in range(c.natom):
        print("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.5f%8.5f%8.5f%6.2f%6.2f          %2s%2s" %
              ("ATOM", i + 1, c.atom_element[i], "", "XXX", "X", 1, "", c.atom_fract[i][0], c.atom_fract[i][1],
               c.atom_fract[i][2], 1.00, 0.00, c.atom_element[i], ""),
              file=ofile)
    return


def write_poscar(ofile, c, fract):
    ''' Write POSCAR file for VASP
    NOTE: the order of atoms changes because of the sorting!
    '''
    # line 1: comment
    print("Written with manage_crystal", file=ofile)
    # line 2: scaling factor for the cell
    print("1", file=ofile)
    # line 3-5: cell matrix
    for k in range(3):
        print("%11.8f %11.8f %11.8f" % (c.matrix[k][0], c.matrix[k][1], c.matrix[k][2]), file=ofile)
    # line 6: list of elements
    for e in c.element:
        print("{} ".format(e), end="", file=ofile)
    print("", file=ofile)
    # line 7: list of elements' count
    for e in c.element:
        print("{} ".format(c.element_count[e]), end="", file=ofile)
    print("", file=ofile)
    # line 8: tell to read fractional (direct) or cartesian coordinates
    if fract:
        print("direct", file=ofile)
    else:
        print("cartesian", file=ofile)
    # line 9+: write coordinates following the order of sorted atoms
    for e in c.element:
        for i, a in enumerate(c.atom_element):
            if e == a:
                if fract:
                    print("%9.5f %9.5f %9.5f" % (c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2]),
                          file=ofile)
                else:
                    print("%9.5f %9.5f %9.5f" % (c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2]),
                          file=ofile)
    return


def write_pwi(ofile, c, pseudopw):
    ''' Write the Crys in .pwi file format for Quantum Espresso '''
    print(" &CONTROL ", file=ofile)
    print("    calculation = 'vc-relax' ", file=ofile)
    print("    verbosity   = 'high' ", file=ofile)
    print("    !restart_mode= 'restart' ", file=ofile)
    print("    wf_collect  = .true. ", file=ofile)
    print("    outdir      = './' ", file=ofile)
    print("    prefix      = 'pwscf' ", file=ofile)
    print("    pseudo_dir  = '/home/ongari/aiida-database/data/qe/%s' " % pseudopw, file=ofile)
    print("      !nstep        = 50", file=ofile)
    print("      !etot_conv_thr= 1.0D-4", file=ofile)
    # Note that etot_conv_thr is extensive: it can be hard to converge for big systems!
    print("      !forc_conv_thr= 1.0D-3", file=ofile)
    # Note that forc_conv_thr is intensive: Ok if max_forc < forc_conv_thr.
    print("    !max_seconds  = 3500", file=ofile)
    print("    !disk_io      = 'none'", file=ofile)
    print(" / ", file=ofile)
    print(" &SYSTEM ", file=ofile)
    print("    ibrav = 0 ", file=ofile)
    print("    nat   = %d " % (c.natom), file=ofile)
    print("    ntyp  = %d " % (c.nelement), file=ofile)
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
    print("    !input_dft   = 'vdw-df2-b86r' ", file=ofile)
    # REMEMBER TO: generate_vdW_kernel_table.x; (vc-relax not implemented for nspin>1)
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
    print("    press_conv_thr = 0.5 ", file=ofile)
    print("    cell_dofree    = 'all' ", file=ofile)
    print(" / ", file=ofile)
    print("ATOMIC_SPECIES ", file=ofile)
    for i, element in enumerate(c.element):
        print("%3s %8.3f  %s" % (element, ptab_mass[element], ptab_qepseudo[pseudopw][element]), file=ofile)
    print("", file=ofile)
    print("K_POINTS gamma ", file=ofile)
    print("", file=ofile)
    print("CELL_PARAMETERS angstrom ", file=ofile)
    #It should be very precise (http://pw_forum.pwscf.narkive.com/26uqaajr/crash-in-routine-set-sym-bl)
    for k in range(3):
        print("     %11.8f %11.8f %11.8f" % (c.matrix[k][0], c.matrix[k][1], c.matrix[k][2]), file=ofile)
    print("", file=ofile)
    print("ATOMIC_POSITIONS angstrom ", file=ofile)
    for i in range(c.natom):
        print("%3s %12.8f %12.8f %12.8f " % (c.atom_element[i], c.atom_xyz[i][0], c.atom_xyz[i][1], c.atom_xyz[i][2]),
              file=ofile)
    return


def write_subsys(ofile, c, bscp2k, potcp2k, fract):
    ''' Write the Crys in .subsys file format for CP2K'''
    print("# Include to the main cp2k.inp using: @INCLUDE 'filename.subsys'", file=ofile)
    print("", file=ofile)
    print("  &SUBSYS", file=ofile)
    print("    &CELL", file=ofile)
    for k, label in enumerate(["A", "B", "C"]):
        print("      %s [angstrom] %8.5f %8.5f %8.5f" % (label, c.matrix[k][0], c.matrix[k][1], c.matrix[k][2]),
              file=ofile)
    print("    &END CELL", file=ofile)
    print("", file=ofile)
    print("    &COORD", file=ofile)

    if fract:
        print("      SCALED .TRUE.", file=ofile)
        for i in range(c.natom):
            print("%3s %9.5f %9.5f %9.5f " %
                  (c.atom_element[i], c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2]),
                  file=ofile)
    else:
        for i in range(c.natom):
            print("%3s %9.5f %9.5f %9.5f " % (c.atom_element[i], c.atom_xyz[i][0], c.atom_xyz[i][1], c.atom_xyz[i][2]),
                  file=ofile)
    print("    &END COORD", file=ofile)
    print("", file=ofile)
    # print elements KIND
    for element in c.element:
        print("    &KIND %s" % element, file=ofile)
        print("      BASIS_SET %s" % bscp2k, file=ofile)
        print("      POTENTIAL %s" % potcp2k, file=ofile)
        print("    &END KIND", file=ofile)
        print("", file=ofile)
    # print element_ghost kind, for BSSE calculation
    for element in c.element:
        print("    &KIND %s_ghost" % element, file=ofile)
        print("      BASIS_SET %s" % bscp2k, file=ofile)
        print("      GHOST", file=ofile)
        print("    &END KIND", file=ofile)
        print("", file=ofile)
    print("  &END SUBSYS", file=ofile)
    return


def write_xyz(ofile, c):
    ''' Write the Crys in .xyz file format'''
    print("%d" % (c.natom), file=ofile)
    print("CELL:  %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
          (c.length[0], c.length[1], c.length[2], c.angle_deg[0], c.angle_deg[1], c.angle_deg[2]),
          file=ofile)
    for i in range(c.natom):
        print("%3s %9.5f %9.5f %9.5f " % (c.atom_element[i], c.atom_xyz[i][0], c.atom_xyz[i][1], c.atom_xyz[i][2]),
              file=ofile)
    return


def write_xyz_tm4(ofile, c):
    ''' Write the Crys in .xyz tilor-made 4 file format for Well's program'''
    print("****PRINTING .xyz TAILOR-MADE4 FOR Qeq program by B.Wells***", file=ofile)
    print("      FRAC       %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
          (c.length[0], c.length[1], c.length[2], c.angle_deg[0], c.angle_deg[1], c.angle_deg[2]),
          file=ofile)
    print("%d" % (c.natom), file=ofile)
    for i in range(c.natom):
        print("%3s %9.5f %9.5f %9.5f " %
              (c.atom_element[i], c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2]),
              file=ofile)
