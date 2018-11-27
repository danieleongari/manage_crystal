# Writing instructions for different formats.
# Functions are sorted in alphabetic order.

from __future__ import print_function
from __future__ import absolute_import
from manage_crystal import Crys
from six.moves import range


def write_axsf(ofile, c):
    ''' Write the Crys in .axsf file format'''
    # Better than .xsf fecause it can be readen (with cell) bt VMD
    print("ANIMSTEPS 1", file=ofile)
    print("CRYSTAL", file=ofile)
    print("PRIMVEC 1", file=ofile)
    for k in range(3):
        print(
            "     %8.5f %8.5f %8.5f" % (c.matrix[k][0], c.matrix[k][1],
                                        c.matrix[k][2]),
            file=ofile)
    print("PRIMCOORD 1", file=ofile)
    print("%d 1" % (c.natom), file=ofile)
    for i in range(c.natom):
        print(
            "%3s %8.3f %8.3f %8.3f " % (c.atom_element[i], c.atom_xyz[i][0],
                                        c.atom_xyz[i][1], c.atom_xyz[i][2]),
            file=ofile)
    return


def write_cif(ofile, c):
    ''' Write the Crys in .cif file format'''
    # BC: don't modify, if possible, to be coherent with git deposited files
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
    print("_atom_site_charge", file=ofile)
    for i in range(c.natom):
        print(
            "{0:10} {1:5} {2:>9.5f} {3:>9.5f} {4:>9.5f} {5:>14.10f}".format(
                c.atom_element[i], c.atom_element[i], c.atom_fract[i][0],
                c.atom_fract[i][1], c.atom_fract[i][2],
                c.atom_charge[0]),  #many decimals to avoid net charge
            file=ofile)
    return


def write_cssr(ofile, c):
    ''' Write the Crys in .cssr file format'''
    print(
        "                               %.3f  %.3f  %.3f" %
        (c.length[0], c.length[1], c.length[2]),
        file=ofile)
    print(
        "                %.3f   %.3f   %.3f   SPGR =  1 P 1         OPT = 1" %
        (c.angle_deg[0], c.angle_deg[1], c.angle_deg[2]),
        file=ofile)
    print("%d   0" % (c.natom), file=ofile)
    print("0 %s       : %s" % ("xxxxxxx", "xxxxxxx"), file=ofile)
    for i in range(c.natom):
        print(
            "%4d %3s %8.5f %8.5f %8.5f    0  0  0  0  0  0  0  0  %7.5f" %
            (i + 1, c.atom_element[i], c.atom_fract[i][0], c.atom_fract[i][1],
             c.atom_fract[i][2], c.atom_charge[i][0]),
            file=ofile)
    return


def write_pdb(file, c):
    ''' Write the Crys in .pdb file format'''
    print(
        "CRYST1{0:>9.3f}{1:>9.3f}{2:>9.3f}{3:>7.2f}{4:>7.2f}{5:>7.2f} P 1           1"
        .format(c.length[0], c.length[1], c.length[2], c.angle_deg[0],
                c.angle_deg[1], c.angle_deg[2]),
        file=ofile)
    for i in range(c.natom):
        print(
            "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.5f%8.5f%8.5f%6.2f%6.2f          %2s%2s"
            % ("ATOM", i + 1, c.atom_element[i], "", "XXX", "X", 1, "",
               c.atom_fract[i][0], c.atom_fract[i][1], c.atom_fract[i][2],
               1.00, 0.00, c.atom_element[i], ""),
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
    print(
        "    pseudo_dir  = '/home/ongari/aiida-database/data/qe/%s' " %
        (args.pseudopw),
        file=ofile)
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
        print(
            "%3s %8.3f  %s" % (element, atomic_mass[c.element_atnum[i]],
                               atomic_pseudo[pseudopw][c.element_atnum[i]]),
            file=ofile)
    print("", file=ofile)
    print("K_POINTS gamma ", file=ofile)
    print("", file=ofile)
    print("CELL_PARAMETERS angstrom ", file=ofile)
    #It should be very precise (http://pw_forum.pwscf.narkive.com/26uqaajr/crash-in-routine-set-sym-bl)
    for k in range(3):
        print(
            "     %11.8f %11.8f %11.8f" % (c.matrix[k][0], c.matrix[k][1],
                                           c.matrix[k][2]),
            file=ofile)
    print("", file=ofile)
    print("ATOMIC_POSITIONS angstrom ", file=ofile)
    for i in range(c.natoms):
        print(
            "%3s %12.8f %12.8f %12.8f " % (c.atom_element[i], c.atom_xyz[i][0],
                                           c.atom_xyz[i][1], c.atom_xyz[i][2]),
            file=ofile)
    return


def write_subsys(ofile, c, outputfilename, bscp2k, potcp2k):
    ''' Write the Crys in .subsys file format for CP2K'''
    print(
        "##### Include it to the main cp2k.inp using: @INCLUDE '%s.subsys'" %
        outputfilename,
        file=ofile)
    print("  &SUBSYS", file=ofile)
    print("    &CELL", file=ofile)
    for k, label in enumerate(["A", "B", "C"]):
        print(
            "      %s [angstrom] %8.5f %8.5f %8.5f" %
            (label, c.matrix[k][0], c.matrix[k][1], c.matrix[k][2]),
            file=ofile)
    print("    &END CELL", file=ofile)
    print("", file=ofile)
    print("    &COORD", file=ofile)
    for i in range(c.natom):
        print(
            "%3s %9.5f %9.5f %9.5f " % (c.atom_element[i], c.aom_xyz[i][0],
                                        c.atom_xyz[i][1], c.atom_xyz[i][2]),
            file=ofile)
    print("    &END COORD", file=ofile)
    print("", file=ofile)
    for element in c.element:
        print("    &KIND %s" % element, file=ofile)
        print("      BASIS_SET %s" % bscp2k, file=ofile)
        print("      POTENTIAL %s" % potcp2k, file=ofile)
        print("    &END KIND", file=ofile)
        print("", file=ofile)
    print("  &END SUBSYS", file=ofile)
    return


def write_xyz(file, c):
    ''' Write the Crys in .xyz file format'''
    print("%d" % (c.natoms), file=ofile)
    print(
        "CELL:  %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
        (c.length[0], c.length[1], c.length[2], c.angle_deg[0], c.angle_deg[1],
         c.angle_deg[2]),
        file=ofile)
    for i in range(c.natom):
        print(
            "%3s %9.5f %9.5f %9.5f " % (c.atom_element[i], c.atom_xyz[i][0],
                                        c.atom_xyz[i][1], c.atom_xyz[i][2]),
            file=ofile)
    return


def write_xyz(file, c):
    ''' Write the Crys in .xyz tilor-made 4 file format for Well's program'''
    print(
        "****PRINTING .xyz TAILOR-MADE4 FOR Qeq program by B.Wells***",
        file=ofile)
    print(
        "      FRAC       %.5f  %.5f  %.5f  %.3f  %.3f  %.3f  " %
        (c.length[0], c.length[1], c.length[2], c.angle_deg[0], c.angle_deg[1],
         c.angle_deg[2]),
        file=ofile)
    print("%d" % (c.natom), file=ofile)
    for i in range(c.natom):
        print(
            "%3s %9.5f %9.5f %9.5f " %
            (c.atom_element[i], c.atom_fract[i][0], c.atom_fract[i][1],
             c.atom_fract[i][2]),
            file=ofile)
