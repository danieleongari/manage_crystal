from __future__ import absolute_import
from manage_crystal import Crys


def parse_cif(file):
    ''' Parse .cif file and return a Crys object '''
    # REQUIREMENTS:
    # - only valid for P1 symmetry
    # - cell data should be specified before the atom data
    # - if some other "_atom_something" are specified before, it does not work
    # - after the atom coordinates it breaks with "loop_" or EOF
    c = Crys()
    c.inp_lengths_angles = True
    while True:
        line = file.readline()
        data = line.split()
        if line == "":
            break
        if len(data) > 0 and (data[0] == "_cell_length_a"):
            c.length[0] = float(data[1])
        if len(data) > 0 and (data[0] == "_cell_length_b"):
            c.length[1] = float(data[1])
        if len(data) > 0 and (data[0] == "_cell_length_c"):
            c.length[2] = float(data[1])
        if len(data) > 0 and (data[0] == "_cell_angle_alpha"):
            c.angle_deg[0] = float(data[1])
        if len(data) > 0 and (data[0] == "_cell_angle_beta"):
            c.angle_deg[1] = float(data[1])
        if len(data) > 0 and (data[0] == "_cell_angle_gamma"):
            c.angle_deg[2] = float(data[1])
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
            break  #go in the loop to read coordinates
    # Read atomic element, coordinates and charges
    if "_atom_site_fract_x" in data_order_dic:
        c.inp_fract = True
    elif "_atom_site_Cartn_x" in data_order_dic:
        c.inp_xyz = True
    while True:
        if line == "" \
           or line =="\n" \
           or data[0] == "loop_" \
           or data[0] == "_loop":
            break
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
            c.atom_charge.append(
                float(data[data_order_dic["_atom_site_charge"]]))
        line = file.readline()
        data = line.split()
    return c


def parse_pdb(file):
    ''' Parse .pdb file and return Crys object '''
    c = Crys()
    c.inp_xyz = True
    while True:
        line = file.readline()
        data = line.split()
        if line == "":
            break
        elif len(data) > 0 and (data[0] == 'END' or data[0] == 'ENDMDL'):
            break
        elif len(data) > 0 and data[0] == 'CRYST1':
            c.inp_lengths_angles = True
            c.length = [
                float(line[0o6:15]),
                float(line[15:24]),
                float(line[24:33])
            ]
            c.angle_deg = [
                float(line[33:40]),
                float(line[40:47]),
                float(line[47:54])
            ]
        elif len(data) > 0 and (data[0] == "ATOM" or data[0] == "HETATM"):
            c.atom_type.append(data[2])  #maybe read to data[-1]
            c.atom_xyz.append(
                [float(line[30:38]),
                 float(line[38:46]),
                 float(line[46:54])])
    return c
