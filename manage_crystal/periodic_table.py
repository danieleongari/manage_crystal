ptab_atnum = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I': 53,
    'Xe': 54,
    'Cs': 55,
    'Ba': 56,
    'La': 57,
    'Ce': 58,
    'Pr': 59,
    'Nd': 60,
    'Pm': 61,
    'Sm': 62,
    'Eu': 63,
    'Gd': 64,
    'Tb': 65,
    'Dy': 66,
    'Ho': 67,
    'Er': 68,
    'Tm': 69,
    'Yb': 70,
    'Lu': 71,
    'Hf': 72,
    'Ta': 73,
    'W': 74,
    'Re': 75,
    'Os': 76,
    'Ir': 77,
    'Pt': 78,
    'Au': 79,
    'Hg': 80,
    'Tl': 81,
    'Pb': 82,
    'Bi': 83,
    'Po': 84,
    'At': 85,
    'Rn': 86,
    'Fr': 87,
    'Ra': 88,
    'Ac': 89,
    'Th': 90,
    'Pa': 91,
    'U': 92,
    'Np': 93,
    'Pu': 94,
    'Am': 95,
    'Cm': 96,
    'Bk': 97,
    'Cf': 98,
    'Es': 99,
    'Fm': 100,
    'Md': 101,
    'No': 102,
    'Lr': 103,
    'Rf': 104,
    'Db': 105,
    'Sg': 106,
    'Bh': 107,
    'Hs': 108,
    'Mt': 109,
    'Ds': 110,
    'Rg': 111,
    'Cn': 112,
    'Nh': 113,
    'Fl': 114,
    'Mc': 115,
    'Lv': 116,
    'Ts': 117,
    'Og': 118,
}

ptab_atnum_inv = {v: k for k, v in ptab_atnum.items()}

ptab_fullname = {
    'H': 'Hydrogen',
    'He': 'Helium',
    'Li': 'Lithium',
    'Be': 'Beryllium',
    'B': 'Boron',
    'C': 'Carbon',
    'N': 'Nitrogen',
    'O': 'Oxygen',
    'F': 'Fluorine',
    'Ne': 'Neon',
    'Na': 'Sodium',
    'Mg': 'Magnesium',
    'Al': 'Aluminium',
    'Si': 'Silicon',
    'P': 'Phosphorus',
    'S': 'Sulfur',
    'Cl': 'Chlorine',
    'Ar': 'Argon',
    'K': 'Potassium',
    'Ca': 'Calcium',
    'Sc': 'Scandium',
    'Ti': 'Titanium',
    'V': 'Vanadium',
    'Cr': 'Chromium',
    'Mn': 'Manganese',
    'Fe': 'Iron',
    'Co': 'Cobalt',
    'Ni': 'Nickel',
    'Cu': 'Copper',
    'Zn': 'Zinc',
    'Ga': 'Gallium',
    'Ge': 'Germanium',
    'As': 'Arsenic',
    'Se': 'Selenium',
    'Br': 'Bromine',
    'Kr': 'Krypton',
    'Rb': 'Rubidium',
    'Sr': 'Strontium',
    'Y': 'Yttrium',
    'Zr': 'Zirconium',
    'Nb': 'Niobium',
    'Mo': 'Molybdenum',
    'Tc': 'Technetium',
    'Ru': 'Ruthenium',
    'Rh': 'Rhodium',
    'Pd': 'Palladium',
    'Ag': 'Silver',
    'Cd': 'Cadmium',
    'In': 'Indium',
    'Sn': 'Tin',
    'Sb': 'Antimony',
    'Te': 'Tellurium',
    'I': 'Iodine',
    'Xe': 'Xenon',
    'Cs': 'Caesium',
    'Ba': 'Barium',
    'La': 'Lanthanum',
    'Ce': 'Cerium',
    'Pr': 'Praseodymium',
    'Nd': 'Neodymium',
    'Pm': 'Promethium',
    'Sm': 'Samarium',
    'Eu': 'Europium',
    'Gd': 'Gadolinium',
    'Tb': 'Terbium',
    'Dy': 'Dysprosium',
    'Ho': 'Holmium',
    'Er': 'Erbium',
    'Tm': 'Thulium',
    'Yb': 'Ytterbium',
    'Lu': 'Lutetium',
    'Hf': 'Hafnium',
    'Ta': 'Tantalum',
    'W': 'Tungsten',
    'Re': 'Rhenium',
    'Os': 'Osmium',
    'Ir': 'Iridium',
    'Pt': 'Platinum',
    'Au': 'Gold',
    'Hg': 'Mercury',
    'Tl': 'Thallium',
    'Pb': 'Lead',
    'Bi': 'Bismuth',
    'Po': 'Polonium',
    'At': 'Astatine',
    'Rn': 'Radon',
    'Fr': 'Francium',
    'Ra': 'Radium',
    'Ac': 'Actinium',
    'Th': 'Thorium',
    'Pa': 'Protactinium',
    'U': 'Uranium',
    'Np': 'Neptunium',
    'Pu': 'Plutonium',
    'Am': 'Americium',
    'Cm': 'Curium',
    'Bk': 'Berkelium',
    'Cf': 'Californium',
    'Es': 'Einsteinium',
    'Fm': 'Fermium',
    'Md': 'Mendelevium',
    'No': 'Nobelium',
    'Lr': 'Lawrencium',
    'Rf': 'Rutherfordium',
    'Db': 'Dubnium',
    'Sg': 'Seaborgium',
    'Bh': 'Bohrium',
    'Hs': 'Hassium',
    'Mt': 'Meitnerium',
    'Ds': 'Darmstadtium',
    'Rg': 'Roentgenium',
    'Cn': 'Copernicium',
    'Nh': 'Nihonium',
    'Fl': 'Flerovium',
    'Mc': 'Moscovium',
    'Lv': 'Livermorium',
    'Ts': 'Tennessine',
    'Og': 'Oganesson',
}

# From UFF: L-J's sigma * 0.5
ptab_vdw_uff = {
    'H': 1.286,
    'He': 1.052,
    'Li': 1.092,
    'Be': 1.223,
    'B': 1.819,
    'C': 1.715,
    'N': 1.630,
    'O': 1.559,
    'F': 1.498,
    'Ne': 1.445,
    'Na': 1.329,
    'Mg': 1.346,
    'Al': 2.004,
    'Si': 1.913,
    'P': 1.847,
    'S': 1.797,
    'Cl': 1.758,
    'Ar': 1.723,
    'K': 1.698,
    'Ca': 1.514,
    'Sc': 1.468,
    'Ti': 1.414,
    'V': 1.400,
    'Cr': 1.347,
    'Mn': 1.319,
    'Fe': 1.297,
    'Co': 1.279,
    'Ni': 1.262,
    'Cu': 1.557,
    'Zn': 1.231,
    'Ga': 1.952,
    'Ge': 1.907,
    'As': 1.884,
    'Se': 1.873,
    'Br': 1.866,
    'Kr': 1.845,
    'Rb': 1.833,
    'Sr': 1.622,
    'Y': 1.490,
    'Zr': 1.392,
    'Nb': 1.410,
    'Mo': 1.360,
    'Tc': 1.335,
    'Ru': 1.320,
    'Rh': 1.305,
    'Pd': 1.291,
    'Ag': 1.402,
    'Cd': 1.269,
    'In': 1.988,
    'Sn': 1.956,
    'Sb': 1.969,
    'Te': 1.991,
    'I': 2.005,
    'Xe': 1.962,
    'Cs': 2.012,
    'Ba': 1.649,
    'La': 1.569,
    'Ce': 1.584,
    'Pr': 1.606,
    'Nd': 1.592,
    'Pm': 1.580,
    'Sm': 1.568,
    'Eu': 1.556,
    'Gd': 1.500,
    'Tb': 1.537,
    'Dy': 1.527,
    'Ho': 1.519,
    'Er': 1.511,
    'Tm': 1.503,
    'Yb': 1.494,
    'Lu': 1.621,
    'Hf': 1.399,
    'Ta': 1.412,
    'W': 1.367,
    'Re': 1.316,
    'Os': 1.390,
    'Ir': 1.265,
    'Pt': 1.227,
    'Au': 1.467,
    'Hg': 1.205,
    'Tl': 1.936,
    'Pb': 1.914,
    'Bi': 1.947,
    'Po': 2.098,
    'At': 2.116,
    'Rn': 2.123,
    'Fr': 2.183,
    'Ra': 1.638,
    'Ac': 1.549,
    'Th': 1.513,
    'Pa': 1.525,
    'U': 1.512,
    'Np': 1.525,
    'Pu': 1.525,
    'Am': 1.506,
    'Cm': 1.482,
    'Bk': 1.487,
    'Cf': 1.476,
    'Es': 1.470,
    'Fm': 1.464,
    'Md': 1.458,
    'No': 1.447,
    'Lr': 1.441,
    'Rf': 2.0,
    'Db': 2.0,
    'Sg': 2.0,
    'Bh': 2.0,
    'Hs': 2.0,
    'Mt': 2.0,
    'Ds': 2.0,
    'Rg': 2.0,
    'Cn': 2.0,
    'Nh': 2.0,
    'Fl': 2.0,
    'Mc': 2.0,
    'Lv': 2.0,
    'Ts': 2.0,
    'Og': 2.0,
}

# From ASE:
ptab_mass = {
    'H': 1.00794,
    'He': 4.00260,
    'Li': 6.94100,
    'Be': 9.01218,
    'B': 10.81100,
    'C': 12.01100,
    'N': 14.00670,
    'O': 15.99940,
    'F': 18.99840,
    'Ne': 20.17970,
    'Na': 22.98977,
    'Mg': 24.30500,
    'Al': 26.98154,
    'Si': 28.08550,
    'P': 30.97376,
    'S': 32.06600,
    'Cl': 35.45270,
    'Ar': 39.94800,
    'K': 39.09830,
    'Ca': 40.07800,
    'Sc': 44.95590,
    'Ti': 47.88000,
    'V': 50.94150,
    'Cr': 51.99600,
    'Mn': 54.93800,
    'Fe': 55.84700,
    'Co': 58.93320,
    'Ni': 58.69340,
    'Cu': 63.54600,
    'Zn': 65.39000,
    'Ga': 69.72300,
    'Ge': 72.61000,
    'As': 74.92160,
    'Se': 78.96000,
    'Br': 79.90400,
    'Kr': 83.80000,
    'Rb': 85.46780,
    'Sr': 87.62000,
    'Y': 88.90590,
    'Zr': 91.22400,
    'Nb': 92.90640,
    'Mo': 95.94000,
    'Tc': 1000000,
    'Ru': 101.07000,
    'Rh': 102.90550,
    'Pd': 106.42000,
    'Ag': 107.86800,
    'Cd': 112.41000,
    'In': 114.82000,
    'Sn': 118.71000,
    'Sb': 121.75700,
    'Te': 127.60000,
    'I': 126.90450,
    'Xe': 131.29000,
    'Cs': 132.90540,
    'Ba': 137.33000,
    'La': 138.90550,
    'Ce': 140.12000,
    'Pr': 140.90770,
    'Nd': 144.24000,
    'Pm': 1000000,
    'Sm': 150.36000,
    'Eu': 151.96500,
    'Gd': 157.25000,
    'Tb': 158.92530,
    'Dy': 162.50000,
    'Ho': 164.93030,
    'Er': 167.26000,
    'Tm': 168.93420,
    'Yb': 173.04000,
    'Lu': 174.96700,
    'Hf': 178.49000,
    'Ta': 180.94790,
    'W': 183.85000,
    'Re': 186.20700,
    'Os': 190.20000,
    'Ir': 192.22000,
    'Pt': 195.08000,
    'Au': 196.96650,
    'Hg': 200.59000,
    'Tl': 204.38300,
    'Pb': 207.20000,
    'Bi': 208.98040,
    'Po': 1000000,
    'At': 1000000,
    'Rn': 1000000,
    'Fr': 1000000,
    'Ra': 226.02540,
    'Ac': 1000000,
    'Th': 232.03810,
    'Pa': 231.03590,
    'U': 238.02900,
    'Np': 237.04820,
    'Pu': 1000000,
    'Am': 1000000,
    'Cm': 1000000,
    'Bk': 1000000,
    'Cf': 1000000,
    'Es': 1000000,
    'Fm': 1000000,
    'Md': 1000000,
    'No': 1000000,
    'Lr': 1000000,
    'Rf': 1000000,
    'Db': 1000000,
    'Sg': 1000000,
    'Bh': 1000000,
    'Hs': 1000000,
    'Mt': 1000000,
    'Ds': 1000000,
    'Rg': 1000000,
    'Cn': 1000000,
    'Nh': 1000000,
    'Fl': 1000000,
    'Mc': 1000000,
    'Lv': 1000000,
    'Ts': 1000000,
    'Og': 1000000,
}

# pseudoopotential: dictionary[pseudo_type][ptab_number]
ptab_qepseudo = {
    'pbe': {
        'H': ' H.pbe-kjpaw_psl.1.0.0.UPF   ! 	 46.Ry	 221.Ry',
        'He': ' ',
        'Li': ' ',
        'Be': ' ',
        'B': ' B.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 43.Ry	 325.Ry',
        'C': ' C.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 40.Ry	 326.Ry',  # C [1]
        'N': ' N.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 320.Ry',  # N [1]
        'O': ' O.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 47.Ry	 323.Ry',  # O [1]
        'F': ' ',
        'Ne': ' ',
        'Na': ' Na.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 66.Ry	 323.Ry',  # Na [1]
        'Mg': ' ',  # Mg [1]
        'Al': ' ',  # Al [5]
        'Si': ' ',  # Si [1]
        'P': ' ',  # P [1]
        'S': ' S.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 181.Ry',  # S [1]
        'Cl': ' Cl.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 223.Ry',  # Cl [1]
        'Ar': ' ',  # Ar [1]
        'K': ' ',  # K [1]
        'Ca': ' Ca.pbe-spn-kjpaw_psl.1.0.0.UPF !   45.Ry	 274.Ry',  # Ca [5]
        'Sc': ' ',  # Sc
        'Ti': ' Ti.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 52.Ry	 575.Ry',  # Ti
        'V': ' ',  # V
        'Cr': ' ',  # Cr
        'Mn': ' Mn.pbe-spn-kjpaw_psl.0.3.1.UPF ! 	 46.Ry	 244.Ry	 (v1.0.0 is not working!)',  # Mn
        'Fe': ' ',  # Fe
        'Co': ' ',  # Co
        'Ni': ' Ni.pbe-n-kjpaw_psl.1.0.0.UPF  ! 	 41.Ry	 236.Ry (spn available)',  # Ni [1]
        'Cu': ' Cu.pbe-dn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 236.Ry (spn available)',  # Cu [1]
        'Zn': ' Zn.pbe-dn-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 276.Ry (spn available)',  # Zn [1]
        'Ga': ' ',  # Ga [1]
        'Ge': ' ',  # Ge [5]
        'As': ' ',  # As [1]
        'Se': ' ',  # Se [1]
        'Br': ' ',  # Br [1]
        'Kr': ' ',  # Kr [1]
        'Rb': ' ',  # Rb [5]
        'Sr': ' Sr.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 262.Ry',  # Sr [5]
        'Y': ' ',  # Y
        'Zr': ' ',  # Zr
        'Nb': ' ',  # Nb
        'Mo': ' Mo.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 49.Ry	 306.Ry',  # Mo
        'Tc': ' ',  # Tc
        'Ru': ' ',  # Ru
        'Rh': ' ',  # Rh
        'Pd': ' ',  # Pd [1]
        'Ag': ' ',  # Ag [1]
        'Cd': ' Cd.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 37.Ry	 638.Ry (spn available)',  # Cd [1]
        'In': ' ',  # In [1]
        'Sn': ' ',  # Sn [1]
        'Sb': ' ',  # Sb [5]
        'Te': ' ',  # Te [1]
        'I': ' ',  # I [1]
        'Xe': ' ',  # Xe [1]
        'Cs': ' ',  # Cs [5]
        'Ba': ' ',  # Ba [5]
        'La': ' ',  # La
        'Ce': ' ',  # Ce
        'Pr': ' ',  # Pr
        'Nd': ' ',  # Nd
        'Pm': ' ',  # Pm
        'Sm': ' ',  # Sm
        'Eu': ' ',  # Eu
        'Gd': ' ',  # Gd
        'Tb': ' ',  # Tb
        'Dy': ' ',  # Dy
        'Ho': ' ',  # Ho
        'Er': ' ',  # Er
        'Tm': ' ',  # Tm
        'Yb': ' ',  # Yb
        'Lu': ' ',  # Lu
        'Hf': ' ',  # Hf
        'Ta': ' ',  # Ta
        'W': ' ',  # W
        'Re': ' ',  # Re
        'Os': ' ',  # Os
        'Ir': ' ',  # Ir
        'Pt': ' ',  # Pt [1]
        'Au': ' ',  # Au [1]
        'Hg': ' ',  # Hg [1]
        'Tl': ' ',  # Tl [1]
        'Pb': ' ',  # Pb [1]
        'Bi': ' ',  # Bi [5]
        'Po': ' ',  # Po [5]
        'At': ' ',  # At [5]
        'Rn': ' ',  # Rn [5]
        'Fr': ' ',  # Fr [5]
        'Ra': ' ',  # Ra [5]
        'Ac': ' ',  # Ac
        'Th': ' ',  # Th
        'Pa': ' ',  # Pa
        'U': ' ',  # U [1]
        'Np': ' ',  # Np
        'Pu': ' ',  # Pu
        'Am': ' ',  # Am
        'Cm': ' ',  # Cm
        'Bk': ' ',  # Bk
        'Cf': ' ',  # Cf
        'Es': ' ',  # Es
        'Fm': ' ',  # Fm
        'Md': ' ',  # Md
        'No': ' ',  # No
        'Lr': ' ',
        'Rf': ' ',
        'Db': ' ',
        'Sg': ' ',
        'Bh': ' ',
        'Hs': ' ',
        'Mt': ' ',
        'Ds': ' ',
        'Rg': ' ',
        'Cn': ' ',
        'Nh': ' ',
        'Fl': ' ',
        'Mc': ' ',
        'Lv': ' ',
        'Ts': ' ',
        'Og': ' ',
    },
    'pbesol': {
        'H': ' H.pbesol-kjpaw_psl.1.0.0.UPF   ! 	 46.Ry	 221.Ry',  # H
        'He': ' ',  # He [1]
        'Li': ' ',  # Li [1]
        'Be': ' ',  # Be [5]
        'B': ' B.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 43.Ry	 324.Ry',  # B [5]
        'C': ' C.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 40.Ry	 325.Ry',  # C [1]
        'N': ' N.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 364.Ry',  # N [1]
        'O': ' O.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 47.Ry	 323.Ry',  # O [1]
        'F': ' ',  # F [1]
        'Ne': ' ',  # Ne [1]
        'Na': ' Na.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 66.Ry	 323.Ry',  # Na [1]
        'Mg': ' ',  # Mg [1]
        'Al': ' ',  # Al [5]
        'Si': ' ',  # Si [1]
        'P': ' ',  # P [1]
        'S': ' S.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 181.Ry',  # S [1]
        'Cl': ' ',  # Cl [1]
        'Ar': ' ',  # Ar [1]
        'K': ' ',  # K [1]
        'Ca': ' Ca.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 274.Ry',  # Ca [5]
        'Sc': ' ',  # Sc
        'Ti': ' ',  # Ti
        'V': ' ',  # V
        'Cr': ' ',  # Cr
        'Mn': ' ',  # Mn
        'Fe': ' ',  # Fe
        'Co': ' Co.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 238.Ry',  # Co
        'Ni': ' Ni.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 41.Ry	 236.Ry (spn available)',  # Ni [1]
        'Cu': ' Cu.pbesol-dn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 236.Ry (spn available)',  # Cu [1]
        'Zn': ' Zn.pbesol-dn-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 276.Ry (spn available)',  # Zn [1]
        'Ga': ' ',  # Ga [1]
        'Ge': ' ',  # Ge [5]
        'As': ' ',  # As [1]
        'Se': ' ',  # Se [1]
        'Br': ' ',  # Br [1]
        'Kr': ' ',  # Kr [1]
        'Rb': ' ',  # Rb [5]
        'Sr': ' Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 262.Ry',  # Sr [5]
        'Y': ' ',  # Y
        'Zr': ' ',  # Zr
        'Nb': ' ',  # Nb
        'Mo': ' ',  # Mo
        'Tc': ' ',  # Tc
        'Ru': ' ',  # Ru
        'Rh': ' ',  # Rh
        'Pd': ' ',  # Pd [1]
        'Ag': ' ',  # Ag [1]
        'Cd': ' Cd.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 37.Ry	 638.Ry (spn available)',  # Cd [1]
        'In': ' ',  # In [1]
        'Sn': ' ',  # Sn [1]
        'Sb': ' ',  # Sb [5]
        'Te': ' ',  # Te [1]
        'I': ' ',  # I [1]
        'Xe': ' ',  # Xe [1]
        'Cs': ' ',  # Cs [5]
        'Ba': ' ',  # Ba [5]
        'La': ' ',  # La
        'Ce': ' ',  # Ce
        'Pr': ' ',  # Pr
        'Nd': ' ',  # Nd
        'Pm': ' ',  # Pm
        'Sm': ' ',  # Sm
        'Eu': ' ',  # Eu
        'Gd': ' ',  # Gd
        'Tb': ' ',  # Tb
        'Dy': ' ',  # Dy
        'Ho': ' ',  # Ho
        'Er': ' ',  # Er
        'Tm': ' ',  # Tm
        'Yb': ' ',  # Yb
        'Lu': ' ',  # Lu
        'Hf': ' ',  # Hf
        'Ta': ' ',  # Ta
        'W': ' ',  # W
        'Re': ' ',  # Re
        'Os': ' ',  # Os
        'Ir': ' ',  # Ir
        'Pt': ' ',  # Pt [1]
        'Au': ' ',  # Au [1]
        'Hg': ' ',  # Hg [1]
        'Tl': ' ',  # Tl [1]
        'Pb': ' ',  # Pb [1]
        'Bi': ' ',  # Bi [5]
        'Po': ' ',  # Po [5]
        'At': ' ',  # At [5]
        'Rn': ' ',  # Rn [5]
        'Fr': ' ',  # Fr [5]
        'Ra': ' ',  # Ra [5]
        'Ac': ' ',  # Ac
        'Th': ' ',  # Th
        'Pa': ' ',  # Pa
        'U': ' ',  # U [1]
        'Np': ' ',  # Np
        'Pu': ' ',  # Pu
        'Am': ' ',  # Am
        'Cm': ' ',  # Cm
        'Bk': ' ',  # Bk
        'Cf': ' ',  # Cf
        'Es': ' ',  # Es
        'Fm': ' ',  # Fm
        'Md': ' ',  # Md
        'No': ' ',  # No
        'Lr': ' ',
        'Rf': ' ',
        'Db': ' ',
        'Sg': ' ',
        'Bh': ' ',
        'Hs': ' ',
        'Mt': ' ',
        'Ds': ' ',
        'Rg': ' ',
        'Cn': ' ',
        'Nh': ' ',
        'Fl': ' ',
        'Mc': ' ',
        'Lv': ' ',
        'Ts': ' ',
        'Og': ' ',
    }
}
