atomic_symbol = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc',
                 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
                 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb',
                 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
                 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
                 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
                 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
                 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
                 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
                 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
                 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og')

atomic_name = (
    'DUMMY', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon',
    'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium', 'Magnesium',
    'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur', 'Chlorine', 'Argon',
    'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium',
    'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium',
    'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium',
    'Ruthenium', 'Rhodium', 'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin',
    'Antimony', 'Tellurium', 'Iodine', 'Xenon', 'Caesium', 'Barium',
    'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium',
    'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum',
    'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium',
    'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium', 'Neptunium',
    'Plutonium', 'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium',
    'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium',
    'Meitnerium', 'Darmstadtium', 'Roentgenium', 'Copernicium', 'Nihonium',
    'Flerovium', 'Moscovium', 'Livermorium', 'Tennessine', 'Oganesson')

atomic_vdw_ase = (
    1000000,  # Dummy
    1.20,  # H
    1.40,  # He [1]
    1.82,  # Li [1]
    1.53,  # Be [5]
    1.92,  # B [5]
    1.70,  # C [1]
    1.55,  # N [1]
    1.52,  # O [1]
    1.47,  # F [1]
    1.54,  # Ne [1]
    2.27,  # Na [1]
    1.73,  # Mg [1]
    1.84,  # Al [5]
    2.10,  # Si [1]
    1.80,  # P [1]
    1.80,  # S [1]
    1.75,  # Cl [1]
    1.88,  # Ar [1]
    2.75,  # K [1]
    2.31,  # Ca [5]
    1000000,  # Sc
    1000000,  # Ti
    1000000,  # V
    1000000,  # Cr
    1000000,  # Mn
    1000000,  # Fe
    1000000,  # Co
    1.63,  # Ni [1]
    1.40,  # Cu [1]
    1.39,  # Zn [1]
    1.87,  # Ga [1]
    2.11,  # Ge [5]
    1.85,  # As [1]
    1.90,  # Se [1]
    1.85,  # Br [1]
    2.02,  # Kr [1]
    3.03,  # Rb [5]
    2.49,  # Sr [5]
    1000000,  # Y
    1000000,  # Zr
    1000000,  # Nb
    1000000,  # Mo
    1000000,  # Tc
    1000000,  # Ru
    1000000,  # Rh
    1.63,  # Pd [1]
    1.72,  # Ag [1]
    1.58,  # Cd [1]
    1.93,  # In [1]
    2.17,  # Sn [1]
    2.06,  # Sb [5]
    2.06,  # Te [1]
    1.98,  # I [1]
    2.16,  # Xe [1]
    3.43,  # Cs [5]
    2.49,  # Ba [5]
    1000000,  # La
    1000000,  # Ce
    1000000,  # Pr
    1000000,  # Nd
    1000000,  # Pm
    1000000,  # Sm
    1000000,  # Eu
    1000000,  # Gd
    1000000,  # Tb
    1000000,  # Dy
    1000000,  # Ho
    1000000,  # Er
    1000000,  # Tm
    1000000,  # Yb
    1000000,  # Lu
    1000000,  # Hf
    1000000,  # Ta
    1000000,  # W
    1000000,  # Re
    1000000,  # Os
    1000000,  # Ir
    1.75,  # Pt [1]
    1.66,  # Au [1]
    1.55,  # Hg [1]
    1.96,  # Tl [1]
    2.02,  # Pb [1]
    2.07,  # Bi [5]
    1.97,  # Po [5]
    2.02,  # At [5]
    2.20,  # Rn [5]
    3.48,  # Fr [5]
    2.83,  # Ra [5]
    1000000,  # Ac
    1000000,  # Th
    1000000,  # Pa
    1.86,  # U [1]
    1000000,  # Np
    1000000,  # Pu
    1000000,  # Am
    1000000,  # Cm
    1000000,  # Bk
    1000000,  # Cf
    1000000,  # Es
    1000000,  # Fm
    1000000,  # Md
    1000000,  # No
    1000000)  # Lr

# From Zeo++:
atomic_vdw_zeopp = (
    1,  #Dummy
    1.09,  #H
    1.4,  #He
    1.82,  #Li
    2,  #Be
    2,  #B
    1.7,  #C
    1.55,  #N
    1.52,  #O
    1.47,  #F
    1.54,  #Ne
    2.27,  #Na
    1.73,  #Mg
    2,  #Al
    2.1,  #Si
    1.8,  #P
    1.8,  #S
    1.75,  #Cl
    1.88,  #Ar
    2.75,  #K
    2,  #Ca
    2,  #Sc
    2,  #Ti
    2,  #V
    2,  #Cr
    2,  #Mn
    2,  #Fe
    2,  #Co
    1.63,  #Ni
    1.4,  #Cu
    1.39,  #Zn
    1.87,  #Ga
    2,  #Ge
    1.85,  #As
    1.9,  #Se
    1.85,  #Br
    2.02,  #Kr
    2,  #Rb
    2,  #Sr
    2,  #Y
    2,  #Zr
    2,  #Nb
    2,  #Mo
    2,  #Tc
    2,  #Ru
    2,  #Rh
    1.63,  #Pd
    1.72,  #Ag
    1.58,  #Cd
    1.93,  #In
    2.17,  #Sn
    2,  #Sb
    2.06,  #Te
    1.98,  #I
    2.16,  #Xe
    2,  #Cs
    2,  #Ba
    2,  #La
    2,  #Ce
    2,  #Pr
    2,  #Nd
    2,  #Pm
    2,  #Sm
    2,  #Eu
    2,  #Gd
    2,  #Tb
    2,  #Dy
    2,  #Ho
    2,  #Er
    2,  #Tm
    2,  #Yb
    2,  #Lu
    2,  #Hf
    2,  #Ta
    2,  #W
    2,  #Re
    2,  #Os
    2,  #Ir
    1.72,  #Pt
    1.66,  #Au
    1.55,  #Hg
    1.96,  #Tl
    2.02,  #Pb
    2,  #Bi
    2,  #Po
    2,  #At
    2,  #Rn
    2,  #Fr
    2,  #Ra
    2,  #Ac
    2,  #Th
    2,  #Pa
    1.86,  #U
    2,  #Np
    2,  #Pu
    2,  #Am
    2,  #Cm
    2,  #Bk
    2,  #Cf
    2,  #Es
    2,  #Fm
    2,  #Md
    2,  #No
    2,  #Lr
    2,  #Rf
    2,  #Db
    2,  #Sg
    2,  #Bh
    2,  #Hs
    2,  #Mt
    2)  #Ds

# From UFF: L-J's sigma * 0.5
atomic_vdw_UFF = (
    1,  #Dummy
    1.286,  #H
    1.052,  #He
    1.092,  #Li
    1.223,  #Be
    1.819,  #B
    1.715,  #C
    1.630,  #N
    1.559,  #O
    1.498,  #F
    1.445,  #Ne
    1.329,  #Na
    1.346,  #Mg
    2.004,  #Al
    1.913,  #Si
    1.847,  #P
    1.797,  #S
    1.758,  #Cl
    1.723,  #Ar
    1.698,  #K
    1.514,  #Ca
    1.468,  #Sc
    1.414,  #Ti
    1.400,  #V
    1.347,  #Cr
    1.319,  #Mn
    1.297,  #Fe
    1.279,  #Co
    1.262,  #Ni
    1.557,  #Cu
    1.231,  #Zn
    1.952,  #Ga
    1.907,  #Ge
    1.884,  #As
    1.873,  #Se
    1.866,  #Br
    1.845,  #Kr
    1.833,  #Rb
    1.622,  #Sr
    1.490,  #Y
    1.392,  #Zr
    1.410,  #Nb
    1.360,  #Mo
    1.335,  #Tc
    1.320,  #Ru
    1.305,  #Rh
    1.291,  #Pd
    1.402,  #Ag
    1.269,  #Cd
    1.988,  #In
    1.956,  #Sn
    1.969,  #Sb
    1.991,  #Te
    2.005,  #I
    1.962,  #Xe
    2.012,  #Cs
    1.649,  #Ba
    1.569,  #La
    1.584,  #Ce
    1.606,  #Pr
    1.592,  #Nd
    1.580,  #Pm
    1.568,  #Sm
    1.556,  #Eu
    1.500,  #Gd
    1.537,  #Tb
    1.527,  #Dy
    1.519,  #Ho
    1.511,  #Er
    1.503,  #Tm
    1.494,  #Yb
    1.621,  #Lu
    1.399,  #Hf
    1.412,  #Ta
    1.367,  #W
    1.316,  #Re
    1.390,  #Os
    1.265,  #Ir
    1.227,  #Pt
    1.467,  #Au
    1.205,  #Hg
    1.936,  #Tl
    1.914,  #Pb
    1.947,  #Bi
    2.098,  #Po
    2.116,  #At
    2.123,  #Rn
    2.183,  #Fr
    1.638,  #Ra
    1.549,  #Ac
    1.513,  #Th
    1.525,  #Pa
    1.512,  #U
    1.525,  #Np
    1.525,  #Pu
    1.506,  #Am
    1.482,  #Cm
    1.487,  #Bk
    1.476,  #Cf
    1.470,  #Es
    1.464,  #Fm
    1.458,  #Md
    1.447,  #No
    1.441,  #Lr
    2,  #Rf
    2,  #Db
    2,  #Sg
    2,  #Bh
    2,  #Hs
    2,  #Mt
    2)  #Ds

# From ASE:
atomic_mass = {
    "H": 1.00794,
    "He": 4.00260,
    "Li": 6.94100,
    "Be": 9.01218,
    "B": 10.81100,
    "C": 12.01100,
    "N": 14.00670,
    "O": 15.99940,
    "F": 18.99840,
    "Ne": 20.17970,
    "Na": 22.98977,
    "Mg": 24.30500,
    "Al": 26.98154,
    "Si": 28.08550,
    "P": 30.97376,
    "S": 32.06600,
    "Cl": 35.45270,
    "Ar": 39.94800,
    "K": 39.09830,
    "Ca": 40.07800,
    "Sc": 44.95590,
    "Ti": 47.88000,
    "V": 50.94150,
    "Cr": 51.99600,
    "Mn": 54.93800,
    "Fe": 55.84700,
    "Co": 58.93320,
    "Ni": 58.69340,
    "Cu": 63.54600,
    "Zn": 65.39000,
    "Ga": 69.72300,
    "Ge": 72.61000,
    "As": 74.92160,
    "Se": 78.96000,
    "Br": 79.90400,
    "Kr": 83.80000,
    "Rb": 85.46780,
    "Sr": 87.62000,
    "Y": 88.90590,
    "Zr": 91.22400,
    "Nb": 92.90640,
    "Mo": 95.94000,
    "Tc": 1000000,
    "Ru": 101.07000,
    "Rh": 102.90550,
    "Pd": 106.42000,
    "Ag": 107.86800,
    "Cd": 112.41000,
    "In": 114.82000,
    "Sn": 118.71000,
    "Sb": 121.75700,
    "Te": 127.60000,
    "I": 126.90450,
    "Xe": 131.29000,
    "Cs": 132.90540,
    "Ba": 137.33000,
    "La": 138.90550,
    "Ce": 140.12000,
    "Pr": 140.90770,
    "Nd": 144.24000,
    "Pm": 1000000,
    "Sm": 150.36000,
    "Eu": 151.96500,
    "Gd": 157.25000,
    "Tb": 158.92530,
    "Dy": 162.50000,
    "Ho": 164.93030,
    "Er": 167.26000,
    "Tm": 168.93420,
    "Yb": 173.04000,
    "Lu": 174.96700,
    "Hf": 178.49000,
    "Ta": 180.94790,
    "W ": 183.85000,
    "Re": 186.20700,
    "Os": 190.20000,
    "Ir": 192.22000,
    "Pt": 195.08000,
    "Au": 196.96650,
    "Hg": 200.59000,
    "Tl": 204.38300,
    "Pb": 207.20000,
    "Bi": 208.98040,
    "Po": 1000000,
    "At": 1000000,
    "Rn": 1000000,
    "Fr": 1000000,
    "Ra": 226.02540,
    "Ac": 1000000,
    "Th": 232.03810,
    "Pa": 231.03590,
    "U ": 238.02900,
    "Np": 237.04820,
    "Pu": 1000000,
    "Am": 1000000,
    "Cm": 1000000,
    "Bk": 1000000,
    "Cf": 1000000,
    "Es": 1000000,
    "Fm": 1000000,
    "Md": 1000000,
    "No": 1000000,
    "Lr": 1000000
}

# pseudoopotential: dictionary[pseudo_type][atomic_number]
atomic_pseudo = {
    'pbe': [
        ' ',  # X
        ' H.pbe-kjpaw_psl.1.0.0.UPF   ! 	 46.Ry	 221.Ry',  # H
        ' ',  # He [1]
        ' ',  # Li [1]
        ' ',  # Be [5]
        ' B.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 43.Ry	 325.Ry',  # B [5]
        ' C.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 40.Ry	 326.Ry',  # C [1]
        ' N.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 320.Ry',  # N [1]
        ' O.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 47.Ry	 323.Ry',  # O [1]
        ' ',  # F [1]
        ' ',  # Ne [1]
        ' Na.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 66.Ry	 323.Ry',  # Na [1]
        ' ',  # Mg [1]
        ' ',  # Al [5]
        ' ',  # Si [1]
        ' ',  # P [1]
        ' S.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 181.Ry',  # S [1]
        ' Cl.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 223.Ry',  # Cl [1]
        ' ',  # Ar [1]
        ' ',  # K [1]
        ' Ca.pbe-spn-kjpaw_psl.1.0.0.UPF !   45.Ry	 274.Ry',  # Ca [5]
        ' ',  # Sc
        ' Ti.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 52.Ry	 575.Ry',  # Ti
        ' ',  # V
        ' ',  # Cr
        ' Mn.pbe-spn-kjpaw_psl.0.3.1.UPF ! 	 46.Ry	 244.Ry	 (v1.0.0 is not working!)',  # Mn
        ' ',  # Fe
        ' ',  # Co
        ' Ni.pbe-n-kjpaw_psl.1.0.0.UPF  ! 	 41.Ry	 236.Ry (spn available)',  # Ni [1]
        ' Cu.pbe-dn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 236.Ry (spn available)',  # Cu [1]
        ' Zn.pbe-dn-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 276.Ry (spn available)',  # Zn [1]
        ' ',  # Ga [1]
        ' ',  # Ge [5]
        ' ',  # As [1]
        ' ',  # Se [1]
        ' ',  # Br [1]
        ' ',  # Kr [1]
        ' ',  # Rb [5]
        ' Sr.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 262.Ry',  # Sr [5]
        ' ',  # Y
        ' ',  # Zr
        ' ',  # Nb
        ' Mo.pbe-spn-kjpaw_psl.1.0.0.UPF ! 	 49.Ry	 306.Ry',  # Mo
        ' ',  # Tc
        ' ',  # Ru
        ' ',  # Rh
        ' ',  # Pd [1]
        ' ',  # Ag [1]
        ' Cd.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 37.Ry	 638.Ry (spn available)',  # Cd [1]
        ' ',  # In [1]
        ' ',  # Sn [1]
        ' ',  # Sb [5]
        ' ',  # Te [1]
        ' ',  # I [1]
        ' ',  # Xe [1]
        ' ',  # Cs [5]
        ' ',  # Ba [5]
        ' ',  # La
        ' ',  # Ce
        ' ',  # Pr
        ' ',  # Nd
        ' ',  # Pm
        ' ',  # Sm
        ' ',  # Eu
        ' ',  # Gd
        ' ',  # Tb
        ' ',  # Dy
        ' ',  # Ho
        ' ',  # Er
        ' ',  # Tm
        ' ',  # Yb
        ' ',  # Lu
        ' ',  # Hf
        ' ',  # Ta
        ' ',  # W
        ' ',  # Re
        ' ',  # Os
        ' ',  # Ir
        ' ',  # Pt [1]
        ' ',  # Au [1]
        ' ',  # Hg [1]
        ' ',  # Tl [1]
        ' ',  # Pb [1]
        ' ',  # Bi [5]
        ' ',  # Po [5]
        ' ',  # At [5]
        ' ',  # Rn [5]
        ' ',  # Fr [5]
        ' ',  # Ra [5]
        ' ',  # Ac
        ' ',  # Th
        ' ',  # Pa
        ' ',  # U [1]
        ' ',  # Np
        ' ',  # Pu
        ' ',  # Am
        ' ',  # Cm
        ' ',  # Bk
        ' ',  # Cf
        ' ',  # Es
        ' ',  # Fm
        ' ',  # Md
        ' ',  # No
        ' '
    ],  # Lr
    'pbesol': [
        ' ',  # X
        ' H.pbesol-kjpaw_psl.1.0.0.UPF   ! 	 46.Ry	 221.Ry',  # H
        ' ',  # He [1]
        ' ',  # Li [1]
        ' ',  # Be [5]
        ' B.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 43.Ry	 324.Ry',  # B [5]
        ' C.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 40.Ry	 325.Ry',  # C [1]
        ' N.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 364.Ry',  # N [1]
        ' O.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 47.Ry	 323.Ry',  # O [1]
        ' ',  # F [1]
        ' ',  # Ne [1]
        ' Na.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 66.Ry	 323.Ry',  # Na [1]
        ' ',  # Mg [1]
        ' ',  # Al [5]
        ' ',  # Si [1]
        ' ',  # P [1]
        ' S.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 181.Ry',  # S [1]
        ' ',  # Cl [1]
        ' ',  # Ar [1]
        ' ',  # K [1]
        ' Ca.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 274.Ry',  # Ca [5]
        ' ',  # Sc
        ' ',  # Ti
        ' ',  # V
        ' ',  # Cr
        ' ',  # Mn
        ' ',  # Fe
        ' Co.pbe-n-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 238.Ry',  # Co
        ' Ni.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 41.Ry	 236.Ry (spn available)',  # Ni [1]
        ' Cu.pbesol-dn-kjpaw_psl.1.0.0.UPF ! 	 45.Ry	 236.Ry (spn available)',  # Cu [1]
        ' Zn.pbesol-dn-kjpaw_psl.1.0.0.UPF ! 	 44.Ry	 276.Ry (spn available)',  # Zn [1]
        ' ',  # Ga [1]
        ' ',  # Ge [5]
        ' ',  # As [1]
        ' ',  # Se [1]
        ' ',  # Br [1]
        ' ',  # Kr [1]
        ' ',  # Rb [5]
        ' Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF ! 	 39.Ry	 262.Ry',  # Sr [5]
        ' ',  # Y
        ' ',  # Zr
        ' ',  # Nb
        ' ',  # Mo
        ' ',  # Tc
        ' ',  # Ru
        ' ',  # Rh
        ' ',  # Pd [1]
        ' ',  # Ag [1]
        ' Cd.pbesol-n-kjpaw_psl.1.0.0.UPF ! 	 37.Ry	 638.Ry (spn available)',  # Cd [1]
        ' ',  # In [1]
        ' ',  # Sn [1]
        ' ',  # Sb [5]
        ' ',  # Te [1]
        ' ',  # I [1]
        ' ',  # Xe [1]
        ' ',  # Cs [5]
        ' ',  # Ba [5]
        ' ',  # La
        ' ',  # Ce
        ' ',  # Pr
        ' ',  # Nd
        ' ',  # Pm
        ' ',  # Sm
        ' ',  # Eu
        ' ',  # Gd
        ' ',  # Tb
        ' ',  # Dy
        ' ',  # Ho
        ' ',  # Er
        ' ',  # Tm
        ' ',  # Yb
        ' ',  # Lu
        ' ',  # Hf
        ' ',  # Ta
        ' ',  # W
        ' ',  # Re
        ' ',  # Os
        ' ',  # Ir
        ' ',  # Pt [1]
        ' ',  # Au [1]
        ' ',  # Hg [1]
        ' ',  # Tl [1]
        ' ',  # Pb [1]
        ' ',  # Bi [5]
        ' ',  # Po [5]
        ' ',  # At [5]
        ' ',  # Rn [5]
        ' ',  # Fr [5]
        ' ',  # Ra [5]
        ' ',  # Ac
        ' ',  # Th
        ' ',  # Pa
        ' ',  # U [1]
        ' ',  # Np
        ' ',  # Pu
        ' ',  # Am
        ' ',  # Cm
        ' ',  # Bk
        ' ',  # Cf
        ' ',  # Es
        ' ',  # Fm
        ' ',  # Md
        ' ',  # No
        ' '
    ]  # Lr
}
