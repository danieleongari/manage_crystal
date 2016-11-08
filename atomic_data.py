################################# atom name

atomic_symbol=('X','H','He',
            'Li','Be','B','C','N','O','F','Ne',
            'Na','Mg','Al','Si','P','S','Cl','Ar',
            'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
            'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
            'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo')

# from ASE
#chemical_symbols = ('X',  'H',  'He', 'Li', 'Be',
#                    'B',  'C',  'N',  'O',  'F',
#                    'Ne', 'Na', 'Mg', 'Al', 'Si',
#                    'P',  'S',  'Cl', 'Ar', 'K',
#                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
#                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
#                    'Zn', 'Ga', 'Ge', 'As', 'Se',
#                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
#                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
#                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
#                    'Sn', 'Sb', 'Te', 'I',  'Xe',
#                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
#                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
#                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
#                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
#                    'Re', 'Os', 'Ir', 'Pt', 'Au',
#                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
#                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
#                    'Th', 'Pa', 'U',  'Np', 'Pu',
#                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
#                    'Fm', 'Md', 'No', 'Lr')


atomic_name = (
    'DUMMY', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium')

################################# atom vdw from ASE
# Van der Waals radii in [A] taken from
#http://www.webelements.com/periodicity/van_der_waals_radius/
#and the references given there.
#Additional source 5 from http://de.wikipedia.org/wiki/Van-der-Waals-Radius
#
#1. A. Bondi, J. Phys. Chem., 1964, 68, 441.
#
#2. L. Pauling, The Nature of the Chemical Bond,
#   Cornell University Press, USA, 1945.
#
#3. J.E. Huheey, E.A. Keiter, and R.L. Keiter in Inorganic Chemistry
#   Principles of Structure and Reactivity, 4th edition, HarperCollins,
#   New York, USA, 1993.W.W. Porterfield in Inorganic chemistry,
#   a unified approach, Addison Wesley Publishing Co.,
#   Reading Massachusetts, USA, 1984.
#
#4. A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data,
#   Macmillan, London, UK, 1992.
#
#5. Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero,
#   Christopher J. Cramer, Donald G. Truhlar Consistent van der Waals Radii
#   for the whole main group in J.phys. Chem A. 2009, 113, 5806-5812,
#   doi:10.1021/jp8111556


atomic_vdw = (
    1000000,  # X
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

################################# atom mass from ASE
atomic_mass = (
   0.00000, # X
   1.00794, # H
   4.00260, # He
   6.94100, # Li
   9.01218, # Be
  10.81100, # B
  12.01100, # C
  14.00670, # N
  15.99940, # O
  18.99840, # F
  20.17970, # Ne
  22.98977, # Na
  24.30500, # Mg
  26.98154, # Al
  28.08550, # Si
  30.97376, # P
  32.06600, # S
  35.45270, # Cl
  39.94800, # Ar
  39.09830, # K
  40.07800, # Ca
  44.95590, # Sc
  47.88000, # Ti
  50.94150, # V
  51.99600, # Cr
  54.93800, # Mn
  55.84700, # Fe
  58.93320, # Co
  58.69340, # Ni
  63.54600, # Cu
  65.39000, # Zn
  69.72300, # Ga
  72.61000, # Ge
  74.92160, # As
  78.96000, # Se
  79.90400, # Br
  83.80000, # Kr
  85.46780, # Rb
  87.62000, # Sr
  88.90590, # Y
  91.22400, # Zr
  92.90640, # Nb
  95.94000, # Mo
    1000000, # Tc
 101.07000, # Ru
 102.90550, # Rh
 106.42000, # Pd
 107.86800, # Ag
 112.41000, # Cd
 114.82000, # In
 118.71000, # Sn
 121.75700, # Sb
 127.60000, # Te
 126.90450, # I
 131.29000, # Xe
 132.90540, # Cs
 137.33000, # Ba
 138.90550, # La
 140.12000, # Ce
 140.90770, # Pr
 144.24000, # Nd
    1000000, # Pm
 150.36000, # Sm
 151.96500, # Eu
 157.25000, # Gd
 158.92530, # Tb
 162.50000, # Dy
 164.93030, # Ho
 167.26000, # Er
 168.93420, # Tm
 173.04000, # Yb
 174.96700, # Lu
 178.49000, # Hf
 180.94790, # Ta
 183.85000, # W
 186.20700, # Re
 190.20000, # Os
 192.22000, # Ir
 195.08000, # Pt
 196.96650, # Au
 200.59000, # Hg
 204.38300, # Tl
 207.20000, # Pb
 208.98040, # Bi
    1000000, # Po
    1000000, # At
    1000000, # Rn
    1000000, # Fr
 226.02540, # Ra
    1000000, # Ac
 232.03810, # Th
 231.03590, # Pa
 238.02900, # U
 237.04820, # Np
    1000000, # Pu
    1000000, # Am
    1000000, # Cm
    1000000, # Bk
    1000000, # Cf
    1000000, # Es
    1000000, # Fm
    1000000, # Md
    1000000, # No
    1000000) # Lw
