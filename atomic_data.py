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

atomic_vdw_zeopp = (
 1,      #Dummy
 1.09,   #H  
 1.4,    #He 
 1.82,   #Li 
 2,      #Be 
 2,      #B  
 1.7,    #C  
 1.55,   #N  
 1.52,   #O  
 1.47,   #F  
 1.54,   #Ne 
 2.27,   #Na 
 1.73,   #Mg 
 2,      #Al 
 2.1,    #Si 
 1.8,    #P  
 1.8,    #S  
 1.75,   #Cl 
 1.88,   #Ar 
 2.75,   #K  
 2,      #Ca 
 2,      #Sc 
 2,      #Ti 
 2,      #V  
 2,      #Cr 
 2,      #Mn 
 2,      #Fe 
 2,      #Co 
 1.63,   #Ni 
 1.4,    #Cu 
 1.39,   #Zn 
 1.87,   #Ga 
 2,      #Ge 
 1.85,   #As 
 1.9,    #Se 
 1.85,   #Br 
 2.02,   #Kr 
 2,      #Rb 
 2,      #Sr 
 2,      #Y  
 2,      #Zr 
 2,      #Nb 
 2,      #Mo 
 2,      #Tc 
 2,      #Ru 
 2,      #Rh 
 1.63,   #Pd 
 1.72,   #Ag 
 1.58,   #Cd 
 1.93,   #In 
 2.17,   #Sn 
 2,      #Sb 
 2.06,   #Te 
 1.98,   #I  
 2.16,   #Xe 
 2,      #Cs 
 2,      #Ba 
 2,      #La 
 2,      #Ce 
 2,      #Pr 
 2,      #Nd 
 2,      #Pm 
 2,      #Sm 
 2,      #Eu 
 2,      #Gd 
 2,      #Tb 
 2,      #Dy 
 2,      #Ho 
 2,      #Er 
 2,      #Tm 
 2,      #Yb 
 2,      #Lu 
 2,      #Hf 
 2,      #Ta 
 2,      #W  
 2,      #Re 
 2,      #Os 
 2,      #Ir 
 1.72,   #Pt 
 1.66,   #Au 
 1.55,   #Hg 
 1.96,   #Tl 
 2.02,   #Pb 
 2,      #Bi 
 2,      #Po 
 2,      #At 
 2,      #Rn 
 2,      #Fr 
 2,      #Ra 
 2,      #Ac 
 2,      #Th 
 2,      #Pa 
 1.86,   #U  
 2,      #Np 
 2,      #Pu 
 2,      #Am 
 2,      #Cm 
 2,      #Bk 
 2,      #Cf 
 2,      #Es 
 2,      #Fm 
 2,      #Md 
 2,      #No 
 2,      #Lr 
 2,      #Rf 
 2,      #Db 
 2,      #Sg 
 2,      #Bh 
 2,      #Hs 
 2,      #Mt 
 2)      #Ds 

atomic_vdw_UFF = (
 1,      #Dummy
   1.286,   #H 
   1.052,  #He
   1.092,  #Li
   1.223,  #Be
   1.819,   #B 
   1.715,   #C 
   1.630,    #N 
   1.559,   #O 
   1.498,   #F 
   1.445,  #Ne
   1.329,  #Na
   1.346,  #Mg
   2.004,  #Al
   1.913,  #Si
   1.847,   #P 
   1.797,   #S 
   1.758,  #Cl
   1.723,  #Ar
   1.698,   #K 
   1.514,  #Ca
   1.468,  #Sc
   1.414,  #Ti
   1.400,   #V 
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
   1.490,   #Y 
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
   2.005,   #I 
   1.962,  #Xe
   2.012,  #Cs
   1.649,  #Ba
   1.569,  #La
   1.584,  #Ce
   1.606,  #Pr
   1.592,  #Nd
   1.580,   #Pm
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
   1.367,   #W 
   1.316,  #Re
   1.390,   #Os
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
   1.512,   #U 
   1.525,  #Np
   1.525,  #Pu
   1.506,  #Am
   1.482,  #Cm
   1.487,  #Bk
   1.476,  #Cf
   1.470,   #Es
   1.464,  #Fm
   1.458,  #Md
   1.447,  #No
   1.441,  #Lr
 2,      #Rf 
 2,      #Db 
 2,      #Sg 
 2,      #Bh 
 2,      #Hs 
 2,      #Mt 
 2)      #Ds 

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
    1000000) # Lr

################################# pseudoopotential: dictionary[pseudo_type][atomic_number]
atomic_pseudo = {'pbe': [   
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
    ' '], # Lr
  
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
    ' ']  # Lr
}
################################# average charges fromm DDEC CoRE MOF
atomic_ddecavgcharges = (
   0.000, #dummy
   0.116, # H
   0.000, # He
   0.890, # Li
   1.219, # Be
   0.354, # B
   0.083, # C
   -0.303, # N
   -0.619, # O
   -0.241, # F
   0.000, # Ne
   0.934, # Na
   1.539, # Mg
   1.912, # Al
   1.448, # Si
   1.609, # P
   0.336, # S
   0.096, # Cl
   0.000, # Ar
   0.929, # K
   1.580, # Ca
   2.066, # Sc
   1.843, # Ti
   1.638, # V
   1.183, # Cr
   1.083, # Mn
   0.963, # Fe
   0.953, # Co
   0.734, # Ni
   0.743, # Cu
   0.945, # Zn
   1.592, # Ga
   1.759, # Ge
   1.109, # As
   0.038, # Se
   -0.252, # Br
   0.000, # Kr
   0.932, # Rb
   1.600, # Sr
   2.047, # Y
   0.000, # Zr
   0.320, # Nb
   0.960, # Mo
   0.000, # Tc
   1.097, # Ru
   0.674, # Rh
   0.457, # Pd
   0.326, # Ag
   0.926, # Cd
   1.568, # In
   1.175, # Sn
   1.989, # Sb
   0.100, # Te
   -0.334, # I
   0.000, # Xe
   0.907, # Cs
   1.642, # Ba
   2.155, # La
   0.000, # Ce
   0.000, # Pr
   0.000, # Nd
   0.000, # Pm
   0.000, # Sm
   0.000, # Eu
   0.000, # Gd
   0.000, # Tb
   0.000, # Dy
   0.000, # Ho
   0.000, # Er
   0.000, # Tm
   0.000, # Yb
   0.000, # Lu
   2.488, # Hf
   0.000, # Ta
   0.913, # W
   -0.107, # Re
   0.000, # Os
   -0.070, # Ir
   0.015, # Pt
   0.017, # Au
   0.439, # Hg
   0.000, # Tl
   0.955, # Pb
   1.216, # Bi
   0.000, # Po
   0.000, # At
   0.000) # Rn


