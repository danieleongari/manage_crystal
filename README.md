# manage_crystal
A tool to convert crystal files (atoms coordinates + unit cell) into common files and extract some useful info

Usage: 

- to get default info about the crystal: 

```
$ manage_crystal.py inputfilename.inputformat [options]`
```

- to convert to another format:

```
$ manage_crystal.py inputfilename.inputformat [options] -o outputfilename.outputformat
```

or

```
$ manage_crystal.py inputfilename.inputformat [options] -o outputformat
```

- to get help:

```
$ manage_crystal.py --help 

usage: manage_crystal.py [-h] [-o OUTPUT] [-silent] [-show]
                         [-showonly SHOWONLY] [-cupw] [-void] [-ovlp]
                         [-pseudopw PSEUDOPW] [-bscp2k BSCP2K] [-resp RESP]
                         [-readcharge READCHARGE]
                         [-readrepeatcharge READREPEATCHARGE] [-x MULTIPL_X]
                         [-y MULTIPL_Y] [-z MULTIPL_Z] [-cutoff CUTOFF]
                         [-chargenull] [-printatoms] [-printatoms_noHCO]
                         [-transl [TRANSL [TRANSL ...]]] [-mol]
                         [-randomize RANDOMIZE] [-chkmetalcharge] [-chkcharge]
                         [-chkdef2] [-chkmepo] [-avgcharges]
                         [-normalizecharges] [-tm1] [-tm2] [-tm3] [-tm4]
                         [-tm5] [-tm6]
                         inputfile

Program to read, extract info and convert crystal files (by Daniele Ongari)

positional arguments:
  inputfile             path to the input file to read
                        IMPLEMENTED: xyz(w/CELL),pdb,cssr,pwi,pwo,cif,xsf,axsf,subsys(CP2K),
                                     restart(CP2K),inp(CP2K),cube,POSCAR(VASP) 
                                     [NEXT: gaussian, dcd+atoms,POSCAR(VASP)]

optional arguments:

  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output filename.extension or just the extension
                        IMPLEMENTED: cif,pdb,cssr,xyz(w/CELL),pwi,subsys(CP2K),axsf,geo(GULP)
  -silent               No output info on the screen
  -show                 Show all the info
                        [skip -silent]
  -showonly SHOWONLY    Show only the required info:
                        cell, CELL, xyz, fract, charge
                        [skip -silent]
  -cupw                 Look for a Copper PaddleWheel
  -void                 Compute void geometrically [NOT WORKING]
  -ovlp                 Look for an overlap and modify the file [WORK IN PROGRESS]
  -pseudopw PSEUDOPW    Pseudo for the .pwi output
  -bscp2k BSCP2K        Gaussian Basis Set for CP2K
  -resp RESP            Read the charges from a cp2k RESP file
                        (also checking if the atoms are the same)
                        BC1: it read the first set of charges
                        BC2: Also a cp2k output file with charges is fine!
  -readcharge READCHARGE
                        Read the charges from a simple list
  -readrepeatcharge READREPEATCHARGE
                        Read the charges from REPEAT output of QE
  -x MULTIPL_X          Extend in the x direction, by the specified times
  -y MULTIPL_Y          Extend in the y direction, by the specified times
  -z MULTIPL_Z          Extend in the z direction, by the specified times
  -cutoff CUTOFF        Automatically extend the UC so that the cutoff is respected
                        (TIP: use -cutoff 0 to just know the perpendicular widths!)
  -chargenull           Delete the charge of the atoms
  -printatoms           Print all atoms types
                        [skip -silent]
  -printatoms_noHCO     Print all atoms types exc. H,C,O
                        [skip -silent]
  -transl [TRANSL [TRANSL ...]]
                        x y z translation in Angs
  -mol                  Considers a molecule for xyz: no cell!
  -randomize RANDOMIZE  Randomize the geometry by a gaussian
                        with the specified delta (angs)
  -chkmetalcharge       Check if the charge on a metal (see list) is neg.
                        [skip -silent]
  -chkcharge            Check if all the charges are zero.
                        [skip -silent]
  -chkdef2              Check if there is a non def2 BS atom (H-La, Hf-Rn).
                        [skip -silent]
  -chkmepo              Check if there is a non MEPO atom (H,V,Cu,Zn,C,N,O,F,Cl,Br,I).
                        [skip -silent]
  -avgcharges           Use average charges from DDEC
  -normalizecharges     Normalize the charges to have a null total charge.
  -tm1                  Tailor-made 1: read  .cif CoRE MOF w/DDEC charges
  -tm2                  Tailor-made 2: print .cif for EQeq
  -tm3                  Tailor-made 3: read  .xyz for B.Wells Qeq
  -tm4                  Tailor-made 4: print .xyz for B.Wells Qeq
  -tm5                  Tailor-made 5: print .xyz for B.Wells Qeq w/zero FC
  -tm6                  Tailor-made 6: read GULP's cif
```
  
Tips:

- you may want to use `-silent` to suppress default verbose output: several options "skip -silent" so that you can print just that information on the screen (e.g., `-printatoms -silent` prints on the screen just the atom types on one line). This make easy to use bash loops for lists of structures. 
  
  
