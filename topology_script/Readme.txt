This is a script that takes as input a Gro structure file, non-bonded itp file, and bonded itp file. It outputs a new Gro structure file, forcefield files that contain all explicit connections (i.e., removes wildcards), and a position restraint file for heavy atoms. Missing bonded interations are printed to a Missing file.

The program will remove box dimensions in the output GRO file that will need to be added back by the user.

usage: makeITP.py [-h] -gro GRO -out OUT [-gout GOUT] -nb NB -bon BON [-ff]
                  [-imp] [-posre] [-gchg GCHG GCHG] [--debug]

Reads in a Gromacs .gro file for a molecule and outputs a molecular .itp file.
Options exist to create posre.itp and topol.top files.

optional arguments:
  -h, --help       show this help message and exit
  -gro GRO         Gromacs .gro Input File.
  -out OUT         Gromacs .itp Output File.
  -gout GOUT       Gromacs .gro Output file.
  -nb NB           Gromacs Non-Bonded .itp Input File.
  -bon BON         Gromacs Bonded .itp Input File.
  -ff              Enable rewriting new force field files that remove
                   wildcards.
  -imp             Turn on finding impropers defined in -bon file. DO NOT USE!
  -posre           Enable the Creation of a posre.itp File.
  -gchg GCHG GCHG  Turns on the extraction of charges from G09 files. Atom
                   order should be identical to GRO file. G09_output_file
                   charge_types[Mulliken, Hirschfeld, CM5, TXT]
  --debug          Enable Debug Information. Does not print out files.
  

  File formats:

  Input GRO file:

  SOL
   12
    1SOL   CHCA    1   0.770   0.927   0.782    ;    C
    1SOL   CHCA    2   0.655   0.847   0.782    ;    C
    1SOL   CHCA    3   0.667   0.707   0.782    ;    C
    1SOL   CHAA    4   0.794   0.648   0.782    ;    C
    1SOL   CHCA    5   0.908   0.728   0.782    ;    C
    1SOL   CHCA    6   0.896   0.867   0.782    ;    C
    1SOL   CHHA    7   0.760   1.035   0.782    ;    H
    1SOL   CHHA    8   0.557   0.893   0.782    ;    H
    1SOL   CHHA    9   1.007   0.683   0.782    ;    H
    1SOL   CHHA   10   0.985   0.930   0.782    ;    H
    1SOL   CHHA   11   0.578   0.646   0.782    ;    H
    1SOL   CHCl   12   0.809   0.473   0.782    ;   Cl



  Input NB file:

[ atomtypes ]
;name  bond_type    at.num       mass         charge        ptype       sigma[nm]         eps[kJ/mol]
 CHCA   CA             6       12.01100    -0.115       A    3.55000e-01  2.92880e-01
 CHHA   HA             1        1.00800     0.115       A    2.42000e-01  1.25520e-01
 CHAA   CA             6       12.01100     0.180       A    3.55000e-01  2.92880e-01
 CHCl   Cl            17       35.45300    -0.180       A    3.40000e-01  1.25520e+00



  Input BON file:

[ bondtypes ]
;  ai   aj    funct      b0[nm]     kb[kJ/mol nm^2]
   CA    HA      1    0.10800   307105.6
   CA    CA      1    0.14000   392459.2
   CA    Cl      1    0.17250   251040.0

[ angletypes ]
;  ai   aj   ak   funct    theta0     k0(kjmol-1 rad-2)
   CA     CA     CA      1   120.000    527.184
   CA     CA     HA      1   120.000    292.880
   CA     CA     Cl      1   120.000    627.600

[ dihedraltypes ]
; proper dihedrals
;  ai   aj   ak    al   funct    c0        c1        c2         c3        c4        c5(kj/mol) 
    X      CA     CA     X       3     30.33400   0.00000 -30.33400   0.00000   0.00000   0.00000

[ dihedraltypes ]
; Improper dihedrals
#define improper_Z_Y_CA_X            180.0        4.6024     2



  Output GRO file:

  SOL
   12    
    1SOL     CA    1   0.770   0.927   0.782
    1SOL     CA    2   0.655   0.847   0.782
    1SOL     CA    3   0.667   0.707   0.782
    1SOL     CA    4   0.794   0.648   0.782
    1SOL     CA    5   0.908   0.728   0.782
    1SOL     CA    6   0.896   0.867   0.782
    1SOL     HA    7   0.760   1.035   0.782
    1SOL     HA    8   0.557   0.893   0.782
    1SOL     HA    9   1.007   0.683   0.782
    1SOL     HA   10   0.985   0.930   0.782
    1SOL     HA   11   0.578   0.646   0.782
    1SOL     Cl   12   0.809   0.473   0.782
