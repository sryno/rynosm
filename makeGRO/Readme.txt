This script takes as input PDB or XYZ files and converts them in to Gromacs GRO files. It also take multiple structures from dihedral scans in XYZ format and create individual GRO files for each structure.


Example: python convertToGromacs.py dpp_trimer_1.xyz XYZ dpp_trimer -scan -NT dpp_atomtypes.txt -name DPP


usage: convertToGromacs.py [-h] [-scan] [--new] [-NT NT] [-name NAME]
                           INPUT_FILE INPUT_FILE GROMACS_FILE

Allows for the conversion of a PDB, XYZ, or G09 formatted file to Gromacs GRO
format.

positional arguments:
  INPUT_FILE    Input File Name and file type [PDB, XYZ, G09]
  GROMACS_FILE  Gromacs Gro Output File Prefix.

optional arguments:
  -h, --help    show this help message and exit
  -scan         Turns on Parsing XYZ files with multiple structures from
                scans.
  --new         Switch to turn on asking for new atom types.
  -NT NT        Extracts atom types from an external file. Required for XYZ
                and G09 files. [Filename]
  -name NAME    System Name. Only used for XYZ and PDB.



