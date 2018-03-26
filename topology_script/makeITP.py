from __future__ import print_function

import argparse
import sys
from scipy import *
from scipy.sparse import *

__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = 'Sean M. Ryno'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def getG09Charges(g09File, chgType):
    fin = open(g09File, 'r')
    fileLine = []
    gaussCharges = []
    for line in fin:
        fileLine.append(line)
    if chgType in ['Hirschfeld', 'CM5']:
        for line in fileLine:
            if 'Hirshfeld charges, spin densities, dipoles, and CM5 charges' in line:
                chargeIndex = fileLine.index(line)
                break
            else:
                pass
    elif chgType == 'Mulliken':
        for line in fileLine:
            if 'Mulliken charges:' in line:
                chargeIndex = fileLine.index(line)
                break
            else:
                pass
    elif chgType == 'TXT':
        chargeIndex = 0
    else:
        print("There is some error in getting G09 charges. (1)", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)
    chargeIndex += 2
    if chgType == 'Hirschfeld':
        for i in range(chargeIndex, len(fileLine)):
            line = fileLine[i].split()
            if len(line) == 8:
                gaussCharges.append(float(line[2]))
            elif line[0] == 'Tot':
                break
            else:
                break
    elif chgType == 'CM5':
        for i in range(chargeIndex, len(fileLine)):
            line = fileLine[i].split()
            if len(line) == 8:
                gaussCharges.append(float(line[7]))
            elif line[0] == 'Tot':
                break
            else:
                break
    elif chgType == 'Mulliken':
        for i in range(chargeIndex, len(fileLine)):
            line = fileLine[i].split()
            if len(line) == 3:
                gaussCharges.append(float(line[2]))
            else:
                break
    elif chgType == 'TXT':
        for i in range(len(fileLine)):
            line = fileLine[i].split()
            if len(line) != 0:
                gaussCharges.append(float(line[0]))
            else:
                pass
    else:
        print("There is some error in getting G09 charges. (2)", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    return gaussCharges


def parseGro(groFile):
    fin = open(groFile, 'r')
    title = fin.readline()
    numAtoms = int(fin.readline().strip())
    resNum, name, atomType, atomNum, x, y, z = [], [], [], [], [], [], []
    for line in fin:
        if len(line.split()) > 3:
            resNum.append(int(line[0:5].strip()))
            name.append(line[5:10].strip())
            atomType.append(line[10:15].strip())
            atomNum.append(int(line[15:20].strip()))
            x.append(float(line[20:28].strip()))
            y.append(float(line[28:36].strip()))
            z.append(float(line[36:44].strip()))
            # elements.append(line[52:].strip())

    return resNum, name, atomType, atomNum, x, y, z


def parseFF(itpFile):
    itp = {}
    fin = open(itpFile, 'r')
    for line in fin:
        if 'atomtypes' in line:
            pass
        elif (line[0] == ';') or (line.strip() == ''):
            pass
        else:
            line = line.split()
            itp[line[0]] = [line[1], int(line[2]), float(line[3]), float(line[4]), line[5], float(line[6]),
                            float(line[7])]
            # itp[line[1]] = [line[0], float(line[2]), float(line[3]), line[4], float(line[5]), float(line[6])]

    return itp


def parseBon(itpFile):
    bonds = []
    angles = []
    dihedrals = []
    proper_dihedrals = []
    impropers = []
    fin = open(itpFile, 'r')
    for line in fin:
        if (line.strip() == '') or (line.strip()[0] == '[') or (line[0] == ';'):
            pass
        elif '#define' in line:
            line = line.split()
            impropers.append([line[0], line[1], line[2], line[3], line[4]])
        elif (len(line.split()) == 5) or ((len(line.split()) >= 6) and (line.split()[5] == ';')):
            line = line.split()
            bonds.append([line[0], line[1], line[2], line[3], line[4]])
        elif (len(line.split()) == 6) or ((len(line.split()) >= 7) and (line.split()[6] == ';')):
            line = line.split()
            angles.append([line[0], line[1], line[2], line[3], line[4], line[5]])
        elif ((len(line.split()) == 11) and (line.split()[4] == '3')) or \
                ((len(line.split()) > 11) and (line.split()[11] == ';')):
            line = line.split()
            dihedrals.append([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9],
                              line[10]])
        elif ((len(line.split()) == 8) and (int(line.split()[4]) == 9)) or ((len(line.split()) > 8) and (line.split()[8] == ';') and (int(line.split()[4]) == 9)):
            line = line.split()
            proper_dihedrals.append([line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]])
        else:
            pass

    return bonds, angles, dihedrals, proper_dihedrals, impropers


def findElements(atomTypes, nb):
    elementNums = []
    radii = []
    elements = []
    for i in atomTypes:
        elementNums.append(nb[i][1])
    for i in range(len(elementNums)):
        elements.append(assign_element(elementNums[i]))
        elementNums[i] = assign_element(elementNums[i])
    for i in elementNums:
        radii.append(assign_radii(i))

    return radii, elements


def assign_element(atomNum):
    """
    Determine element based on atom number.
    """
    if atomNum == 1:
        element = 'H'
    elif atomNum == 6:
        element = 'C'
    elif atomNum == 7:
        element = 'N'
    elif atomNum == 8:
        element = 'O'
    elif atomNum == 9:
        element = 'F'
    elif atomNum == 16:
        element = 'S'
    elif atomNum == 14:
        element = 'Si'
    elif atomNum == 17:
        element = 'Cl'
    else:
        print("Error in assigning element", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    return element


def assign_radii(atomType):
    """
    Determines the vdW radii of atoms based on Element type
    """
    if atomType == 'X':
        atomRadius = 0.023
    elif atomType == 'H':
        atomRadius = 0.023
    elif atomType == 'C':
        atomRadius = 0.068
    elif atomType == 'F':
        atomRadius = 0.064
    elif atomType == 'Si':
        atomRadius = 0.080
    elif atomType == 'S':
        atomRadius = 0.068
    elif atomType == 'N':
        atomRadius = 0.065
    elif atomType == 'O':
        atomRadius = 0.060
    elif atomType == 'Cl':
        atomRadius = 0.099
    else:
        print("Error in assigning radius", file=sys.stderr)
        print("Offending atomType: ", atomType, file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    return atomRadius


def AtomConnections(natoms, x, y, z, atom_radii):
    """
        Creates sparse matrix that stores atom connections
        Atoms are defined to be connected to themselves
    """
    connected = lil_matrix((natoms, natoms))

    for i in range(natoms):
        for j in range(i):
            temp_distance = sqrt(
                (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]))
            if temp_distance <= (1.5 * (atom_radii[i] + atom_radii[j])):
                connected[i, j] = 1
                connected[j, i] = 1
            else:
                pass
        connected[i, i] = 0

    return connected


def findBonds(connections):
    size = connections[0].shape[1]
    bonds = []
    for i in range(size):
        for j in range(size):
            if (i == j):
                pass
            elif connections[i, j] == 1:
                if ([i, j] in bonds) or ([j, i] in bonds):
                    pass
                else:
                    bonds.append([i, j])

    return bonds


def findAnglesOrig(connections):
    size = connections[0].shape[1]
    angles = []
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if (i == j) or (i == k) or (j == k):
                    pass
                elif (connections[i, j] == 1):
                    if (connections[i, k] == 1) or (connections[j, k] == 1):
                        if ([i, j, k] in angles) or ([k, i, j] in angles) or ([j, i, k] in angles) or (
                                    [k, j, i] in angles):
                            pass
                        else:
                            angles.append([i, j, k])

    return angles


def findAngles(connections, bonds):
    size = connections[0].shape[1]
    angles = []
    for i in bonds:
        for j in range(size):
            if (i[0] == j) or (i[1] == j):
                pass
            elif (connections[i[0], j] == 1):
                if ([j, i[0], i[1]] in angles) or ([j, i[1], i[0]] in angles):
                    pass
                else:
                    angles.append([j, i[0], i[1]])
            elif (connections[i[1], j] == 1):
                if ([i[0], i[1], j] in angles) or ([i[1], i[0], j] in angles):
                    pass
                else:
                    angles.append([i[0], i[1], j])

    return angles


def findDihedralsOrig(connections):
    size = connections[0].shape[1]
    dihedrals = []
    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
                    if (i == j) or (i == k) or (i == l) or (j == k) or (j == l) or (k == l):
                        pass
                    elif ([k, i, j, l] in dihedrals) or ([l, i, j, k] in dihedrals) or ([k, j, i, l] in dihedrals) or (
                                [l, j, i, k] in dihedrals):
                        pass
                    elif (connections[i, j] == 1):
                        if (connections[i, k] == 1) or (connections[i, l] == 1):
                            if (connections[j, k] == 1) or (connections[j, l] == 1):
                                dihedrals.append([k, i, j, l])

    return dihedrals


def findDihedrals(connections, angles):
    size = connections[0].shape[1]
    dihedrals = []
    for i in angles:
        for j in range(size):
            if (i[0] == j) or (i[1] == j) or (i[2] == j):
                pass
            elif (connections[i[0], j] == 1):
                if ([i[2], i[1], i[0], j] in dihedrals) or ([j, i[0], i[1], i[2]] in dihedrals):
                    pass
                else:
                    dihedrals.append([j, i[0], i[1], i[2]])
            elif (connections[i[2], j] == 1):
                if ([i[0], i[1], i[2], j] in dihedrals) or ([j, i[2], i[1], i[0]] in dihedrals):
                    pass
                else:
                    dihedrals.append([i[0], i[1], i[2], j])

    return dihedrals


def findImpropers(nb, improperParams, elements, atomType, bonds):
    improperSites = []
    improperParamList = []
    impropBonds = []
    impropDihedrals = []
    for i in range(len(improperParams)):
        for j in range(len(atomType)):
            # print(nb[atomType[j]][0])
            if improperParams[i][1].split('_')[3] == nb[atomType[j]][0]:
                tempImproper = improperParams[i][1].split('_')
                improperSites.append([tempImproper[1], tempImproper[2], j, tempImproper[4]])
                improperParamList.append(improperParams[i])
    # print(improperParamList)
    # print(improperParams)
    for i in range(len(improperSites)):
        impropBonds.append([])
        for j in range(len(bonds)):
            if improperSites[i][2] in bonds[j]:
                impropBonds[i].append(bonds[j])
    for i in range(len(impropBonds)):
        for j in range(len(impropBonds[i])):
            list = []
            if (elements[impropBonds[i][j][0]] == 'H') or (elements[impropBonds[i][j][1]] == 'H'):
                list.append(j)
            else:
                pass
            if len(list) > 0:
                for k in sorted(list, reverse=True):
                    impropBonds[i].pop(k)
    # print(impropBonds)
    for i in range(len(impropBonds)):
        if len(impropBonds[i]) == 3:
            # print(impropBonds[i])
            # print(impropBonds[i][0])
            # print(impropBonds[i][0][0])
            a = improperSites[i][2]
            if impropBonds[i][0][0] != a:
                b = impropBonds[i][0][0]
            elif impropBonds[i][0][1] != a:
                b = impropBonds[i][0][1]
            else:
                print("These is something wrong with impropers. (1)", file=sys.stderr)
                sys.exit(1)
            if impropBonds[i][1][0] != a:
                c = impropBonds[i][1][0]
            elif impropBonds[i][1][1] != a:
                c = impropBonds[i][1][1]
            else:
                print("These is something wrong with impropers. (2)", file=sys.stderr)
                sys.exit(1)
            if impropBonds[i][2][0] != a:
                d = impropBonds[i][2][0]
            elif impropBonds[i][2][1] != a:
                d = impropBonds[i][2][1]
            else:
                print("These is something wrong with impropers. (3)", file=sys.stderr)
                sys.exit(1)
            # impropDihedrals.append([nb[atomType[b]][0],nb[atomType[c]][0],nb[atomType[a]][0],nb[atomType[d]][0]])
            impropDihedrals.append([b, c, a, d])

    # for i in range(len)
    # print(impropDihedrals)
    # print(improperParamList)
    # print(impropBonds)
    # print(improperSites)

    return impropDihedrals, improperParamList


def assignTypes(nb, atomType, bonds, angles, dihedral):
    bondTypes = []
    angleTypes = []
    dihedralTypes = []
    for i in range(len(bonds)):
        bondTypes.append([nb[atomType[bonds[i][0]]][0], nb[atomType[bonds[i][1]]][0]])
    for i in range(len(angles)):
        angleTypes.append([nb[atomType[angles[i][0]]][0], nb[atomType[angles[i][1]]][0], nb[atomType[angles[i][2]]][0]])
    for i in range(len(dihedral)):
        dihedralTypes.append(
            [nb[atomType[dihedral[i][0]]][0], nb[atomType[dihedral[i][1]]][0], nb[atomType[dihedral[i][2]]][0],
             nb[atomType[dihedral[i][3]]][0]])

    return bondTypes, angleTypes, dihedralTypes


def printTopol(outFile, resNum, name, atomType, atomNum, nb, bonds, angles, dihedrals, bondParams, angleParams,
               properDihedralParams, dihedralParams, improperPrint, gaussCharges):
    fout = open(outFile, 'w')
    mout = open(outFile + '.missing', 'w')

    # Append Proper Dihedrals to RB Dihedrals
    allDihedralParams = []
    for i in properDihedralParams:
        allDihedralParams.append(i)
    for i in dihedralParams:
        allDihedralParams.append(i)

    # Assign Atom Types
    bondTypes, angleTypes, dihedralTypes = assignTypes(nb, atomType, bonds, angles, dihedrals)
    missingBonds, missingAngles, missingDihedrals = [], [], []

    # Begin Writing Topology File
    print('; Topology file for ', outFile, file=fout)
    print('', file=fout)
    print('[ moleculetype ]', file=fout)
    print('; name     nrexcl', file=fout)
    print('  ', name[0], '     3', file=fout)
    print('', file=fout)

    # Print Atoms
    print('[ atoms ]', file=fout)
    print(';   nr     type  resnr residue        atom     cgnr    charge     mass', file=fout)
    if len(gaussCharges) == 0:
        for i in range(len(atomNum)):
            print(
                '{0:>6d}{1:>8s}{2:>6d}{3:>7s}{4:>15s}{5:>9d}{6:>14.6f}{7:>10.3f}'.format(atomNum[i], atomType[i],
                                                                                         resNum[i],
                                                                                         name[i], nb[atomType[i]][0],
                                                                                         atomNum[i], nb[atomType[i]][3],
                                                                                         nb[atomType[i]][2]), file=fout)
    else:
        for i in range(len(atomNum)):
            print(
                '{0:>6d}{1:>8s}{2:>6d}{3:>7s}{4:>15s}{5:>9d}{6:>14.6f}{7:>10.3f}'.format(atomNum[i], atomType[i],
                                                                                         resNum[i],
                                                                                         name[i], nb[atomType[i]][0],
                                                                                         atomNum[i], gaussCharges[i],
                                                                                         nb[atomType[i]][2]), file=fout)
    print('', file=fout)

    # Print Bonds
    print('[ bonds ]', file=fout)
    print(';   ai    aj   funct         c0            c1', file=fout)
    for i in range(len(bonds)):
        missing = 1
        for j in bondParams:
            if ([bondTypes[i][0], bondTypes[i][1]] == [j[0], j[1]]) or (
                        [bondTypes[i][1], bondTypes[i][0]] == [j[0], j[1]]):
                missing = 0
                break
            elif ((j[0] == 'X') and (bondTypes[i][1] == j[1])) or ((j[1] == 'X') and (bondTypes[i][1] == j[0])):
                missing = 0
                break
            elif ((j[0] == 'X') and (bondTypes[i][0] == j[1])) or ((j[1] == 'X') and (bondTypes[i][0] == j[0])):
                missing = 0
                break
            else:
                pass
        if missing == 0:
            print('{0:>6d}{1:>6d}{2:>6d}{3:>15.4f}{4:>15.1f}  ; {5:>6s}{6:>6s}'.format(bonds[i][0] + 1, bonds[i][1] + 1,
                                                                                       int(j[2]), float(j[3]),
                                                                                       float(j[4]), bondTypes[i][0],
                                                                                       bondTypes[i][1]), file=fout)
        elif missing == 1:
            if ([bondTypes[i][0], bondTypes[i][1]] not in missingBonds) and (
                        [bondTypes[i][1], bondTypes[i][0]] not in missingBonds):
                missingBonds.append([bondTypes[i][0], bondTypes[i][1]])
                print("Missing Bonds: {0:d} {1:d}".format(bonds[i][0] + 1, bonds[i][1] + 1))
            else:
                pass
        else:
            pass
    for i in missingBonds:
        print('Missing Bond Parameters: ', i[0], i[1], file=mout)
    print('', file=fout)

    # Print Angles
    print('[ angles ]', file=fout)
    print(';   ai    aj    ak   funct       theta0          k0', file=fout)
    for i in range(len(angles)):
        missing = 1
        for j in angleParams:
            if ([angleTypes[i][0], angleTypes[i][1], angleTypes[i][2]] == [j[0], j[1], j[2]]) or \
                    ([angleTypes[i][2], angleTypes[i][1], angleTypes[i][0]] == [j[0], j[1], j[2]]):
                missing = 0
                break
            elif ((j[0] == 'X') and ([angleTypes[i][1], angleTypes[i][2]] == [j[1], j[2]])) or \
                    ((j[2] == 'X') and ([angleTypes[i][0], angleTypes[i][1]] == [j[0], j[1]])) or \
                    ((j[0] == 'X') and (j[2] == 'X') and (angleTypes[i][1] == j[1])):
                missing = 0
                break
            elif ((j[2] == 'X') and ([angleTypes[i][1], angleTypes[i][0]] == [j[1], j[0]])) or \
                    ((j[0] == 'X') and ([angleTypes[i][2], angleTypes[i][1]] == [j[0], j[1]])) or \
                    ((j[0] == 'X') and (j[2] == 'X') and (angleTypes[i][1] == j[1])):
                missing = 0
                break
            elif ((j[1] == 'X') and ([angleTypes[i][0], angleTypes[i][2]] == [j[0], j[2]])) or \
                    ((j[1] == 'X') and ([angleTypes[i][0], angleTypes[i][2]] == [j[2], j[0]])):
                missing = 0
                break
            else:
                pass
        if missing == 0:
            print('{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>15.2f}{5:>15.3f}  ; {6:>6s}{7:>6s}{8:>6s}'.format(angles[i][0] + 1,
                                                                                                     angles[i][1] + 1,
                                                                                                     angles[i][2] + 1,
                                                                                                     int(j[3]),
                                                                                                     float(j[4]),
                                                                                                     float(j[5]),
                                                                                                     angleTypes[i][0],
                                                                                                     angleTypes[i][1],
                                                                                                     angleTypes[i][2]),
                  file=fout)
        elif missing == 1:
            if ([angleTypes[i][0], angleTypes[i][1], angleTypes[i][2]] not in missingAngles) and (
                        [angleTypes[i][2], angleTypes[i][1], angleTypes[i][0]] not in missingAngles):
                missingAngles.append([angleTypes[i][0], angleTypes[i][1], angleTypes[i][2]])
                print("Missing Angles: {0:d} {1:d} {2:d}".format(angles[i][0] + 1, angles[i][1] + 1, angles[i][2] + 1))
            else:
                pass
        else:
            pass
    for i in missingAngles:
        print('Missing Angle Parameters: ', i[0], i[1], i[2], file=mout)
    print('', file=fout)

    # Print Dihedrals
    print('[ dihedrals ]', file=fout)
    print(
        ';   ai    aj    ak    al   funct           c0             c1             c2             c3             c4             c5',
        file=fout)
    for i in range(len(dihedrals)):
        missing = 1
        for j in allDihedralParams:
            if ([dihedralTypes[i][0], dihedralTypes[i][1], dihedralTypes[i][2], dihedralTypes[i][3]] == [j[0], j[1],
                                                                                                         j[2], j[3]]) \
                    or ([dihedralTypes[i][0], dihedralTypes[i][1], dihedralTypes[i][2], dihedralTypes[i][3]] == [j[3],
                                                                                                                 j[2],
                                                                                                                 j[1],
                                                                                                                 j[0]]):
                missing = 0
                break
            elif ((j[0] == 'X') and (
                        [dihedralTypes[i][1], dihedralTypes[i][2], dihedralTypes[i][3]] == [j[1], j[2], j[3]])) \
                    or ((j[0] == 'X') and (
                                [dihedralTypes[i][2], dihedralTypes[i][1], dihedralTypes[i][0]] == [j[1], j[2], j[3]])):
                missing = 0
                break
            elif ((j[3] == 'X') and (
                        [dihedralTypes[i][0], dihedralTypes[i][1], dihedralTypes[i][2]] == [j[0], j[1], j[2]])) \
                    or ((j[3] == 'X') and (
                                [dihedralTypes[i][3], dihedralTypes[i][2], dihedralTypes[i][1]] == [j[0], j[1], j[2]])):
                missing = 0
                break
            elif ((j[0] == 'X') and (j[3] == 'X') and ([dihedralTypes[i][1], dihedralTypes[i][2]] == [j[1], j[2]])) \
                    or ((j[0] == 'X') and (j[3] == 'X') and (
                                [dihedralTypes[i][1], dihedralTypes[i][2]] == [j[2], j[1]])):
                missing = 0
                break
            else:
                pass
        if missing == 0:
            if len(j) == 11:
                print('{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>6d}{5:>15.6f}{6:>15.6f}{7:>15.6f}{8:>15.6f}{9:>15.6f}{10:>15.6f}  ; {11:>6s}{12:>6s}{13:>6s}{14:>6s}'.format(
                        dihedrals[i][0] + 1, dihedrals[i][1] + 1, dihedrals[i][2] + 1, dihedrals[i][3] + 1, int(j[4]),
                        float(j[5]), float(j[6]), float(j[7]), float(j[8]), float(j[9]), float(j[10]), dihedralTypes[i][0],
                        dihedralTypes[i][1], dihedralTypes[i][2], dihedralTypes[i][3]), file=fout)
            elif len(j) == 8:
                print('{0:>6d}{1:>6d}{2:>6d}{3:>6d}{4:>6d}{5:>15.6f}{6:>15.6f}{7:>15d}  ; {8:>6s}{9:>6s}{10:>6s}{11:>6s}'.format(
                        dihedrals[i][0] + 1, dihedrals[i][1] + 1, dihedrals[i][2] + 1, dihedrals[i][3] + 1, int(j[4]),
                        float(j[5]), float(j[6]), int(j[7]), dihedralTypes[i][0], dihedralTypes[i][1],
                        dihedralTypes[i][2], dihedralTypes[i][3]), file=fout)
        elif missing == 1:
            if ([dihedralTypes[i][0], dihedralTypes[i][1], dihedralTypes[i][2],
                 dihedralTypes[i][3]] not in missingDihedrals) and (
                        [dihedralTypes[i][3], dihedralTypes[i][2], dihedralTypes[i][1],
                         dihedralTypes[i][0]] not in missingDihedrals):
                missingDihedrals.append(
                    [dihedralTypes[i][0], dihedralTypes[i][1], dihedralTypes[i][2], dihedralTypes[i][3]])
                print("Missing dihedral: {0:d} {1:d} {2:d} {3:d}".format(dihedrals[i][0] + 1, dihedrals[i][1] + 1,
                                                                         dihedrals[i][2] + 1, dihedrals[i][3] + 1))
            else:
                pass
        else:
            pass
    for i in missingDihedrals:
        print('Missing Dihedral Parameters: ', i[0], i[1], i[2], i[3], file=mout)
    print('', file=fout)

    # Print Impropers
    if improperPrint == True:
        impropDihedrals, improperParamList = findImpropers(nb, improperParams, elements, atomType, bonds)
        print("[ dihedrals ]", file=fout)
        print(";  ai    aj    ak    al    funct            c0", file=fout)
        for i in range(len(impropDihedrals)):
            print(
                "{0:>6d}{1:>6d}{2:>6d}{3:>6d}     1  {4:s}".format(impropDihedrals[i][0] + 1, impropDihedrals[i][1] + 1,
                                                                   impropDihedrals[i][2] + 1, impropDihedrals[i][3] + 1,
                                                                   improperParamList[i][1]), file=fout)
        print("", file=fout)
    else:
        pass

    print('; Include Position restraint file', '#ifdef POSRES', '#include "posre.itp"', '#endif', sep='\n',
          file=fout)

    return bondTypes, angleTypes, dihedralTypes


def printTopol_noparams(outFile, resNum, name, atomType, atomNum, nb, bonds, angles, dihedrals, gaussCharges):
    fout = open(outFile, 'w')

    # Begin Writing Topology File
    print('; Topology file for ', outFile, file=fout)
    print('', file=fout)
    print('[ moleculetype ]', file=fout)
    print('; name     nrexcl', file=fout)
    print('  ', name[0], '     3', file=fout)
    print('', file=fout)

    # Print Atoms
    print('[ atoms ]', file=fout)
    print(';   nr     type  resnr residue        atom     cgnr    charge     mass', file=fout)
    if len(gaussCharges) == 0:
        for i in range(len(atomNum)):
            print(
                '{0:>6d}{1:>8s}{2:>6d}{3:>7s}{4:>15s}{5:>9d}{6:>14.6f}{7:>10.3f}'.format(atomNum[i], atomType[i],
                                                                                         resNum[i],
                                                                                         name[i], nb[atomType[i]][0],
                                                                                         atomNum[i], nb[atomType[i]][3],
                                                                                         nb[atomType[i]][2]), file=fout)
    else:
        for i in range(len(atomNum)):
            print(
                '{0:>6d}{1:>8s}{2:>6d}{3:>7s}{4:>15s}{5:>9d}{6:>14.6f}{7:>10.3f}'.format(atomNum[i], atomType[i],
                                                                                         resNum[i],
                                                                                         name[i], nb[atomType[i]][0],
                                                                                         atomNum[i], gaussCharges[i],
                                                                                         nb[atomType[i]][2]), file=fout)
    print('', file=fout)

    # Print Bonds
    print('[ bonds ]', file=fout)
    print(';   ai    aj   funct         c0            c1', file=fout)
    for i in range(len(bonds)):
        print('{0:>6d}{1:>6d}'.format(bonds[i][0] + 1, bonds[i][1] + 1), file=fout)

    # Print Angles
    print('[ angles ]', file=fout)
    print(';   ai    aj    ak   funct       theta0          k0', file=fout)
    for i in range(len(angles)):
        print('{0:>6d}{1:>6d}{2:>6d}     1'.format(angles[i][0] + 1, angles[i][1] + 1, angles[i][2] + 1), file=fout)

    # Print Dihedrals
    print('[ dihedrals ]', file=fout)
    print(
        ';   ai    aj    ak    al   funct           c0             c1             c2             c3             c4             c5',
        file=fout)
    for i in range(len(dihedrals)):
        print(
                '{0:>6d}{1:>6d}{2:>6d}{3:>6d}     3'.format(
                    dihedrals[i][0] + 1, dihedrals[i][1] + 1, dihedrals[i][2] + 1, dihedrals[i][3] + 1), file=fout)

    print('; Include Position restraint file', '#ifdef POSRES', '#include "posre.itp"', '#endif', sep='\n',
          file=fout)


def printGro(grOut, resNum, name, atomType, atomNum, x, y, z, nb):
    gout = open(grOut, 'w')
    print(name[0], file=gout)
    print("   {0:<6d}".format(len(resNum)), file=gout)
    for i in range(len(atomType)):
        print("{0:>5d}{1:<5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}".format(resNum[i], name[i], nb[atomType[i]][0],
                                                                           atomNum[i], x[i], y[i], z[i]), file=gout)


def printPosRes(atomNum, elements):
    fout = open('posre.itp', 'w')
    print('; In this topology include file, you will find position restraint',
          '; entries for all heavy atoms in your original .gro file.',
          '; This means that all protons are not restrained.', '', '[ position_restraints ]',
          '; atom  funct  fx    fy    fz', sep='\n', file=fout)
    for i in range(len(atomNum)):
        if elements[i] != 'H':
            print('{0:>6d}     1  1000  1000  1000'.format(atomNum[i]), file=fout)
        else:
            pass


def printFF(nb, bondParams, angleParams, properDihedralParams, dihedralParams, improperParams, atomType, bonds, angles,
            dihedrals):
    bondTypes, angleTypes, dihedralTypes = assignTypes(nb, atomType, bonds, angles, dihedrals)
    newBonds = []
    newAngles = []
    newDihedrals = []
    for i in bondTypes:
        if ([i[0], i[1]] in newBonds) or ([i[1], i[0]] in newBonds):
            pass
        else:
            newBonds.append([i[0], i[1]])
    for i in angleTypes:
        if ([i[0], i[1], i[2]] in newAngles) or ([i[2], i[1], i[0]] in newAngles):
            pass
        else:
            newAngles.append([i[0], i[1], i[2]])
    for i in dihedralTypes:
        if ([i[0], i[1], i[2], i[3]] in newDihedrals) or ([i[3], i[2], i[1], i[0]] in newDihedrals):
            pass
        else:
            newDihedrals.append([i[0], i[1], i[2], i[3]])

    # Append Proper Dihedrals to RB Dihedrals
    allDihedralParams = []
    for i in properDihedralParams:
        allDihedralParams.append(i)
    for i in dihedralParams:
        allDihedralParams.append(i)
    # print(allDihedralParams)

    nbout = open('ffnb.itp', 'w')
    bonout = open('ffbon.itp', 'w')

    print("; Generated with makeITP.py by Sean M. Ryon", file=nbout)
    print("[ atomtypes ]", file=nbout)
    print(
        ";       name     bondType   at.num       mass         charge        ptype       sigma[nm]         eps[kJ/mol]",
        file=nbout)
    for i in sorted(nb):
        print("{0:>12s}{1:>10s}{2:>10d}{3:>14.3f}{4:>14.6f}{5:>11s}{6:>16.4f}{7:>20.6f}".format(i, nb[i][0], nb[i][1],
                                                                                                nb[i][2], nb[i][3],
                                                                                                nb[i][4], nb[i][5],
                                                                                                nb[i][6]), file=nbout)
    print('', file=nbout)

    print("; Generated with makeITP.py by Sean M. Ryon", file=bonout)
    print("[ bondtypes ]", file=bonout)
    print(";  ai    aj    funct      b0[nm]     kb[kJ/mol nm^2]", file=bonout)
    for i in newBonds:
        found = False
        for j in bondParams:
            if ([i[0], i[1]] == [j[0], j[1]]) or ([i[1], i[0]] == [j[0], j[1]]):
                print(
                    "{0:>6s}{1:>6s}{2:>6d}{3:>14.4f}{4:>16.2f}".format(i[0], i[1], int(j[2]), float(j[3]), float(j[4])),
                    file=bonout)
                found = True
                break
            elif ((j[0] == 'X') and (i[1] == j[1])) or ((j[1] == 'X') and (i[1] == j[0])):
                print(
                    "{0:>6s}{1:>6s}{2:>6d}{3:>14.4f}{4:>16.2f}".format(i[0], i[1], int(j[2]), float(j[3]), float(j[4])),
                    file=bonout)
                found = True
                break
            elif ((j[0] == 'X') and (i[0] == j[1])) or ((j[1] == 'X') and (i[0] == j[0])):
                print(
                    "{0:>6s}{1:>6s}{2:>6d}{3:>14.4f}{4:>16.2f}".format(i[0], i[1], int(j[2]), float(j[3]), float(j[4])),
                    file=bonout)
                found = True
                break
            else:
                found = False
        if found == True:
            pass
        elif found == False:
            print("Match not found:", i)
    print('', file=bonout)

    print("[ angletypes ]", file=bonout)
    print(";  ai    aj    ak   funct    theta0     k0(kjmol-1 rad-2)", file=bonout)
    for i in newAngles:
        found = False
        for j in angleParams:
            if ([i[0], i[1], i[2]] == [j[0], j[1], j[2]]) or ([i[2], i[1], i[0]] == [j[0], j[1], j[2]]):
                print(
                    "{0:>6s}{1:>6s}{2:>6s}{3:>6d}{4:>12.2f}{5:>16.3f}".format(i[0], i[1], i[2], int(j[3]), float(j[4]),
                                                                              float(j[5])),
                    file=bonout)
                found = True
                break
            elif ((j[0] == 'X') and ([i[1], i[2]] == [j[1], j[2]])) or \
                    ((j[2] == 'X') and ([i[0], i[1]] == [j[0], j[1]])) or \
                    ((j[0] == 'X') and (j[2] == 'X') and (i[1] == j[1])):
                print(
                    "{0:>6s}{1:>6s}{2:>6s}{3:>6d}{4:>12.2f}{5:>16.3f}".format(i[0], i[1], i[2], int(j[3]), float(j[4]),
                                                                              float(j[5])),
                    file=bonout)
                found = True
                break
            elif ((j[2] == 'X') and ([i[1], i[0]] == [j[1], j[0]])) or \
                    ((j[0] == 'X') and ([i[2], i[1]] == [j[0], j[1]])) or \
                    ((j[0] == 'X') and (j[2] == 'X') and (i[1] == j[1])):
                print(
                    "{0:>6s}{1:>6s}{2:>6s}{3:>6d}{4:>12.2f}{5:>16.3f}".format(i[0], i[1], i[2], int(j[3]), float(j[4]),
                                                                              float(j[5])),
                    file=bonout)
                found = True
                break
            elif ((j[1] == 'X') and ([i[0], i[2]] == [j[0], j[2]])) or \
                    ((j[1] == 'X') and ([i[0], i[2]] == [j[2], j[0]])):
                print(
                    "{0:>6s}{1:>6s}{2:>6s}{3:>6d}{4:>12.2f}{5:>16.3f}".format(i[0], i[1], i[2], int(j[3]), float(j[4]),
                                                                              float(j[5])),
                    file=bonout)
                found = True
                break
            else:
                found = False
        if found == True:
            pass
        elif found == False:
            print("Match not found:", i)
    print('', file=bonout)

    print("[ dihedraltypes ]", file=bonout)
    print(";  ai    aj    ak     al   funct    c0           c1          c2          c3          c4          c5(kj/mol)",
          file=bonout)
    for i in newDihedrals:
        found = False
        for j in allDihedralParams:
            if ([i[0], i[1], i[2], i[3]] == [j[0], j[1], j[2], j[3]]) or (
                        [i[3], i[2], i[1], i[0]] == [j[0], j[1], j[2], j[3]]):
                if len(j) == 11:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:>12.6f}{9:>12.6f}{10:>12.6f}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), float(j[7]), float(j[8]),
                            float(j[9]), float(j[10])), file=bonout)
                elif len(j) == 8:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12d}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), int(j[7])), file=bonout)
                found = True
                break
            elif ((j[0] == 'X') and ([i[1], i[2], i[3]] == [j[1], j[2], j[3]])) or (
                        (j[0] == 'X') and ([i[2], i[1], i[0]] == [j[1], j[2], j[3]])):
                if len(j) == 11:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:>12.6f}{9:>12.6f}{10:>12.6f}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), float(j[7]), float(j[8]),
                            float(j[9]), float(j[10])), file=bonout)
                elif len(j) == 8:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12d}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), int(j[7])), file=bonout)
                found = True
                break
            elif ((j[3] == 'X') and ([i[0], i[1], i[2]] == [j[0], j[1], j[2]])) or (
                        (j[3] == 'X') and ([i[3], i[2], i[1]] == [j[0], j[1], j[2]])):
                if len(j) == 11:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:>12.6f}{9:>12.6f}{10:>12.6f}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), float(j[7]), float(j[8]),
                            float(j[9]), float(j[10])), file=bonout)
                elif len(j) == 8:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12d}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), int(j[7])), file=bonout)
                found = True
                break
            elif ((j[0] == 'X') and (j[3] == 'X') and ([i[1], i[2]] == [j[1], j[2]])) or (
                            (j[0] == 'X') and (j[3] == 'X') and ([i[2], i[1]] == [j[1], j[2]])):
                if len(j) == 11:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12.6f}{8:>12.6f}{9:>12.6f}{10:>12.6f}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), float(j[7]), float(j[8]),
                            float(j[9]), float(j[10])), file=bonout)
                elif len(j) == 8:
                    print(
                        "{0:>6s}{1:>6s}{2:>6s}{3:>6s}{4:>6d}{5:>12.6f}{6:>12.6f}{7:>12d}".format(
                            i[0], i[1], i[2], i[3], int(j[4]), float(j[5]), float(j[6]), int(j[7])), file=bonout)
                found = True
                break
            else:
                found = False
        if found == True:
            pass
        elif found == False:
            print("Match not found:", i, j)
    print('', file=bonout)

    print("[ dihedraltypes ]", file=bonout)
    print("; Improper dihedrals", file=bonout)
    for i in improperParams:
        print("{0:7s}{1:>20s}{2:>16.1f}{3:>14.3f}{4:>6d}".format(i[0], i[1], float(i[2]), float(i[3]), int(i[4])),
              file=bonout)
    print('', file=bonout)


if __name__ == '__main__':
    """
    Reads in a Gromacs .gro file for a molecule and outputs a molecular .itp file.
    """

    # Parse Command-line Input
    parser = argparse.ArgumentParser(description='Reads in a Gromacs .gro file for a molecule and outputs a molecular '
                                                 '.itp file. Options exist to create posre.itp and topol.top files.')
    parser.add_argument('-gro', nargs=1, help='Gromacs .gro Input File.', required=True)
    parser.add_argument('-out', nargs=1, help='Gromacs .itp Output File.', required=True)
    parser.add_argument('-gout', nargs=1, help='Gromacs .gro Output file.', default='NULL')
    parser.add_argument('-nb', nargs=1, help='Gromacs Non-Bonded .itp Input File.', required=True)
    parser.add_argument('-bon', nargs=1, help='Gromacs Bonded .itp Input File.', default=['NULL'])
    parser.add_argument('-ff', action='store_true',
                        help='Enable rewriting new force field files that remove wildcards.')
    parser.add_argument('-imp', action='store_true', help='Turn on finding impropers defined in -bon file. DO NOT USE!')
    parser.add_argument('-posre', action='store_true', help='Enable the Creation of a posre.itp File.')
    parser.add_argument('-gchg', nargs=2,
                        help='Turns on the extraction of charges from G09 files. Atom order should be identical to GRO file. G09_output_file charge_types[Mulliken, Hirschfeld, CM5, TXT]',
                        default=['False', 'Null'])
    parser.add_argument('--debug', action='store_true', help='Enable Debug Information. Does not print out files.')
    args = parser.parse_args()

    # Get charges from G09 File
    if vars(args)['gchg'][0] != 'False':
        if vars(args)['gchg'][1] in ['Mulliken', 'Hirschfeld', 'CM5', 'TXT']:
            gaussCharges = getG09Charges(vars(args)['gchg'][0], vars(args)['gchg'][1])
        else:
            print("The charge type input is not valid.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
    else:
        gaussCharges = []

    # Sanity Check
    if (vars(args)['bon'][0] == 'NULL') and (vars(args)['ff'] == True):
        print("You must define a bonded forcefield to print new formatted forcefield.", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    # Get Parameter Data, Bond, Angles, and Dihedrals
    resNum, name, atomType, atomNum, x, y, z = parseGro(vars(args)['gro'][0])
    nb = parseFF(vars(args)['nb'][0])
    if vars(args)['bon'][0] != 'NULL':
        bondParams, angleParams, dihedralParams, properDihedralParams, improperParams = parseBon(vars(args)['bon'][0])
    else:
        bondParams, angleParams, dihedralParams, properDihedralParams, improperParams = [], [], [], [], []
    radii, elements = findElements(atomType, nb)
    atomConnects = AtomConnections(len(atomType), x, y, z, radii)
    bonds = findBonds(atomConnects)
    angles = findAngles(atomConnects, bonds)
    dihedrals = findDihedrals(atomConnects, angles)

    # Determine If Debug Information Should Be Printed
    if vars(args)['debug'] == True:
        debugOut = open('debugITP.out', 'w')
        print('resNum', resNum, '', 'name', name, '', 'atomType', atomType, '', 'atomNum', atomNum, '', 'x', x, '', 'y',
              y, '', 'z', z, '', 'elements', elements, '', 'nb', nb, '', 'bondParams', bondParams, '', 'angleParams',
              angleParams, '', 'dihedralParams', dihedralParams, '', 'properDihedralParams', properDihedralParams, '',
              'improperParams', improperParams, '', 'radii', radii, '', 'bonds', bonds, '', 'angles', angles, '',
              'dihedrals', dihedrals, sep='\n', file=debugOut)
        sys.exit(0)

    # Print Output Files
    if vars(args)['bon'][0] != 'NULL':
        printTopol(vars(args)['out'][0], resNum, name, atomType, atomNum, nb, bonds,
                    angles, dihedrals, bondParams, angleParams, properDihedralParams, dihedralParams,
                   vars(args)['imp'], gaussCharges)
    else:
        printTopol_noparams(vars(args)['out'][0], resNum, name, atomType, atomNum, nb, bonds,
                   angles, dihedrals, gaussCharges)

    # Print Warning About Using the Improper Dihedral Algorithm
    if vars(args)['imp'] == True:
        print('\n')
        print("WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!")
        print('')
        print(
            "Improper Dihedrals are an experimental feature and should not be used as they fail to construct dihedrals")
        print("for rings correctly. YOU HAVE BEEN WARNED!")
        print("Make certain to check that Improper Dihedrals are correctly represented in the .itp file!")
        print('')
        print("WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!  WARNING!")
    else:
        pass

    # If 'gout' is given print a new gromacs file
    if vars(args)['gout'][0] != 'NULL':
        printGro(vars(args)['gout'][0], resNum, name, atomType, atomNum, x, y, z, nb)
    else:
        pass

    # Print position restraints if desired
    if vars(args)['posre'] == True:
        printPosRes(atomNum, elements)
    else:
        pass

    # Print new force field files if desired
    if vars(args)['ff'] == True:
        printFF(nb, bondParams, angleParams, properDihedralParams, dihedralParams, improperParams, atomType, bonds, angles, dihedrals)
    else:
        pass
