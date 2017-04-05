from __future__ import print_function

import argparse
import sys

__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = 'Sean M. Ryno'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def parsePDB(pdb, residueName):
    atomNum = []
    element = []
    atomType = []
    resName = []
    x = []
    y = []
    z = []

    fin = open(pdb, 'r')

    for line in fin:
        if (line[:6].strip() == 'HETATM') or (line[:6].strip() == 'ATOM'):
            atomNum.append(int(line[6:11].strip()))
            atomType.append(line[12:16].strip())
            if residueName == 'Null':
                resName.append(line[17:20].strip())
            else:
                resName.append(residueName)
            x.append(float(line[30:38].strip()))
            y.append(float(line[38:46].strip()))
            z.append(float(line[46:54].strip()))
            element.append(line[76:78].strip())

    return atomNum, element, atomType, resName, x, y, z


def parseXYZ(xyz, residueName):
    atomNum = []
    element = []
    atomType = []
    resName = []
    x = []
    y = []
    z = []
    finLines = []

    #residueName = input('That is the res name: ')

    fin = open(xyz, 'r')
    numLines = int(fin.readline().strip())
    fin.readline()
    for line in fin:
        finLines.append(line)
    for line in range(len(finLines)):
        if finLines[line].strip() == '':
            break
        else:
            cols = finLines[line].split()
            atomNum.append(line+1)
            element.append(cols[0])
            x.append(float(cols[1]))
            y.append(float(cols[2]))
            z.append(float(cols[3]))
            resName.append(residueName)

    return atomNum, element, atomType, resName, x, y, z


def parseXYZScan(xyz, residueName):
    atomNum = []
    element = []
    atomType = []
    resName = []
    x = []
    y = []
    z = []
    finLines = []

    fin = open(xyz, 'r')
    for line in fin:
        finLines.append(line)

    numLines = int(finLines[0].strip())
    molcNum = []

    for i in range(len(finLines)):
        if len(finLines[i].split()) > 0:
            if finLines[i].split()[0] == str(numLines):
                molcNum.append(i)
            else: pass
        else: pass

    for i in range(len(molcNum)):
        x.append([])
        y.append([])
        z.append([])
        if i == 0:
            for j in range(molcNum[i] + 2, molcNum[i] + 2 + numLines):
                cols = finLines[j].split()
                atomNum.append(j - 1)
                element.append(cols[0])
                x[i].append(float(cols[1]))
                y[i].append(float(cols[2]))
                z[i].append(float(cols[3]))
                resName.append(residueName)
        else:
            for j in range(molcNum[i] + 2, molcNum[i] + 2 + numLines):
                cols = finLines[j].split()
                x[i].append(float(cols[1]))
                y[i].append(float(cols[2]))
                z[i].append(float(cols[3]))

    return atomNum, element, atomType, resName, x, y, z


def readAtomTypes(txtfile):
    fin = open(txtfile, 'r')
    atomTypes = []
    for line in fin:
        atomTypes.append(line.strip())

    fin.close()

    return atomTypes


def changeAtomTypes(oldAtomTypes):
    """
    Converts old Atom Types to new Atom Types. Uses Dictionary data types to prevent asking the same question twice.
    :param oldAtomTypes:
    :return:
    """
    typesDict = {}
    newAtomTypes = []
    for i in oldAtomTypes:
        if i not in typesDict:
            aType = input('old Type: ${0:s}   new Name: '.format(i))
            typesDict[i] = aType
        else:
            pass
    for i in oldAtomTypes:
        newAtomTypes.append(typesDict[i])

    return newAtomTypes


def angToNm(x):
    """
    Converts Angstroms to Nanometers.
    :param x:
    :return:
    """
    for i in range(len(x)):
        x[i] = x[i] * 0.1

    return x


def writeGro(groFile, numAtoms, atomNums, elements, x, y, z, atomTypes, sysName):
    """
    Writes a .gro formatted coordinate file. Note that .gro has limited accuracy.
    :param groFile:
    :param numAtoms:
    :param atomNums:
    :param elements:
    :param x:
    :param y:
    :param z:
    :param atomTypes:
    :param sysName:
    :return:
    """
    fout = open(groFile, 'w')
    resID = 1
    print(sysName[0], file=fout)
    print('', numAtoms, file=fout)
    for i in range(len(atomNums)):
        print('{0:>5d}{1:<5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}    ;{7:>5s}'.format(resID, sysName[i], atomTypes[i],
                                                                                       atomNums[i],
                                                                                       x[i], y[i], z[i], elements[i]),
              file=fout)


if __name__ == '__main__':
    """
    Allows for the conversion of a PDB, XYZ, or G09 formatted file to Gromacs GRO format.
    """

    parser = argparse.ArgumentParser(description='Allows for the conversion of a PDB, XYZ, or G09 formatted file to '
                                                 'Gromacs GRO format.')
    parser.add_argument('INPUT_FILE', nargs=2, help='Input File Name and file type [PDB, XYZ, G09]')
    parser.add_argument('GROMACS_FILE', nargs=1, help='Gromacs Gro Output File Prefix.')
    parser.add_argument('-scan', action='store_true', help='Turns on Parsing XYZ files with multiple structures from scans.')
    parser.add_argument('--new', action='store_true', help='Switch to turn on asking for new atom types.')
    parser.add_argument('-NT', nargs=1, help='Extracts atom types from an external file. Required for XYZ and G09 files. [Filename]',
                        default=['Null'])
    parser.add_argument('-name', nargs=1, help='System Name. Only used for XYZ and PDB.', default=['Null'])
    args = parser.parse_args()

    if vars(args)['scan'] == False:
        if vars(args)['INPUT_FILE'][1] == 'PDB':
            atomNum, element, atomType, resName, x, y, z = parsePDB(vars(args)['INPUT_FILE'][0], vars(args)['name'][0])
        elif vars(args)['INPUT_FILE'][1] == 'XYZ':
            atomNum, element, atomType, resName, x, y, z = parseXYZ(vars(args)['INPUT_FILE'][0], vars(args)['name'][0])
            atomType = readAtomTypes(vars(args)['NT'][0])
        elif vars(args)['INPUT_FILE'][1] == 'G09':
            print("This does not yet function.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
            atomNum, element, atomType, resName, x, y, z = parseG09(vars(args)['INPUT_FILE'][0], vars(args)['name'][0])
            atomType = readAtomTypes(vars(args)['NT'][0])
        x = angToNm(x)
        y = angToNm(y)
        z = angToNm(z)
        if vars(args)['new'] == True:
            atomType = changeAtomTypes(atomType)
        elif vars(args)['NT'][0] != 'Null':
            atomType = readAtomTypes(vars(args)['NT'][0])
        else:
            pass

        writeGro(vars(args)['GROMACS_FILE'][0] + '.gro', len(atomNum), atomNum, element, x, y, z, atomType, resName)

    elif vars(args)['scan'] == True:
        atomNum, element, atomType, resName, X, Y, Z = parseXYZScan(vars(args)['INPUT_FILE'][0], vars(args)['name'][0])
        for i in range(len(X)):
            X[i] = angToNm(X[i])
            Y[i] = angToNm(Y[i])
            Z[i] = angToNm(Z[i])
        atomType = readAtomTypes(vars(args)['NT'][0])
        for i in range(len(X)):
            j = i * 10
            writeGro(vars(args)['GROMACS_FILE'][0] + '_' + str(j) + '.gro', len(atomNum), atomNum, element, X[i], Y[i], Z[i], atomType, resName)