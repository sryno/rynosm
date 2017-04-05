from __future__ import print_function

import sys, os
import numpy as np
from scipy import *
from scipy.sparse import lil_matrix
import argparse
import nm_conn_molID


__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = 'Sean M. Ryno'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def parseGro(inFile):
    fin = open(inFile, 'r')

    mName, aName, x, y, z, cAxes = [], [], [], [], [], 0

    fin.readline()
    nAtoms = int(fin.readline().strip())
    for line in fin:
        try:
            float(line.split()[0])
            cAxes = [line.split()[0], line.split()[1], line.split()[2], line.split()[3], line.split()[4],
                     line.split()[5], line.split()[6], line.split()[7], line.split()[8]]
            break
        except (ValueError, IndexError):
            try:
                temp = line.split()[1]
                mName.append(line[5:10].strip())
                aName.append(line[10:15].strip())
                x.append(float(line[20:28].strip()))
                y.append(float(line[28:36].strip()))
                z.append(float(line[36:44].strip()))
            except (ValueError, IndexError):
                pass

    return nAtoms, mName, aName, x, y, z, cAxes


def assignElements(aNames):
    elmnts = np.ndarray.tolist(np.zeros(len(aNames)).astype(str))

    for i in range(len(atomNames)):
        if aNames[i][:2].upper() == 'CL':
            elmnts[i] = 'Cl'
        elif aNames[i][0].upper() == 'H':
            elmnts[i] = 'H'
        elif aNames[i][0].upper() == 'C':
            elmnts[i] = 'C'
        elif aNames[i][0].upper() == 'S':
            elmnts[i] = 'S'
        elif aNames[i][0].upper() == 'N':
            elmnts[i] = 'N'
        elif aNames[i][0].upper() == 'O':
            elmnts[i] = 'O'
        elif aNames[i][0].upper() == 'F':
            elmnts[i] = 'F'
        else:
            print("The following atom name was not recognized:", aNames[i], file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)

    return elmnts


def assignMass(elmnts):
    """
    Determines the mass of atoms based on Element type
    """
    aMass = np.ndarray.tolist(np.zeros(len(elmnts)).astype(str))
    for i in range(len(elmnts)):
        if elmnts[i] == 'X':
            aMass[i] = 0.0
        elif elmnts[i] == 'H':
            aMass[i] = 1.00794
        elif elmnts[i] == 'C':
            aMass[i] = 12.011
        elif elmnts[i] == 'F':
            aMass[i] = 18.998403
        elif elmnts[i] == 'Si':
            aMass[i] = 28.0855
        elif elmnts[i] == 'S':
            aMass[i] = 32.066
        elif elmnts[i] == 'N':
            aMass[i] = 14.0067
        elif elmnts[i] == 'O':
            aMass[i] = 15.9994
        elif elmnts[i] == 'Cl':
            aMass[i] = 35.4527
        else:
            print("Error in assigning masses:", i, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)

    return aMass


def assignRadii(elmnts):
    """
    Determines the mass of atoms based on Element type. Some values have been modofied from their reported covalent
    radii to not give unrealistic bond distances. E.g., C-Cl would allow bonds of 2.5 Angstroms.
    """
    aRadii = np.ndarray.tolist(np.zeros(len(elmnts)).astype(str))
    for i in range(len(elmnts)):
        if elmnts[i] == 'X':
            aRadii[i] = 0.23
        elif elmnts[i] == 'H':
            aRadii[i] = 0.23
        elif elmnts[i] == 'C':
            aRadii[i] = 0.68
        elif elmnts[i] == 'F':
            aRadii[i] = 0.64
        elif elmnts[i] == 'Si':
            #aRadii[i] = 1.2
            aRadii[i] = 0.80
        elif elmnts[i] == 'S':
            #aRadii[i] = 1.02
            aRadii[i] = 0.68
        elif elmnts[i] == 'N':
            #aRadii[i] = 0.68
            aRadii[i] = 0.65
        elif elmnts[i] == 'O':
            #aRadii[i] = 0.68
            aRadii[i] = 0.60
        elif elmnts[i] == 'Cl':
            #aRadii[i] = 0.99
            aRadii[i] = 0.68
        else:
            print("Error in assigning radius:", i, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)

    return aRadii


def assignDonorAccept(mNames):
    donAccpt = np.ndarray.tolist(np.zeros(len(mNames)).astype(str))
    uniqueNames = []

    for i in mNames:
        if i not in uniqueNames:
            uniqueNames.append(i)

    userNames = np.ndarray.tolist(np.zeros(len(uniqueNames)).astype(str))
    userNames = {}
    for i in range(len(uniqueNames)):
        userNames[uniqueNames[i]] = input("Is " + uniqueNames[i] + " donor [D] or acceptor [A]?").upper()

    return userNames


def cy_nmToAng(x1, y1, z1):
    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    z1 = np.asarray(z1, dtype=np.float64)

    nm_conn_molID.cyc_nm_to_ang(x1, y1, z1, len(x1))

    x1 = x1.tolist()
    y1 = y1.tolist()
    z1 = z1.tolist()

    return x1, y1, z1


def cy_angToNm(x1, y1, z1):
    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    z1 = np.asarray(z1, dtype=np.float64)

    nm_conn_molID.cyc_ang_to_nm(x1, y1, z1, len(x1))

    x1 = x1.tolist()
    y1 = y1.tolist()
    z1 = z1.tolist()

    return x1, y1, z1


def cy_connections(natoms, x1, y1, z1, atom_radii):

    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    z1 = np.asarray(z1, dtype=np.float64)
    atom_radii = np.asarray(atom_radii, dtype=np.float64)

    q = np.zeros((natoms, natoms), dtype=np.int32)

    connected = np.asarray(q, dtype=np.int32)

    nm_conn_molID.cyc_connections(natoms, x1, y1, z1, atom_radii, max_atoms_molecule, connected, natoms)

    return connected


def cy_new_connections(natoms, x1, y1, z1, atom_radii):

    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    z1 = np.asarray(z1, dtype=np.float64)
    atom_radii = np.asarray(atom_radii, dtype=np.float64)
    connected = np.asarray(np.zeros((natoms, max_bonds+1), dtype=np.int32))

    nm_conn_molID.cyc_new_connections(natoms, x1, y1, z1, atom_radii, max_atoms_molecule, connected, max_bonds+1)

    return connected


def connections_old(natoms, x, y, z, atom_radii):
    """
        Creates sparse matrix that stores atom connections
        Atoms are defined to be connected to themselves
    """

    connected = lil_matrix((natoms, natoms))
    #connected = np.zeros((natoms, natoms), dtype=np.int32)

    for i in range(natoms):
        for j in range(i):
            temp_distance = 0.0
            if (abs(i - j) < max_atoms_molecule):
                temp_distance = sqrt(
                    (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]))
                if temp_distance <= (1.5 * (atom_radii[i] + atom_radii[j])):
                    connected[i, j] = 1
                    connected[j, i] = 1
            else:
                pass
        connected[i, i] = 1

    return connected


def cy_molecule_ID(connections):

    nmolecules = int(0)
    natoms = connections.shape[0]
    q = np.zeros(natoms, dtype=np.int32)

    moleculeID = np.asarray(q, dtype=np.int32)

    nm_conn_molID.cyc_molecule_ID(connections, natoms, max_atoms_molecule, moleculeID, nmolecules)

    return moleculeID, nmolecules


def cy_new_molecule_ID(connections, natoms):

    nmolecules = int(0)
    iwidth = connections.shape[1]
    moleculeID = np.asarray(np.zeros(natoms, dtype=np.int32), dtype=np.int32)

    nm_conn_molID.cyc_new_molecule_ID(connections, natoms, iwidth, moleculeID, nmolecules)

    return moleculeID, nmolecules


def new_molecule_ID(connections, natoms):

    nmolecules = int(0)
    iwidth = connections.shape[1]
    moleculeID = np.asarray(np.zeros(natoms, dtype=np.int32), dtype=np.int32)

    for i in range(natoms):
        if (moleculeID[i] == 0):
            nmolecules += 1
        moleculeID[i] = nmolecules
        for j in range(iwidth):
            if connections[i][j] != 0:
                moleculeID[connections[i][j] - 1] = nmolecules
                for k in range(iwidth):
                    tmpndx = connections[i][j] - 1
                    if connections[connections[i][j] - 1][k] != 0:
                        print(tmpndx)
                        moleculeID[connections[tmpndx][k] - 1] = nmolecules

    return moleculeID, nmolecules


def writeConn(connct, fileOut):
    fout = open(fileOut + '.conn', 'w')
    d = lambda n, minn, maxn: max(min(maxn, n), minn)
    f = lambda a: (abs(a) + a) / 2

    nAtm = connct.shape[0]
    print("Atom connections. Number of Atoms:", nAtm, file=fout)

    for i in range(nAtm):
        print(i+1, sep=' ', end=' ', file=fout)
        for j in range(int(f(i - max_atoms_molecule)), int(d((i + max_atoms_molecule), 0, nAtm))):
            if abs(i-j) <= max_atoms_molecule:
                if (connct[i, j] != 0) and (i != j):
                    print(j+1, sep=' ', end=' ', file=fout)
        print('', file=fout)

    fout.close()

    return


def writeConn_new(nAtm, connct, fileOut):
    fout = open(fileOut + '.conn', 'w')
    d = lambda n, minn, maxn: max(min(maxn, n), minn)
    f = lambda a: (abs(a) + a) / 2

    print("Atom connections. Number of Atoms:", nAtm, file=fout)

    for i in range(nAtm):
        print(i + 1, sep=' ', end=' ', file=fout)
        for j in range(len(connct[i])):
            if (connct[i][j] != 0) and (connct[i][j] != i+1):
                print(connct[i][j], sep=' ', end=' ', file=fout)
        print('', file=fout)

    fout.close()

    return


def writeMolcID(molID, fileOut):
    fout = open(fileOut + '.mol', 'w')

    nAtm = molID.shape[0]
    print("Molecule IDs by atom. Number of Atoms:", nAtm, file=fout)

    for i in range(nAtm):
        print(i+1, molID[i], file=fout)

    fout.close()

    return


def writeStruct(molID, x1, y1, z1, dA, mNames, fileOut):
    fout = open(fileOut + '.struct', 'w')

    nAtm = molID.shape[0]
    print("Structure file for", fileOut, "     Number of Atoms:", nAtm, file=fout)

    for i in range(nAtm):
        print("{0:<5d} {1:>10.5f} {2:>10.5f} {3:>10.5f} {4:>2s}".format(molID[i], x1[i], y1[i], z1[i], dA[mNames[i]]),
              file=fout)

    fout.close()

    return


if __name__ == '__main__':
    """
    Reads in a Gromacs GRO file and outputs XYZ, Connections, and Molecule IDs.
    """

    parser = argparse.ArgumentParser(description='Reads in a Gromacs GRO file and outputs XYZ, Connections, and Molecule IDs.')
    parser.add_argument('INPUT_FILE', nargs=1, help='Input Gromacs File Name')
    parser.add_argument('OUTPUT_NAME', nargs=1, help='Name Prefix for Output Files')
    parser.add_argument('-maxatm', nargs=1, help="Sets the maximum molecule size. Default: 150", default=['150'])
    parser.add_argument('-maxbnd', nargs=1, help="Sets the maximum number of bonds an atom can have. Default: 4", default=['4'])
    parser.add_argument('-conn', action='store_true', help="Turn on outputting of atom connections.")
    parser.add_argument('-mol', action='store_true', help="Turn on outputting of molecule ID to separate file.")
    args = parser.parse_args()

    max_atoms_molecule = int(vars(args)['maxatm'][0])
    max_bonds = int(vars(args)['maxbnd'][0])
    numAtoms, molNames, atomNames, x, y, z, cryAxes = parseGro(vars(args)['INPUT_FILE'][0])
    elements = assignElements(atomNames)
    atomMasses = assignMass(elements)
    atomRadii = assignRadii(elements)
    donorAcceptor = assignDonorAccept(molNames)

    x, y, z = cy_nmToAng(x, y, z)

    print("Determining atom connections.", file=sys.stdout)
    #connections = cy_connections(len(x), x, y, z, atomRadii)
    #connections = []
    #connections = connections_old(len(x), x, y, z, atomRadii)

    connections = cy_new_connections(len(x), x, y, z, atomRadii)

    print("Connections Done.", file=sys.stdout)

    print("Determining molcule IDs.", file=sys.stdout)
    #moleculeIDs, numMolecules = cy_molecule_ID(connections)
    moleculeIDs, numMolecules = cy_new_molecule_ID(connections, len(x))
    #moleculeIDs, numMolecules = new_molecule_ID(connections, len(x))
    print("Molecule IDs Done.", file=sys.stdout)

    x, y, z = cy_angToNm(x, y, z)

    if vars(args)['conn'] == True:
        print("Writing atom connections.", file=sys.stdout)
        #writeConn(connections, vars(args)['OUTPUT_NAME'][0])
        writeConn_new(len(x), connections, vars(args)['OUTPUT_NAME'][0])

    if vars(args)['mol'] == True:
        print("Writing molecule IDs.", file=sys.stdout)
        writeMolcID(moleculeIDs, vars(args)['OUTPUT_NAME'][0])

    print("Writing output structure file.", file=sys.stdout)
    writeStruct(moleculeIDs, x, y, z, donorAcceptor, molNames, vars(args)['OUTPUT_NAME'][0])

