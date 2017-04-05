from __future__ import print_function

import sys, os
import numpy as np
from scipy import *
from scipy.sparse import lil_matrix
import argparse
import conn_molID


__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = 'Sean M. Ryno'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def assignElements(aNames):
    """
    Assigns elements based on common GROMACS naming nomenclature.
    """

    elmnts = np.ndarray.tolist(np.zeros(len(aNames)).astype(str))

    for i in range(len(aNames)):
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
    Determines the mass of atoms based on Element type. Some values have been modified from their reported covalent
    radii to not give unrealistic bond distances. E.g., C-Cl would allow bonds of 2.5 Angstroms. Reference values
    are commented out.
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


def cy_connections(natoms, x1, y1, z1, atom_radii, max_bonds, max_atoms_molecule):
    """
    Creates an [natoms, max_bonds+1] array that stores all connected atoms for a given index. All numbering
    starts at 0.
    """

    x1 = np.asarray(x1, dtype=np.float64)
    y1 = np.asarray(y1, dtype=np.float64)
    z1 = np.asarray(z1, dtype=np.float64)
    atom_radii = np.asarray(atom_radii, dtype=np.float64)
    connected = np.asarray(np.zeros((natoms, max_bonds+1), dtype=np.int32))

    conn_molID.cyc_connections(natoms, x1, y1, z1, atom_radii, max_atoms_molecule, connected, max_bonds+1)

    return connected


def connections_old(natoms, x, y, z, atom_radii, max_atoms_molecule):
    """
        Creates sparse matrix that stores atom connections
        Atoms are defined to be connected to themselves
    """

    connected = lil_matrix((natoms, natoms))

    #The below variable uses a non-sparse matrix that is faster, but can use huge amounts of memory.
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


def cy_molecule_ID(connections, natoms):
    """
    Creates a 1-D array of length [natoms] that contains the molecule number to which each atom belongs. Number of
    molecules starts at 1 while atom numbering (array index) starts at 0.
    """

    nmolecules = int(0)
    iwidth = connections.shape[1]
    moleculeID = np.asarray(np.zeros(natoms, dtype=np.int32), dtype=np.int32)

    conn_molID.cyc_molecule_ID(connections, natoms, iwidth, moleculeID, nmolecules)

    return moleculeID, nmolecules


def molecule_ID_cyc_connections(connections, natoms):
    """
    Same as cy_molecule_ID, but written in Python for easier trouble shooting.
    """

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


def old_molecule_ID(natoms, connections, max_atoms_molecule):
    """
        Returns the modecule ID of each atom
        Very time consuming and should not be used as not always accurate.
        Only works with connections_old
    """

    nmolecules = 0
    moleculeID = zeros(natoms)
    for i in range(natoms):
        if moleculeID[i] == 0:
            nmolecules += 1
            molecule_string = zeros(natoms)
            molecule_string[i] = 1
        for j in range(natoms):
            if (abs(i - j) < max_atoms_molecule) and (connections[i, j] == 1):
                for k in range(natoms):
                    if (connections[j, k] == 1):
                        molecule_string[k] = 1
                molecule_string[j] = 1
        for j in range(natoms):
            if (molecule_string[j] == 1):
                moleculeID[j] = nmolecules

    return moleculeID, nmolecules


if __name__ == '__main__':
    """
    Provides an importable module for python to determine atomic bonds (connections) and which series of bonds form
    molecules. This should not be used as a standalone program, but rather as a framework for other programs to link to.
    """

    parser = argparse.ArgumentParser(description='Provides an importable module for python to determine atomic bonds '
                                                 '(connections) and which series of bonds form molecules. This should '
                                                 'not be used as a standalone program, but rather as a framework for '
                                                 'other programs to link to.')
    parser.add_argument('NULL', nargs='*', help='NULL')
    args = parser.parse_args()

    # The below provide examples of how to assign elements, masses, and radii. Each takes as input a 1-D array.
    # elements = assignElements(atomNames)
    # atomMasses = assignMass(elements)
    # atomRadii = assignRadii(elements)

    # The below provide examples of how to call the connection routines.
    # connections = connections_old(len(x), x, y, z, atomRadii, max_atoms_molecule)
    # connections = cy_connections(len(x), x, y, z, atomRadii, maxBonds, max_atoms_molecule)

    # The blow provide examples of how to call the MoleculeID routines.
    # moleculeIDs, numMolecules = cy_molecule_ID(connections, len(x))
    # moleculeIDs, numMolecules = molecule_ID_cyc_connections(connections, len(x))
    # moleculeIDs, numMolecules = old_molecule_ID(len(x), connections, max_atoms_molecule)

    print("This program is not intended to be used as a standalone. Please import the provided functions into your"
          "programs as necessary.", file=sys.stderr)
    print("Exiting...", file=sys.stderr)
    sys.exit(1)

