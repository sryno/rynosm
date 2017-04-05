import cython
cimport cython
import numpy as np
cimport numpy as np


cdef extern from "c_conn_molID.c":
    void c_connections(int natoms, double *x, double *y, double *z, double *atom_radii,
            int max_atoms_molecule, int *connected, int iwidth)


@cython.boundscheck(False)
def cyc_connections(int natoms,
                    np.ndarray[np.float64_t, ndim=1] x not None,
                    np.ndarray[np.float64_t, ndim=1] y not None,
                    np.ndarray[np.float64_t, ndim=1] z not None,
                    np.ndarray[np.float64_t, ndim=1] atom_radii not None,
                    int max_atoms_molecule,
                    np.ndarray[np.int32_t, ndim=2] connected not None,
                    int iwidth):
    c_connections(<int> natoms, <double*> x.data, <double*> y.data, <double*> z.data, <double *> atom_radii.data,
                  <int> max_atoms_molecule, <int*> connected.data, <int> iwidth)


cdef extern from "c_conn_molID.c":
    void c_molecule_ID(int *connections, int natoms, int iwidth, int *moleculeID, int nmolecules)


def cyc_molecule_ID(
        np.ndarray[np.int32_t, ndim=2] connections not None,
        int natoms, int iwidth,
        np.ndarray[np.int32_t, ndim=1] moleculeID not None,
        int nmolecules):
    c_molecule_ID(<int*> connections.data, <int> natoms, <int> iwidth, <int*> moleculeID.data, <int> nmolecules)
