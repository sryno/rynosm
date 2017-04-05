import cython
cimport cython

import numpy as np
cimport numpy as np


cdef extern from "c_nm_conn_molID.c":
    void c_nm_to_ang(double *x, double *y, double *z, int natoms)


@cython.boundscheck(False)
def cyc_nm_to_ang(np.ndarray[np.float64_t, ndim=1] x not None, np.ndarray[np.float64_t, ndim=1] y not None,
                  np.ndarray[np.float64_t, ndim=1] z not None, int natoms):
    c_nm_to_ang(<double*> x.data, <double*> y.data, <double*> z.data, <int> natoms)


cdef extern from "c_nm_conn_molID.c":
    void c_ang_to_nm(double *x, double *y, double *z, int natoms)


@cython.boundscheck(False)
def cyc_ang_to_nm(np.ndarray[np.float64_t, ndim=1] x not None, np.ndarray[np.float64_t, ndim=1] y not None,
                  np.ndarray[np.float64_t, ndim=1] z not None, int natoms):
    c_ang_to_nm(<double*> x.data, <double*> y.data, <double*> z.data, <int> natoms)


cdef extern from "c_nm_conn_molID.c":
    void c_new_connections(int natoms, double *x, double *y, double *z, double *atom_radii,
            int max_atoms_molecule, int *connected, int iwidth)


@cython.boundscheck(False)
def cyc_new_connections(int natoms,
                    np.ndarray[np.float64_t, ndim=1] x not None,
                    np.ndarray[np.float64_t, ndim=1] y not None,
                    np.ndarray[np.float64_t, ndim=1] z not None,
                    np.ndarray[np.float64_t, ndim=1] atom_radii not None,
                    int max_atoms_molecule,
                    np.ndarray[np.int32_t, ndim=2] connected not None,
                    int iwidth):
    c_new_connections(<int> natoms, <double*> x.data, <double*> y.data, <double*> z.data, <double *> atom_radii.data,
                  <int> max_atoms_molecule, <int*> connected.data, <int> iwidth)



cdef extern from "c_nm_conn_molID.c":
    void c_connections(int natoms, double *x, double *y, double *z, double *atom_radii, int max_atoms_molecule,
                       int *connected, int iwidth)


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


cdef extern from "c_nm_conn_molID.c":
    void c_new_molecule_ID(int *connections, int natoms, int iwidth, int *moleculeID, int nmolecules)


def cyc_new_molecule_ID(
        np.ndarray[np.int32_t, ndim=2] connections not None,
        int natoms, int iwidth,
        np.ndarray[np.int32_t, ndim=1] moleculeID not None,
        int nmolecules):
    c_new_molecule_ID(<int*> connections.data, <int> natoms, <int> iwidth, <int*> moleculeID.data, <int> nmolecules)


cdef extern from "c_nm_conn_molID.c":
    void c_molecule_ID(int *connections, int iwidth, int max_atoms_molecule, int *moleculeID, int nmolecules)


def cyc_molecule_ID(
        np.ndarray[np.int32_t, ndim=2] connections not None,
        int iwidth,
        int max_atoms_molecule,
        np.ndarray[np.int32_t, ndim=1] moleculeID not None,
        int nmolecules):
    c_molecule_ID(<int*> connections.data, <int> iwidth, <int> max_atoms_molecule, <int*> moleculeID.data,
        <int> nmolecules)
