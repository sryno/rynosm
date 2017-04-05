#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void c_nm_to_ang(double *x, double *y, double *z, int natoms){
    double nmToAng = 10.0;
    int i;
#pragma omp parallel for ordered schedule(dynamic)
    for (i=0; i<natoms; i++){
        x[i] = x[i] * nmToAng;
        y[i] = y[i] * nmToAng;
        z[i] = z[i] * nmToAng;
    }
    return;
}


void c_ang_to_nm(double *x, double *y, double *z, int natoms){
    double nmToAng = 0.1;
    int i;
#pragma omp parallel for ordered schedule(dynamic)
    for (i=0; i<natoms; i++){
        x[i] = x[i] * nmToAng;
        y[i] = y[i] * nmToAng;
        z[i] = z[i] * nmToAng;
    }
    return;
}


void c_new_connections(int natoms, double *x, double *y, double *z, double *atom_radii, int max_atoms_molecule,
                int *connected, int iwidth){
    double temp_distance;
    int i, j;
    int bondNum = 0;
    for (i=0; i<natoms; i++){
        bondNum = 0;
        #pragma omp parallel for ordered
        for (j=0; j<natoms; j++){
            if (abs(i - j) < max_atoms_molecule){
                temp_distance = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]));
                if (temp_distance <= (1.5 * (atom_radii[i] + atom_radii[j]))){
                    connected[i*iwidth + bondNum] = j+1;
                    bondNum++;
                }
            }
        }
    }
    return;
}


void c_connections(int natoms, double *x, double *y, double *z, double *atom_radii,
            int max_atoms_molecule, int *connected, int iwidth){
    double temp_distance;
    int i, j;
    //#pragma omp parallel for ordered schedule(dynamic)//#pragma omp simd
    for (i=0; i<natoms; i++){
        #pragma omp parallel for ordered
        for (j=0; j<i; j++){
            if (abs(i - j) < max_atoms_molecule){
                temp_distance = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]));
                if (temp_distance <= (1.5 * (atom_radii[i] + atom_radii[j]))){
                    connected[j*iwidth + i] = 1;
                    connected[i*iwidth + j] = 1;
                }
            }
        }
        connected[i*iwidth + i] = 1;
    }
    return;
}


void c_new_molecule_ID(int *connections, int natoms, int iwidth, int *moleculeID, int nmolecules){
    int i, j, k, tmpndx;
    for (i=0; i<natoms; i++){
        if (moleculeID[i] == 0){
            nmolecules += 1;
        }
        moleculeID[i] = nmolecules;
        for (j=0; j<iwidth; j++){
            if (connections[i*iwidth +j] != 0){
                moleculeID[connections[i*iwidth + j]-1] = nmolecules;
                for (k=0; k<iwidth; k++){
                    tmpndx = connections[i*iwidth+j]-1;
                    if (connections[tmpndx*iwidth + k] != 0){
                        moleculeID[connections[tmpndx*iwidth + k]-1] = nmolecules;
                    }
                }
            }
        }
    }
    return;
}


void c_molecule_ID(int *connected, int iwidth, int max_atoms_molecule,
        int *moleculeID, int nmolecules){
    int i, j, k;
    int natoms = iwidth;
    for (i=0; i<natoms; i++){
        if (moleculeID[i] == 0){
            nmolecules += 1;
        }
        moleculeID[i] = nmolecules;
        for (j=(i-max_atoms_molecule); j<(i+1+max_atoms_molecule); j++){
            if ((j >= 0) && (j < natoms)){
                if ((abs(i-j) < max_atoms_molecule) && (connected[j*iwidth + i] == 1)){
                    moleculeID[j] = nmolecules;
                    for (k=(i-max_atoms_molecule); k<(i+1+max_atoms_molecule); k++){
                        if ((k >= 0) && (k < natoms)){
                            if (connected[k*iwidth + j] == 1){
                                moleculeID[k] = nmolecules;
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}
