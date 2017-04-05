#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


void c_connections(int natoms, double *x, double *y, double *z, double *atom_radii, int max_atoms_molecule,
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
                    connected[i*iwidth + bondNum] = j;
                    bondNum++;
                }
            }
        }
    }
    return;
}


void c_molecule_ID(int *connections, int natoms, int iwidth, int *moleculeID, int nmolecules){
    int i, j, k, tmpndx;
    for (i=0; i<natoms; i++){
        if (moleculeID[i] == 0){
            nmolecules += 1;
        }
        moleculeID[i] = nmolecules;
        for (j=0; j<iwidth; j++){
            if (connections[i*iwidth +j] != 0){
                moleculeID[connections[i*iwidth + j]] = nmolecules;
                for (k=0; k<iwidth; k++){
                    tmpndx = connections[i*iwidth+j];
                    if (connections[tmpndx*iwidth + k] != 0){
                        moleculeID[connections[tmpndx*iwidth + k]] = nmolecules;
                    }
                }
            }
        }
    }
    return;
}
