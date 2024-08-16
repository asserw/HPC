//
// Created by asser on 6/19/2024.
//

#ifndef HPC_PROJECT_LJ_DIRECT_SUMMATION_NEIGHBORS_H
#define HPC_PROJECT_LJ_DIRECT_SUMMATION_NEIGHBORS_H
#include "atoms.h"
#include "neighbors.h"

double lj_direct_summation_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0, double cutoff = 2.5);

#endif // HPC_PROJECT_LJ_DIRECT_SUMMATION_NEIGHBORS_H
