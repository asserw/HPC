//
// Created by asser on 5/10/2024.
//

#ifndef VERLET_H
#define VERLET_H
#include "atoms.h" // Include the Atoms structure definition

// Function prototypes for Verlet integration steps using the Atoms structure
void verlet_step1(Atoms &atoms, double timestep);

void verlet_step2(Atoms &atoms, double timestep);

#endif  // VERLET_H
