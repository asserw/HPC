//
// Created by asser on 5/10/2024.
//
#include "../header/verlet.h"

// First step of the Velocity-Verlet integration: Update positions and half-step velocities
void verlet_step1(Atoms &atoms, double timestep) {
    // Update velocities at the half timestep and positions using Eigen's vectorized operations
    atoms.velocities += 0.5 * atoms.forces * timestep;
    atoms.positions += atoms.velocities * timestep;
}

// Second step of the Velocity-Verlet integration: Complete the velocity update
void verlet_step2(Atoms &atoms, double timestep) {
    // Complete the update of velocities with the new forces using Eigen's vectorized operations
    atoms.velocities += 0.5 * atoms.forces * timestep;
}


