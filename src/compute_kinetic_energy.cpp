//
// Created by asser on 6/18/2024.
//
#include "../header/compute_kinetic_energy.h"

double compute_kinetic_energies(Atoms &atoms) {
    double kinetic_energy = 0.0;
    atoms.kinetic_energies.setZero();
    for (size_t i = 0; i < atoms.nb_atoms(); ++i) {
        atoms.kinetic_energies(i) = 0.5 * atoms.masses(i) * atoms.velocities.col(i).matrix().squaredNorm();
        kinetic_energy += atoms.kinetic_energies(i);
    }
    return kinetic_energy;
}

