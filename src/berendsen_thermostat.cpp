//
// Created by asser on 7/6/2024.
//
#include "../header/berendsen_thermostat.h"
#include "../header/compute_kinetic_energy.h"

const double k_B = 8.6173303e-5; // Boltzmann constant in eV/K

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time) {
    double kinetic_energy = compute_kinetic_energies(atoms);
    double current_temperature = 2.0 * kinetic_energy / (3.0 * atoms.nb_atoms() * k_B);
    double lambda = sqrt(1 + (timestep / relaxation_time) * (target_temperature / current_temperature - 1));
    atoms.velocities *= lambda;
}