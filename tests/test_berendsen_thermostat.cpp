//
// Created by asser on 7/6/2024.
//
#include <gtest/gtest.h>
#include "../header/berendsen_thermostat.h"
#include "../header/compute_kinetic_energy.h"

TEST(BerendsenThermostatTest, TemperatureControl) {
    constexpr int nb_atoms = 1;
    constexpr double target_temperature = 300.0; // K
    constexpr double timestep = 0.001;
    constexpr double relaxation_time = 0.1;
    constexpr double total_time = 100.0;
    constexpr int num_steps = static_cast<int>(total_time / timestep);
    const double k_B = 8.6173303e-5; // Boltzmann constant in eV/K


    Atoms atoms(nb_atoms);
    atoms.velocities.setRandom();  // Initialize random velocities

    // Apply Berendsen thermostat for a sufficient number of steps
    for (int i = 0; i < num_steps; ++i) {
        berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);
    }

    // Compute final temperature
    double final_temperature = 2.0 * compute_kinetic_energies(atoms) / (3.0 * nb_atoms * k_B);

    // Check that the final temperature is close to the target temperature
    EXPECT_NEAR(final_temperature, target_temperature, 1.0);
}
