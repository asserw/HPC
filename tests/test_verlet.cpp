#include "../header/atoms.h"
#include "../header/verlet.h"
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

void run_verlet(Atoms &atoms, int num_steps, double timestep) {
    for (int i = 0; i < num_steps; ++i) {
        verlet_step1(atoms, timestep);
        // Forces are constant, no need to recompute here if using a simple constant force
        verlet_step2(atoms, timestep);
    }
}

void calculate_expected_results(const Atoms &atoms, double timestep, int num_steps, Eigen::Array3Xd &expected_positions, Eigen::Array3Xd &expected_velocities) {
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        // Time
        double time = num_steps * timestep;

        // Expected positions and velocities considering constant force
        expected_positions(0, i) = atoms.positions(0, i) + atoms.velocities(0, i) * time + 0.5 * time * time; // 0.5 * a * t^2 where a = 1
        expected_positions(1, i) = atoms.positions(1, i) + atoms.velocities(1, i) * time + 0.5 * time * time;
        expected_positions(2, i) = atoms.positions(2, i) + atoms.velocities(2, i) * time + 0.5 * time * time;

        expected_velocities(0, i) = atoms.velocities(0, i) + time; // a * t where a = 1
        expected_velocities(1, i) = atoms.velocities(1, i) + time;
        expected_velocities(2, i) = atoms.velocities(2, i) + time;
    }
}

void check_results(const Atoms &atoms, const Eigen::Array3Xd &expected_positions, const Eigen::Array3Xd &expected_velocities) {
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        // Validate positions and velocities
        ASSERT_NEAR(atoms.positions(0, i), expected_positions(0, i), 1e-6);
        ASSERT_NEAR(atoms.positions(1, i), expected_positions(1, i), 1e-6);
        ASSERT_NEAR(atoms.positions(2, i), expected_positions(2, i), 1e-6);

        ASSERT_NEAR(atoms.velocities(0, i), expected_velocities(0, i), 1e-6);
        ASSERT_NEAR(atoms.velocities(1, i), expected_velocities(1, i), 1e-6);
        ASSERT_NEAR(atoms.velocities(2, i), expected_velocities(2, i), 1e-6);
    }
}

TEST(VerletTest, SingleAtom) {
    int nbAtoms = 1;
    Atoms atoms(nbAtoms);

    // Set initial conditions to random values
    atoms.positions.setRandom();
    atoms.velocities.setRandom();
    atoms.forces.setConstant(1.0); // Apply a force of 1.0 in all directions

    double timestep = 1.0;
    int num_steps = 100;

    // Calculate expected results before running the integration
    Eigen::Array3Xd expected_positions(3, nbAtoms);
    Eigen::Array3Xd expected_velocities(3, nbAtoms);
    calculate_expected_results(atoms, timestep, num_steps, expected_positions, expected_velocities);

    // Run the Velocity-Verlet integration for 10 timesteps
    run_verlet(atoms, num_steps, timestep);

    // Check the results
    check_results(atoms, expected_positions, expected_velocities);
}

TEST(VerletTest, MultipleAtoms) {
    int nbAtoms = 100;
    Atoms atoms(nbAtoms);

    // Set initial conditions to random values
    atoms.positions.setRandom();
    atoms.velocities.setRandom();
    atoms.forces.setConstant(1.0);  // Set all forces to 1.0

    double timestep = 1.0;
    int num_steps = 100;

    // Calculate expected results before running the Verlet integration
    Eigen::Array3Xd expected_positions(3, nbAtoms);
    Eigen::Array3Xd expected_velocities(3, nbAtoms);
    calculate_expected_results(atoms, timestep, num_steps, expected_positions, expected_velocities);

    // Run the Velocity-Verlet integration
    run_verlet(atoms, num_steps, timestep);

    // Check positions and velocities after Verlet integration
    check_results(atoms, expected_positions, expected_velocities);
}
