//
// Created by asser on 6/18/2024.
//
#include "../header/atoms.h"
#include "../header/compute_kinetic_energy.h"
#include "../header/lj_direct_summation_neighbors.h"
#include "../header/neighbors.h"
#include "../header/verlet.h"
#include "../header/xyz.h"
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
    auto [names,positions, velocities] = read_xyz_with_velocities("lj54.xyz");

    Atoms atoms(positions, velocities);
    atoms.masses.setOnes();

    double timestep = 0.001;
    double total_time = 100.0;
    int num_steps = static_cast<int>(total_time / timestep);

    NeighborList neighbor_list;
    double cutoff = 2.5;

    std::ofstream traj("traj.xyz");

    for (int step = 0; step < num_steps; ++step) {
        verlet_step1(atoms, timestep);
        double potential_energy = lj_direct_summation_neighbors(atoms, neighbor_list, 1.0, 1.0, cutoff);
        verlet_step2(atoms, timestep);

        // Compute kinetic energy
        double kinetic_energy = compute_kinetic_energies(atoms);

        double total_energy = kinetic_energy + potential_energy;
        std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;

        // Write trajectory at every step
        write_xyz(traj, atoms);
    }

    traj.close();

    return 0;
}
