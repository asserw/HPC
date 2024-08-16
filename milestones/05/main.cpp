//
// Created by asser on 7/6/2024.
//
#include "../header/atoms.h"
#include "../header/berendsen_thermostat.h"
#include "../header/lj_direct_summation.h"
#include "../header/verlet.h"
#include "../header/xyz.h"
#include <chrono>
#include <iostream>

void equilibrate_system(Atoms &atoms, double target_temperature, double timestep, double relaxation_time, int num_steps) {
    for (int i = 0; i < num_steps; ++i) {
        verlet_step1(atoms, timestep);
        lj_direct_summation(atoms, 1.0, 1.0);  // Compute forces using normal LJ summation
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);  // Apply thermostat
    }
}

void measure_performance(int nx, int ny, int nz, double total_time, double timestep, double lattice_constant) {
    int num_steps = static_cast<int>(total_time / timestep);
    Atoms atoms(nx, ny, nz, lattice_constant);

    auto start = std::chrono::high_resolution_clock::now();
    equilibrate_system(atoms, 300.0, timestep, 0.1, num_steps);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    int nb_atoms = nx * ny * nz;
    std::cout << "Simulation time for " << nb_atoms << " atoms: " << elapsed.count() << " seconds.\n";
}

int main() {
    double total_time = 100.0;  // Total simulation time
    double timestep = 0.001;  // Timestep for the simulation

    measure_performance(1,1,1,total_time,timestep,0.8);   // 1 atom
    measure_performance(5, 5, 5, total_time, timestep, 0.8);  // 125 atoms
    measure_performance(6, 6, 6, total_time, timestep, 0.8);  // 216 atoms
    measure_performance(8, 8, 8, total_time, timestep, 0.8);  // 512 atoms
    measure_performance(10, 10, 10, total_time, timestep, 0.8);  // 1000 atoms
    return 0;
}
