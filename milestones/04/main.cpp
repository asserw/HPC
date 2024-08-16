//
// Created by asser on 6/18/2024.
//
#include "../header/atoms.h"
#include "../header/compute_kinetic_energy.h"
#include "../header/lj_direct_summation.h"
#include "../header/verlet.h"
#include "../header/xyz.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

void run_simulation(double timestep, double total_time, const std::string &output_file, std::vector<double> &energy_deviation) {
    auto [names, positions, velocities] = read_xyz_with_velocities("lj54.xyz");

    Atoms atoms(positions, velocities);
    atoms.masses.setOnes();

    int num_steps = static_cast<int>(total_time / timestep);

    std::ofstream traj(output_file);
    std::vector<double> total_energies(num_steps);

    double potential_energy_i = lj_direct_summation(atoms, 1.0, 1.0);
    double kinetic_energy_i = compute_kinetic_energies(atoms);
    double initial_total_energy = kinetic_energy_i + potential_energy_i;
    double tolerance = 0.2; // < 0.2% of initial_total_energy

    std::cout << "Initial Total Energy: " << initial_total_energy << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        verlet_step1(atoms, timestep);
        double potential_energy = lj_direct_summation(atoms, 1.0, 1.0);
        verlet_step2(atoms, timestep);

        // Compute kinetic energy
        double kinetic_energy = compute_kinetic_energies(atoms);

        double total_energy = kinetic_energy + potential_energy; // energy in eV
        total_energies[step] = total_energy;
        std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;

        // Store energy deviation
        energy_deviation.push_back(std::abs(total_energy - initial_total_energy));

        // Check energy conservation
        if (std::abs(total_energy - initial_total_energy) > tolerance) {
            std::cerr << "Warning: Energy difference exceeded tolerance at step " << step << std::endl;
        }

        if (step % 100 == 0) {
            write_xyz(traj, atoms);
        }
    }
    traj.close();

    // Write total energy to a file for plotting
    std::ofstream energy_file("energy_" + std::to_string(timestep) + ".txt");
    for (const auto &energy : total_energies) {
        energy_file << energy << "\n";
    }
    energy_file.close();
}

int main() {
    double total_time = 100.0; // time in ps
    std::vector<double> energy_deviation_001, energy_deviation_004, energy_deviation_008, energy_deviation_010;

    run_simulation(0.001, total_time, "traj_001.xyz", energy_deviation_001); // timestep in ps
    run_simulation(0.004, total_time, "traj_004.xyz", energy_deviation_004);
    run_simulation(0.008, total_time, "traj_008.xyz", energy_deviation_008);
    run_simulation(0.010, total_time, "traj_010.xyz", energy_deviation_010);
    // Analysis: Print the maximum energy deviation for each time step
    auto max_dev_001 = *std::max_element(energy_deviation_001.begin(), energy_deviation_001.end());
    auto max_dev_004 = *std::max_element(energy_deviation_004.begin(), energy_deviation_004.end());
    auto max_dev_008 = *std::max_element(energy_deviation_008.begin(), energy_deviation_008.end());
    auto max_dev_010 = *std::max_element(energy_deviation_010.begin(), energy_deviation_010.end());

    std::cout << "Max energy deviation for time step 0.001: " << max_dev_001 << std::endl;
    std::cout << "Max energy deviation for time step 0.004: " << max_dev_004 << std::endl;
    std::cout << "Max energy deviation for time step 0.008: " << max_dev_008 << std::endl;
    std::cout << "Max energy deviation for time step 0.010: " << max_dev_010 << std::endl;

    return 0;
}