//
// Created by asser on 6/19/2024.
//
#include "../header/lj_direct_summation_neighbors.h"
#include <cmath>

double lj_direct_summation_neighbors(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma, double cutoff) {
    double potential_energy = 0.0;
    atoms.forces.setZero();

    double shift = 4 * epsilon * (std::pow(sigma / cutoff, 12) - std::pow(sigma / cutoff, 6));

    neighbor_list.update(atoms, cutoff);

    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            // Calculate the distance vector between atoms i and j
            Eigen::Array3d r_ij = atoms.positions.col(i) - atoms.positions.col(j);
            double r = r_ij.matrix().norm(); // Magnitude of the distance vector

            if (r < cutoff) {
                // Calculate Lennard-Jones potential components
                double sr = sigma / r;
                double sr6 = std::pow(sr, 6);
                double sr12 = std::pow(sr, 12);

                // Calculate the Lennard-Jones potential energy with the shift
                double lj_potential = 4 * epsilon * (sr12 - sr6) - shift;
                potential_energy += lj_potential;

                // Calculate the magnitude of the force
                double force_magnitude = 24 * epsilon * (2 * sr12 - sr6) / r;
                Eigen::Array3d force = force_magnitude * (r_ij / r);

                // Update forces for atoms i and j
                atoms.forces.col(i) += force;
                atoms.forces.col(j) -= force;
            }
        }
    }

    return potential_energy;
}
