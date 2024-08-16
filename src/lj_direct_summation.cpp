        //
        // Created by asser on 6/18/2024.
        //
        #include "../header/lj_direct_summation.h"
#include <cmath>

        double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
            double potential_energy = 0.0;
            atoms.forces.setZero();

            for (size_t i = 0; i < atoms.nb_atoms(); ++i) {
                for (size_t j = i + 1; j < atoms.nb_atoms(); ++j) {
                    // Calculate the distance vector between atoms i and j
                    Eigen::Array3d r_ij = atoms.positions.col(i) - atoms.positions.col(j);
                    double r = r_ij.matrix().norm(); // Magnitude of the distance vector

                    // Calculate Lennard-Jones potential components
                    double sr = sigma / r;
                    double sr6 = std::pow(sr, 6);
                    double sr12 = std::pow(sr, 12);

                    // Calculate the Lennard-Jones potential energy
                    double lj_potential = 4 * epsilon * (sr12 - sr6);
                    potential_energy += lj_potential;

                    // Calculate the magnitude of the force
                    double force_magnitude = 24 * epsilon * (2 * sr12 - sr6) / r; // Or (r*r)
                    Eigen::Array3d force = force_magnitude * (r_ij/r);            //    r_ij

                    // Update forces for atoms i and j
                    atoms.forces.col(i) += force;
                    atoms.forces.col(j) -= force;
                }
            }

            return potential_energy;
        }
