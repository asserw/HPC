#ifndef HPC_PROJECT_ATOMS_H
#define HPC_PROJECT_ATOMS_H

#include <Eigen/Dense>
#include <cassert>
#include <vector>

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
using Energies_t = Eigen::ArrayXd;

class Atoms {
  public:
    Masses_t masses;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Energies_t potential_energies;
    Energies_t kinetic_energies;

    // Constructor to initialize the number of atoms
    Atoms(size_t num_atoms) :
          masses(Masses_t(num_atoms)),
          positions(Positions_t(3, num_atoms)),
          velocities(Positions_t(3, num_atoms)),
          forces(Positions_t(3, num_atoms)),
          potential_energies(Energies_t(num_atoms)),
          kinetic_energies(Energies_t(num_atoms)) {
        masses.setOnes();
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        potential_energies.setZero();
        kinetic_energies.setZero();
    }

    // Constructor with initial positions
    Atoms(const Positions_t &p) :
          masses(Masses_t(p.cols())),
          positions(p),
          velocities(3, p.cols()),
          forces(3, p.cols()),
          potential_energies(Energies_t(p.cols())),
          kinetic_energies(Energies_t(p.cols())) {
        masses.setOnes();
        velocities.setZero();
        forces.setZero();
        potential_energies.setZero();
        kinetic_energies.setZero();
    }

    // Constructor with initial positions and velocities
    Atoms(const Positions_t &p, const Velocities_t &v) :
          masses(Masses_t(p.cols())),
          positions(p),
          velocities(v),
          forces(3, p.cols()),
          potential_energies(Energies_t(p.cols())),
          kinetic_energies(Energies_t(p.cols())) {
        assert(p.cols() == v.cols());
        masses.setOnes();
        forces.setZero();
        potential_energies.setZero();
        kinetic_energies.setZero();
    }
    // Constructor to initialize the number of atoms with a cubic lattice
    Atoms(int nx, int ny, int nz, double lattice_constant) :
          masses(Masses_t(nx * ny * nz)),
          positions(Positions_t(3, nx * ny * nz)),
          velocities(Positions_t(3, nx * ny * nz)),
          forces(Positions_t(3, nx * ny * nz)),
          potential_energies(Energies_t(nx * ny * nz)),
          kinetic_energies(Energies_t(nx * ny * nz)) {
        masses.setOnes();
        velocities.setRandom();
        forces.setZero();
        potential_energies.setZero();
        kinetic_energies.setZero();
        create_cubic_lattice(nx, ny, nz, lattice_constant);
    }

    // Method to create a cubic lattice
    void create_cubic_lattice(int nx, int ny, int nz, double lattice_constant) {
        int index = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    positions(0, index) = i * lattice_constant;
                    positions(1, index) = j * lattice_constant;
                    positions(2, index) = k * lattice_constant;
                    ++index;
                }
            }
        }
    }

    // Get the number of atoms
    size_t nb_atoms() const {
        return positions.cols();
    }
};

#endif // HPC_PROJECT_ATOMS_H
