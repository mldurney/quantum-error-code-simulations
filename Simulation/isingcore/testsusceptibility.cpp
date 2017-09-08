#include <cassert>
#include <iostream>
#include "isinghelpers.h"
#include "lattices.h"

using namespace ising;

void printSusceptibility(Lattice *lattice) {
    double q = 2 * ising::PI / lattice->getSize();
    imap3 runningCorr;

    auto configs = lattice->getConfigs();
    auto displacements = lattice->getXDisplacements();
    auto indices = lattice->getIndices();

    for (auto &replicas : configs) {
        uint index = replicas[0]->getReplicaIndex();

        auto spins = replicas[0]->getSpins();
        for (auto &i : indices) {
            for (auto &j : indices) {
                runningCorr[index][i][j] += spins[i] * spins[j];
            }
        }
    }

    cdmap sumCorrK0, sumCorrKq;
    for (auto &replicas : configs) {
        for (auto &i : indices) {
            for (auto &j : indices) {
                uint index = replicas[0]->getReplicaIndex();
                sumCorrK0[index] += runningCorr[index][i][j];

                if (lattice->getProperties().rows % 2 == 0 &&
                    displacements[i][j] == lattice->getProperties().rows / 2) {
                    sumCorrKq[index] +=
                        ((cdouble(runningCorr[index][i][j]) *
                          std::exp(cdouble(0, q * displacements[i][j]))) +
                         (cdouble(runningCorr[index][i][j]) *
                          std::exp(cdouble(0, q * -displacements[i][j])))) /
                        cdouble(2);
                } else {
                    sumCorrKq[index] +=
                        cdouble(runningCorr[index][i][j]) *
                        std::exp(cdouble(0, q * displacements[i][j]));
                }
            }
        }
    }

    for (auto &replicas : configs) {
        uint index = replicas[0]->getReplicaIndex();
        sumCorrK0[index] /= cdouble(lattice->getNumIndices());
        sumCorrKq[index] /= cdouble(lattice->getNumIndices());

        std::cout << "k = 0:\t" << sumCorrK0[index] << std::endl;
        std::cout << "k = q:\t" << sumCorrKq[index] << std::endl;
    }

    std::cout << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s name_of_hamiltonian_file\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    std::ifstream file(argv[1]);

    if (!file) {
        printf("Invalid file name. %s does not exist!\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    char shape;
    Hamiltonian h = readHamiltonian(file, shape);

    double t = .1;
    double dt = 0;
    uint n = 1;
    char m = 'p';

    Lattice *lattice = chooseLattice(shape, h, t, dt, n, m);

    for (uint i = 0; i < 500; ++i) {
        lattice->monteCarloSweep();
    }

    std::cout << std::endl;
    std::cout << "Printing lattice:\n";
    lattice->getReplicaCopy(0).print();

    std::cout << std::endl;
    std::cout << "Priting susceptibility:\n";
    printSusceptibility(lattice);
}