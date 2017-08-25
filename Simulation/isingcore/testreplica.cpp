#include <cassert>
#include <iostream>
#include "lattices.h"

using namespace ising;

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

    char shape = '\0';
    int rows = -1;
    int cols = -1;
    char c = (char)file.peek();

    if (isalpha(c)) {
        std::string line;
        getline(file, line);
        std::istringstream lineStream(line);

        lineStream >> shape;
        lineStream.ignore();
        lineStream >> rows;
        lineStream.ignore();
        lineStream >> cols;
    }

    Hamiltonian h =
        Hamiltonian(importHamiltonianVector(file), shape, rows, cols);

    file.close();

    double t = 1;
    double dt = .1;
    uint n = 10;
    char m = 'p';
    int num = 5;

    Lattice lattice(h, t, dt, n, m);
    Replica replica(lattice.getProperties(), num);

    std::cout << std::endl;

    // Test that properties match input values

    auto prop = replica.getProperties();
    std::cout << "Lattice type: " << prop.type << std::endl;
    assert(prop.shape == shape && "Shape incongruent!\n");
    assert(prop.minT == t && "Min temperature incongurent!\n");
    assert(prop.dT == dt && "Change in temperature incongruent! \n");
    assert(prop.numT == n && "Number of temperatures incongruent!\n");
    assert(prop.mode == m && "Mode incongruent!\n");
    std::cout << "Lattice size: " << prop.size << std::endl;
    assert(prop.rows == h.getRows() && "Rows incongruent!\n");
    assert(prop.cols == h.getCols() && "Cols incongruent!\n");
    assert(prop.numIndices == h.getNumIndices() &&
           "Num of indices incongruent!\n");

    auto h1 = prop.hamiltonian.getHamiltonian(), h2 = h.getHamiltonian();
    for (auto it1 = h1.begin(), it2 = h2.begin();
         it1 != h1.end() && it2 != h2.end(); ++it1, ++it2) {
        for (auto it11 = it1->begin(), it22 = it2->begin();
             it11 != it1->end() && it22 != it2->end(); ++it11, ++it22) {
            assert(*it11 == *it22 && "Hamiltonian incongruent!\n");
        }
    }

    std::cout << std::endl;

    // Test replica member functions

    assert(replica.getTemperature() == t + num * dt &&
           "Replica temperature incorrect!\n");
    double newT = (replica.getTemperature() + t) / 2;
    replica.setTemperature(newT);
    assert(replica.getTemperature() == newT &&
           "New replica temperature incorrect!\n");

    std::cout << "Printing spins:" << std::endl;
    auto spins = replica.getSpins();
    for (auto &s : spins) {
        std::cout << (int)s;
    }

    std::cout << std::endl << std::endl;

    replica.print();
    std::cout << std::endl;

    replica.flipSpins();
    auto flipped = replica.getSpins();
    for (int i = (int)flipped.size() - 1; i >= 0; --i) {
        assert(flipped[i] == -1 * spins[i] && "Spins not flipped!\n");
    }

    replica.flipSpin(0);
    flipped = replica.getSpins();
    assert(flipped[0] == spins[0] && "Spin not flipped!\n");

    std::cout << "Printing lattice:" << std::endl;
    replica.print();
    std::cout << std::endl;

    std::cout << "Printing updated lattice:" << std::endl;
    replica.update();
    replica.print();
    std::cout << std::endl;

    std::cout << "Printing reinitialized lattice:" << std::endl;
    replica.reinit();
    replica.print();
    std::cout << std::endl;

    std::cout << "Total energy: " << replica.getTotalEnergy() << std::endl;
    std::cout << "Magnetization: " << replica.getMagnetization() << std::endl;
    std::cout << std::endl;

    uint updates = 1000;
    std::cout << "Printing lattice after " << updates
              << " updates:" << std::endl;
    for (unsigned i = 0; i < updates; ++i) {
        replica.update();
    }
    replica.print();
    std::cout << std::endl;
    std::cout << "Total energy: " << replica.getTotalEnergy() << std::endl;
    std::cout << "Magnetization: " << replica.getMagnetization() << std::endl;

    std::cout << std::endl;
}