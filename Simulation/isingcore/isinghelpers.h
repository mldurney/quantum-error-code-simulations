#ifndef ISINGHELPERS_H_
#define ISINGHELPERS_H_

#include <stdio.h>
#include <experimental/filesystem>
#include "lattices.h"

namespace fs = std::experimental::filesystem;

namespace ising {
const int MAX_FILENAME_SIZE = 255;

Hamiltonian readHamiltonian(std::ifstream& file, char& shape);
Lattice* chooseLattice(char shape, const Hamiltonian& hamiltonian, double t,
                       double dt, uint n, char m);
std::string getOutFilename(const std::string& inFilename,
                           const std::string& newDir);
std::string getOutFilename(const std::string& inFilename,
                           const std::string& oldDir,
                           const std::string& newDir);
void writeOutput(const std::string& filename, const dmap& temperatures,
                 const dmap& results);
}

#endif /* ISINGHELPERS_H_ */
