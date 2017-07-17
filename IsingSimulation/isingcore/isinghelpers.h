#ifndef ISINGHELPERS_H_
#define ISINGHELPERS_H_

#include <stdio.h>
#include "lattices.h"

namespace ising {
const int MAX_FILENAME_SIZE = 255;

Hamiltonian readHamiltonian(std::ifstream& file, char& shape);
Lattice* chooseLattice(const char shape, const Hamiltonian& hamiltonian,
                       const double temp, const char mode);
std::string getOutFilename(const std::string& inFilename,
                           const std::string& oldDir,
                           const std::string& newDir);
void writeOutput(const std::string& filename, const std::vector<double>& temp,
                 const std::vector<double>& results);
}

#endif /* ISINGHELPERS_H_ */
