#ifndef ISINGHELPERS_H_
#define ISINGHELPERS_H_

#include <experimental/filesystem>>
#include <stdio.h>
#include "lattices.h"

namespace fs = std::experimental::filesystem;

namespace ising {
const int MAX_FILENAME_SIZE = 255;

Hamiltonian readHamiltonian(std::ifstream& file, char& shape);
Lattice* chooseLattice(const char shape, const Hamiltonian& hamiltonian,
                       const double temp, const char mode, const bool init = true);
std::string getOutFilename(const std::string& inFilename,
						   const std::string& newDir);
std::string getOutFilename(const std::string& inFilename,
                           const std::string& oldDir,
                           const std::string& newDir);
void writeOutput(const std::string& filename, const dvector& temp,
                 const dvector& results);
}

#endif /* ISINGHELPERS_H_ */
