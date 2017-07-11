#ifndef ISINGHELPERS_H_
#define ISINGHELPERS_H_

#include <stdio.h>
#include "lattices.h"

namespace ising {
const int MAX_FILENAME_SIZE = 255;

Hamiltonian readHamiltonian(std::ifstream &file, char &shape);
std::string getOutFilename(std::string inFilename, std::string oldDir,
                           std::string newDir);
void writeOutput(std::string filename, std::vector<double> temp,
                 std::vector<double> results);
}

#endif /* ISINGHELPERS_H_ */
