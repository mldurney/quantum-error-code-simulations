#ifndef ISING_H_
#define ISING_H_

#include "isinghelpers.h"

string receiveFilename(int argc, char* argv[]);
Hamiltonian readHamiltonian(ifstream& file, char& shape);
void runLattice(Lattice* lattice);
void runLatticeHelp();

#endif /* ISING_H_ */
