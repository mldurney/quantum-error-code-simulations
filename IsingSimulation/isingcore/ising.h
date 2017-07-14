#ifndef ISING_H_
#define ISING_H_

#include "isinghelpers.h"

namespace ising {
std::string receiveHamiltonianFile(int argc, char *argv[]);
void runLattice(Lattice *lattice);
void runLatticeHelp();
}

#endif /* ISING_H_ */
