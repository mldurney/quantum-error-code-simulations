#ifndef ISING_H_
#define ISING_H_

#include "isinghelpers.h"

namespace ising {
std::string receiveHamiltonianFile(int argc, char *argv[]);
void runReplica(Replica replica, unsigned int preupdates = 0);
void runReplicaHelp();
}

#endif /* ISING_H_ */
