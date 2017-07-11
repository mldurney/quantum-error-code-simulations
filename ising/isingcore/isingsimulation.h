#ifndef ISINGSIMULATION_H_
#define ISINGSIMULATION_H_

#include <thread>
#include "isinghelpers.h"
#include "simulatedlattice.h"

namespace ising {
const int PREUPDATES = 500;

void manageSimulations(std::string inFilename, double t, double dt, int n,
                       int updates, char mode);
void receiveInputCMD(int argc, char *argv[], std::string &filename, double &t,
                     double &dt, int &n, int &updates, char &mode);
int getNumThreads(int remaining);
int getMaxThreads();
void runLatticeSimulation(Lattice *lattice, int updates,
                          double &indMagnetization, double &indBinderCumulant);
}

#endif /* ISINGSIMULATION_H_ */
