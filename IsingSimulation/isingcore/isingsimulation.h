#ifndef ISINGSIMULATION_H_
#define ISINGSIMULATION_H_

#include <thread>
#include "isinghelpers.h"
#include "simulatedlattice.h"

namespace ising {
const int PREUPDATES = 500;

void manageSimulations(const std::string& inFilename, const double t, const double dt, const int n,
                       const int updates, const char mode);
void receiveInputCMD(int argc, char *argv[], std::string &filename, double &t,
                     double &dt, int &n, int &updates, char &mode);
int getNumThreads(const int remaining);
int getMaxThreads();
void runLatticeSimulation(Lattice *lattice, const int updates,
                          double &indMagnetization, double &indBinderCumulant);
}

#endif /* ISINGSIMULATION_H_ */