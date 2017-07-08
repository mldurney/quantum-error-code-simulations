#ifndef ISINGSIMULATION_H_
#define ISINGSIMULATION_H_

#include <thread>
#include "ising.h"
#include "simulatedlattice.h"

const int PREUPDATES = 500;

void receiveInput(int argc, char* argv[], string& filename, double& t,
        double& dt, int& n, int& updates, char& mode);
int getNumThreads(int remaining);
int getMaxThreads();
void runLatticeSimulation(Lattice* lattice, int updates,
        double* magnetizationPtr, double* binderCumulantPtr);

#endif /* ISINGSIMULATION_H_ */
