#ifndef ISINGSIMULATION_H_
#define ISINGSIMULATION_H_

#include <thread>
#include "simulation.h"

namespace ising {
void receiveSimulationInput(int argc, char *argv[], std::string &filename,
                            double &t, double &dt, int &n, int &updates,
                            int &trials, char &mode);
void manageSimulation(const std::string &inFilename, const double t,
                      const double dt, const int n, const int updates,
                      const int trials, const char mode);
}

#endif /* ISINGSIMULATION_H_ */
