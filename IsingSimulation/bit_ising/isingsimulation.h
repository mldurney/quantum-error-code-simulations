#ifndef ISINGSIMULATION_H_
#define ISINGSIMULATION_H_

#include <thread>
#include "isinghelpers.h"
#include "simulatedlattice.h"
#include "threadpool.h"

namespace ising {

void manageSimulations(const std::string &inFilename, const double t,
                       const double dt, const int n, const int updates,
                       const int trials, const char mode);
void receiveInputCMD(int argc, char *argv[], std::string &filename, double &t,
                     double &dt, int &n, int &updates, int &trials, char &mode);
int getMaxThreads();
int getNumThreads(const int remaining);

template <typename Iter>
void runInPool(Iter begin, Iter end, int threadCount);
}

#endif /* ISINGSIMULATION_H_ */
