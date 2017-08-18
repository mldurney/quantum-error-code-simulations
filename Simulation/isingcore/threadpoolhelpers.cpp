#include "threadpoolhelpers.h"

using namespace ising;

unsigned int ising::getMaxThreads() {
    unsigned int maxThreads = std::thread::hardware_concurrency();

    if (maxThreads != 0) {
        return maxThreads;
    } else {
        return 4;
    }
}

unsigned int ising::getNumThreads(unsigned int remaining) {
    unsigned int numThreads = getMaxThreads();

    if (remaining < numThreads) {
        numThreads = remaining;
    }

    if (numThreads > MAXTHREADS) {
        numThreads = MAXTHREADS;
    }

    return numThreads;
}