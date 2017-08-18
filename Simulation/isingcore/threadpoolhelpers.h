#ifndef THREADPOOLHELPERS_H_
#define THREADPOOLHELPERS_H_

#include <thread>
#include "threadpool.h"

namespace ising {
const uint MAXTHREADS = 32;

template <typename Iter>
void runInPool(Iter begin, Iter end, unsigned int threadCount) {
    ThreadPool pool(threadCount);
    for (; begin != end; begin = std::next(begin)) pool.enqueue(*begin);
}

unsigned int getMaxThreads();
unsigned int getNumThreads(unsigned int remaining);
}

#endif /* THREADPOOLHELPERS_H_ */