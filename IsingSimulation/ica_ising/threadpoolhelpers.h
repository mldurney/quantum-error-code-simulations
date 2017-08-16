#ifndef THREADPOOLHELPERS_H_
#define THREADPOOLHELPERS_H_

#include <thread>

namespace ising {
template <typename Iter>
void runInPool(Iter begin, Iter end, unsigned int threadCount) {
    ThreadPool pool(threadCount);
    for (; begin != end; begin = std::next(begin)) pool.enqueue(*begin);
}

unsigned int getMaxThreads() {
    unsigned int maxThreads = std::thread::hardware_concurrency();

    if (maxThreads != 0) {
        return maxThreads;
    } else {
        return 4;
    }
}

unsigned int getNumThreads(unsigned int remaining) {
    unsigned int numThreads = getMaxThreads();

    if (remaining < numThreads) {
        numThreads = remaining;
    }

    return numThreads;
}
}

#endif /* THREADPOOLHELPERS_H_ */