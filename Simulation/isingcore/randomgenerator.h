#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

#include <stdlib.h>
#include <time.h>

namespace ising {
class RandomGenerator {
   public:
    RandomGenerator() {
        unsigned int zFactor = (unsigned int)reinterpret_cast<intptr_t>(&zSeed);
        unsigned int wFactor = (unsigned int)reinterpret_cast<intptr_t>(&wSeed);

        srand((unsigned int)time(NULL) + (zFactor ^ wFactor) / zFactor);
        zSeed = rand() + (zFactor & rand());
        wSeed = rand() - (wFactor | rand());
    }

    inline unsigned int MWC() { return (zNew() << 16) + wNew(); }
    float randFloatCO() { return asFloat(0x3F800000U | (MWC() >> 9)) - 1.0f; }

   private:
    unsigned int zSeed;
    unsigned int wSeed;

    inline unsigned int zNew() {
        return zSeed = 36969 * (zSeed & 65535) + (zSeed >> 16);
    }

    inline unsigned int wNew() {
        return wSeed = 18000 * (wSeed & 65535) + (wSeed >> 16);
    }

    float asFloat(unsigned int i) {
        union {
            unsigned int i;
            float f;
        } pun = {i};
        return pun.f;
    }
};
}

#endif /* RANDOMGENERATOR_H_ */

/**
        Credits to Andy Gainey
   (https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers)
        and George Marsaglia (http://www.cse.yorku.ca/~oz/marsaglia-rng.html)
        for algoritms behind randFloatCO and MWC, respectively, used in random
        number generation.
*/

/**
        Variables for better random generator:
        #include <random>
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_int_distribution<> randSpin;
        std::uniform_real_distribution<> randProb;
        std::uniform_int_distribution<> randInd;
*/