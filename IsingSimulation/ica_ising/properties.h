#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include "common.h"
#include "hamiltonian.h"

namespace ising {
struct LatticeProperties {
    LatticeProperties(Hamiltonian h, double t, double dt, unsigned int n,
                      char m)
        : hamiltonian(h),
          numIndices(h.getNumIndices()),
          shape(h.getShape()),
          rows(h.getRows()),
          cols(h.getCols()),
          minT(t),
          dT(dt),
          numT(n),
          mode(m) {}
    ~LatticeProperties() {}

    const Hamiltonian hamiltonian;
    const unsigned int numIndices;
    ivector indices;
    ivector2 hFunction;
    ivector2 localTerms;
    i2arrayvector locations;
    ivector3 indInteractions;
    ivector2 xDisplacements;
    ivector2 yDisplacements;
    dvector2 distances;

    std::string type;
    const char shape;
    int rows;
    int cols;
    unsigned int size;
    const double minT;
    const double dT;
    const unsigned int numT;
    char mode;
};
}

#endif /* PROPERTIES_H_ */