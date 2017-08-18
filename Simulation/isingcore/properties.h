#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include "common.h"
#include "hamiltonian.h"

namespace ising {
struct LatticeProperties {
    LatticeProperties(Hamiltonian h, ldouble t, ldouble dt, uint n, char m)
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
    const uint numIndices;
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
    uint size;
    const ldouble minT;
    const ldouble dT;
    const uint numT;
    char mode;
};
}

#endif /* PROPERTIES_H_ */