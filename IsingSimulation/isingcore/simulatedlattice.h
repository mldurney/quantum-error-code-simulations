#ifndef SIMULATEDLATTICE_H_
#define SIMULATEDLATTICE_H_

#include "lattices.h"

namespace ising {
class SimulatedLattice {
   public:
    SimulatedLattice(Lattice *lattice, const int preupdates = 0);
    ~SimulatedLattice();
    Lattice *getLattice() const { return lattice; }
    void setAveMag(double mag);
    void setAveMag2(double mag2);
    void setAveMag4(double mag4);
    double getAveMag() const { return aveMag; }
    double getAveMag2() const { return aveMag2; }
    double getAveMag4() const { return aveMag4; }
    double getBinderCumulant();

   private:
    Lattice *lattice;
    double aveMag;
    double aveMag2;
    double aveMag4;
};
}

#endif /* SIMULATEDLATTICE_H_ */
