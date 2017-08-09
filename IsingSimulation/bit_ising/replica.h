#ifndef REPLCIA_H_
#define REPLICA_H_

#include "common.h"

namespace ising {
class Replica {
   public:
    Replica();
    ~Replica() {}

    cvector getSpins() const { return spins; }
    Lattice::LatticeProperties getProperties() const { return properties; }

    void metropolisUpdate();
    int getTotalEnergy() { return findTotalEnergy(); }
    double getMagnetization() { return findMagnetization(); }
    void reinit() { initSpins(); }

   protected:
   private:
    void updateAll();
    void updatePseudo();
    void updateRandom();
    double findProbability(int index);
    int findTotalEnergy();
    inline int findIndexEnergy(int index);
    double findMagnetization();

    void initSpins();

    cvector spins;
    Lattice::LatticeProperties properties;
};
}

#endif /* REPLICA_H_ */