#ifndef REPLCIA_H_
#define REPLICA_H_

#include <cmath>
#include "common.h"
#include "properties.h"

namespace ising {
class Replica {
   public:
    Replica(const LatticeProperties& properties, uint n);
    ~Replica() {}

    cvector getSpins() const { return spins; }
    const LatticeProperties& getProperties() const { return prop; }
    uint getReplicaIndex() { return replicaIndex; }
    ldouble getTemperature() { return temperature; }
    void setTemperature(ldouble t);
    int getTotalEnergy() { return findTotalEnergy(); }
    ldouble getMagnetization() { return findMagnetization(); }

    void update();
    void reinit() { initSpins(); }
    void flipSpins();
    void flipSpin(int index) { spins[index] *= -1; }
    void print() const;

   private:
    void updateAll();
    void updatePseudo();
    void updateRandom();
    ldouble findProbability(int index);
    int findTotalEnergy();
    inline int findIndexEnergy(int index);
    ldouble findMagnetization();

    void initSpins();

    cvector spins;
    const LatticeProperties& prop;
    uint replicaIndex;
    ldouble temperature;
    ivector randomizedIndices;
    RandomGenerator gen;
};
}

#endif /* REPLICA_H_ */