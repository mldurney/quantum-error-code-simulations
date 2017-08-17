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
    double getTemperature() { return temperature; }
    void setTemperature(double t);
    int getTotalEnergy() { return findTotalEnergy(); }
    double getMagnetization() { return findMagnetization(); }

    void update();
    void reinit() { initSpins(); }
    void flipSpins();
    void flipSpin(int index) { spins[index] *= -1; }
    void print() const;

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
    const LatticeProperties& prop;
    uint replicaIndex;
    double temperature;
    ivector randomizedIndices;
    RandomGenerator gen;
};
}

#endif /* REPLICA_H_ */