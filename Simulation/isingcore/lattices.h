#ifndef LATTICES_H_
#define LATTICES_H_

#include <stdlib.h>
#include <memory>
#include "hamiltonian.h"
#include "replica.h"

typedef std::vector<std::shared_ptr<ising::Replica>> replicavector;
typedef std::vector<replicavector> replicavector2;

namespace ising {
const uint REPLICAS = 2;

class Lattice {
   public:
    Lattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m = 'p');
    virtual ~Lattice() = default;

    std::string getType() const { return prop.type; }
    char getShape() const { return prop.shape; }
    ldouble getMinTemperature() const { return prop.minT; }
    ldouble getChangeTemperature() const { return prop.dT; }
    int getNumTemperatures() const { return prop.numT; }
    char getMode() const { return prop.mode; }
    int getSize() const { return prop.size; }
    int getRows() const { return prop.rows; }
    int getCols() const { return prop.cols; }
    ldouble getJTemperature() const { return jTemperature; }

    const Hamiltonian& getHamiltonian() const { return prop.hamiltonian; }
    const ivector2& getHFunction() const { return prop.hFunction; }
    const ivector& getIndices() const { return prop.indices; }
    int getNumIndices() const { return prop.numIndices; }
    const i2arrayvector& getLocations() const { return prop.locations; }
    const ivector3& getIndInteractions() const { return prop.indInteractions; }
    const ivector2& getXDisplacements() const { return prop.xDisplacements; }
    const ivector2& getYDisplacements() const { return prop.yDisplacements; }
    const dvector2& getDistances() const { return prop.distances; }

    const LatticeProperties& getProperties() const { return prop; }
    const replicavector2& getConfigs() const { return configs; }
    const ivector getReplicaIndices() const { return replicaIndices; }
    Replica getReplicaCopy(unsigned i, unsigned j = 0);

    void monteCarloSweep();
    void houdayerClusterMove();
    void houdayerClusterMove(uint index);
    void parallelTemperingUpdate();
    void HCA();
    void ICA();

    virtual int findXDisplacement(int, int) { return 1; };
    virtual int findYDisplacement(int, int) { return 1; };
    virtual ldouble findDistance(int, int) { return 1; }
    void switchMode(char m);
    void setTemperature(ldouble t);
    void setJTemperature(ldouble t) { jTemperature = t; }

   protected:
    void setType(std::string t) { prop.type = t; }
    void setMode(char m) { prop.mode = m; }
    void setSize(int s) { prop.size = s; }
    void setRows(int r) { prop.rows = r; }
    void setCols(int c) { prop.cols = c; }

    virtual void checkShape() const {}
    void shapeError() const;
    void generateDistances();

    LatticeProperties prop;
    RandomGenerator gen;

   private:
    replicavector2 configs;
    ivector replicaIndices;
    ldouble jTemperature;

    void mapsToSequences();
    void swapConfigs(uint i, uint j);
};

class RectangularLattice : public virtual Lattice {
   public:
    RectangularLattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m,
                       int r = -1, int c = -1);
    int findXDisplacement(int i, int j);
    int findYDisplacement(int i, int j);
    ldouble findDistance(int i, int j);

   protected:
    void checkShape() const;

   private:
    void guessRowsCols();
};

class SquareLattice : public RectangularLattice {
   public:
    SquareLattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m = 'p',
                  int s = -1);
    int findXDisplacement(int i, int j) {
        return RectangularLattice::findXDisplacement(i, j);
    }
    int findYDisplacement(int i, int j) {
        return RectangularLattice::findYDisplacement(i, j);
    }
    ldouble findDistance(int i, int j) {
        return RectangularLattice::findDistance(i, j);
    }

   protected:
    void checkShape() const;

   private:
    void guessSide();
};

class TriangularLattice : public virtual Lattice {
   public:
    TriangularLattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m = 'p',
                      int r = -1, int c = -1);
    int findXDisplacement(int i, int j);
    int findYDisplacement(int i, int j);
    ldouble findDistance(int i, int j);

   protected:
    void checkShape() const;

   private:
    void guessRowsCols();
};

class STriangularLattice : public TriangularLattice {
   public:
    STriangularLattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m = 'p',
                       int s = -1);
    int findXDisplacement(int i, int j) {
        return TriangularLattice::findXDisplacement(i, j);
    }
    int findYDisplacement(int i, int j) {
        return TriangularLattice::findYDisplacement(i, j);
    }
    ldouble findDistance(int i, int j) {
        return TriangularLattice::findDistance(i, j);
    }

   protected:
    void checkShape() const;

   private:
    void guessSide();
};
}

#endif /* LATTICES_H_ */
