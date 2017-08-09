#ifndef LATTICES_H_
#define LATTICES_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <numeric>
#include <random>
#include <utility>
#include "hamiltonian.h"

namespace ising {
const int REPLICAS = 2;

class Lattice {
   public:
    Lattice(Hamiltonian h, double t, char m = 'p', bool init = true);
    virtual ~Lattice() = default;

    std::string getType() const { return type; }
    char getShape() const { return shape; }
    double getTemp() const { return temp; }
    char getMode() const { return mode; }
    int getSize() const { return size; }
    const Hamiltonian getHamiltonian() const { return hamiltonian; }
    const ivector2 getHFunction() const { return hFunction; }
    const ivector getIndices() const { return indices; }
    const int getNumIndices() const { return numIndices; }
    const i2arrayvector getLocations() const { return locations; }
    const ivector3 getIndInteractions() const { return indInteractions; }
    cvector2 getReplicas() const { return spinReplicas; }
    cvector getSpins(int replicaNum) const { return spinReplicas[replicaNum]; }
    ivector2 getXDisplacements() const { return xDisplacements; }
    ivector2 getYDisplacements() const { return yDisplacements; }
    dvector2 getDistances() const { return distances; }

    void metropolisUpdate(int replicaNum);
    void switchMode(char m);
    int getTotalEnergy(int replicaNum) { return findTotalEnergy(replicaNum); }
    double getMagnetization(int replicaNum) {
        return findMagnetization(replicaNum);
    }

    virtual int findXDisplacement(int, int) { return 1; };
    virtual int findYDisplacement(int, int) { return 1; };
    virtual double findDistance(int, int) { return 1; }
    void reinit() { initSpins(); }
    void flipSpins(int replicaNum);
    virtual void print(int replicaNum, int cols = -1) const;

   protected:
    virtual void checkShape() const {}
    void shapeError() const;
    void setType(std::string t) { type = t; }
    void setMode(char m) { mode = m; }
    void setSize(int s) { size = s; }
    void generateDistances();

    void updateAll(int replicaNum);
    void updatePseudo(int replicaNum);
    void updateRandom(int replicaNum);
    double findProbability(int replicaNum, int index);
    int findTotalEnergy(int replicaNum);
    inline int findIndexEnergy(int replicaNum, int index);
    double findMagnetization(int replicaNum);

    float asFloat(uint32_t i);
    float randFloatCO();
    inline unsigned int MWC() { return (zNew() << 16) + wNew(); }

    class LatticeProperties {
       protected:
        const Hamiltonian hamiltonian;
        ivector2 hFunction;
        ivector indices;
        ivector randomizedIndices;
        const int numIndices;
        i2arrayvector locations;
        ivector3 indInteractions;
        cvector2 spinReplicas;
        ivector2 xDisplacements;
        ivector2 yDisplacements;
        dvector2 distances;

       private:
        std::string type;
        const char shape;
        const double temp;
        char mode;
        int size;
    } properties;

   private:
    unsigned int zSeed;
    unsigned int wSeed;

    void mapsToSequences();
    void initSpins();

    inline unsigned int zNew() {
        return zSeed = 36969 * (zSeed & 65535) + (zSeed >> 16);
    }
    inline unsigned int wNew() {
        return wSeed = 18000 * (wSeed & 65535) + (wSeed >> 16);
    }
};

class RectangularLattice : public virtual Lattice {
   public:
    RectangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
                       int r = -1, int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void print() const { Lattice::print(getCols()); }
    int findXDisplacement(int i, int j);
    int findYDisplacement(int i, int j);
    double findDistance(int i, int j);

   protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }
    void checkShape() const;

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class SquareLattice : public RectangularLattice {
   public:
    SquareLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
                  int s = -1);
    int getSide() const { return side; }
    void print() const { Lattice::print(getSide()); }
    int findXDisplacement(int i, int j) {
        return RectangularLattice::findXDisplacement(i, j);
    }
    int findYDisplacement(int i, int j) {
        return RectangularLattice::findYDisplacement(i, j);
    }
    double findDistance(int i, int j) {
        return RectangularLattice::findDistance(i, j);
    }

   protected:
    void setSide(int s) { side = s; }
    void checkShape() const;

   private:
    int side;

    void guessSide();
};

class TriangularLattice : public virtual Lattice {
   public:
    TriangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
                      int r = -1, int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void print() const { Lattice::print(getCols()); }
    int findXDisplacement(int i, int j);
    int findYDisplacement(int i, int j);
    double findDistance(int i, int j);

   protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }
    void checkShape() const;

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class STriangularLattice : public TriangularLattice {
   public:
    STriangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
                       int s = -1);
    int getSide() const { return side; }
    void print() const { Lattice::print(getSide()); }
    int findXDisplacement(int i, int j) {
        return TriangularLattice::findXDisplacement(i, j);
    }
    int findYDisplacement(int i, int j) {
        return TriangularLattice::findYDisplacement(i, j);
    }
    double findDistance(int i, int j) {
        return TriangularLattice::findDistance(i, j);
    }

   protected:
    void setSide(int s) { side = s; }
    void checkShape() const;

   private:
    int side;

    void guessSide();
};
}

#endif /* LATTICES_H_ */

/**
        Credits to Andy Gainey
   (https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers)
        and George Marsaglia (http://www.cse.yorku.ca/~oz/marsaglia-rng.html)
        for algoritms behind randFloatCO and MWC, respectively, used in random
        number generation.
*/

/**
        Variables for better random generator:
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_int_distribution<> randSpin;
        std::uniform_real_distribution<> randProb;
        std::uniform_int_distribution<> randInd;
*/