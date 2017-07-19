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
enum { PLUS = 0, MINUS = 1 };
enum { ALL = 'a', PSEUDO = 'p', RANDOM = 'r' };
enum { RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't', STRIANGLE = 'v' };
const double E = 2.71828182845904523536;

class Lattice {
   public:
    Lattice(Hamiltonian h, double t, char m = 'p');
    virtual ~Lattice() = default;

    char getShape() const { return shape; }
    double getTemp() const { return temp; }
    char getMode() const { return mode; }

    void updateLattice();
    void switchMode(char m);
    int getTotalEnergy() { return findTotalEnergy(); }
    double getMagnetism() { return findMagnetism(); }
    virtual void printLattice(int cols = -1) const;

   protected:
    const Hamiltonian hamiltonian;
    const std::vector<std::vector<int>> hFunction;
    const std::vector<int> indices;
    std::vector<int> randomizedIndices;
    const int numIndices;
    const std::map<int, std::vector<int>> localTerms;
    const std::map<int, std::vector<std::vector<int>>> indInteractions;
    std::vector<int> spins;
	const std::map<int, std::map<int, double>> distances;

    virtual void checkShape() const {}
    void shapeError() const;
    void setMode(char m) { mode = m; };
	std::map<int, std::map<int, double>> generateDistances() {}
	virtual double findDistance(int i, int j) { return 1; }

    void updateAll();
    void updatePseudo();
    void updateRandom();
    double findProbability(int index);
    int findTotalEnergy();
    inline int findIndexEnergy(int index);
    double findMagnetism();

    float asFloat(uint32_t i);
    float randFloatCO();
    inline unsigned int MWC() { return (zNew() << 16) + wNew(); }

   private:
    const char shape;
    const double temp;
    char mode;
    unsigned int zSeed;
    unsigned int wSeed;

    void initSpins();

    inline unsigned int zNew() {
        return zSeed = 36969 * (zSeed & 65535) + (zSeed >> 16);
    }
    inline unsigned int wNew() {
        return wSeed = 18000 * (wSeed & 65535) + (wSeed >> 16);
    }
};

class LatticeFast : public virtual Lattice {
   public:
    LatticeFast(Hamiltonian h, double t, char m = 'p')
        : Lattice(h, t, m), coupling((char)h.getHamiltonian()[0][0]){};
    char getCoupling() const { return coupling; }

   protected:
    void updateAll();
    void updateRandom();
    int findTotalEnergy();
    inline int findIndexEnergy(int index);
    double findMagnetism();

   private:
    const char coupling;
};

class RectangularLattice : public virtual Lattice {
   public:
    RectangularLattice(Hamiltonian h, double t, char m = 'p', int r = -1,
                       int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice() const { Lattice::printLattice(getCols()); }

   protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }
    void checkShape() const;
	double findDistance(int i, int j);

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class RectangularLatticeFast : public RectangularLattice, public LatticeFast {
   public:
    RectangularLatticeFast(Hamiltonian h, double t, char m = 'p', int r = -1,
                           int c = -1)
        : Lattice(h, t, m),
          RectangularLattice(h, t, m, r, c),
          LatticeFast(h, t, m){};

   protected:
    void checkShape() const { RectangularLattice::checkShape(); }
};

class SquareLattice : public RectangularLattice {
   public:
    SquareLattice(Hamiltonian h, double t, char m = 'p', int s = -1);
    int getSide() const { return side; }
    void printLattice() const { Lattice::printLattice(getSide()); }

   protected:
    void setSide(int s) { side = s; }
    void checkShape() const;

   private:
    int side;

    void guessSide();
};

class SquareLatticeFast : public SquareLattice, public LatticeFast {
   public:
    SquareLatticeFast(Hamiltonian h, double t, char m = 'p', int s = -1)
        : Lattice(h, t, m), SquareLattice(h, t, m, s), LatticeFast(h, t, m){};

   protected:
    void checkShape() const { SquareLattice::checkShape(); }
};

class TriangularLattice : public virtual Lattice {
   public:
    TriangularLattice(Hamiltonian h, double t, char m = 'p', int r = -1,
                      int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice() const { Lattice::printLattice(getCols()); }

   protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }
    void checkShape() const;
	double findDistance(int i, int j);

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class TriangularLatticeFast : public TriangularLattice, public LatticeFast {
   public:
    TriangularLatticeFast(Hamiltonian h, double t, char m = 'p', int r = -1,
                          int c = -1)
        : Lattice(h, t, m),
          TriangularLattice(h, t, m, r, c),
          LatticeFast(h, t, m){};

   protected:
    void checkShape() const { TriangularLattice::checkShape(); }
};

class STriangularLattice : public TriangularLattice {
   public:
    STriangularLattice(Hamiltonian h, double t, char m = 'p', int s = -1);
    int getSide() const { return side; }
    void printLattice() const { Lattice::printLattice(getSide()); }

   protected:
    void setSide(int s) { side = s; }
    void checkShape() const;

   private:
    int side;

    void guessSide();
};

class STriangularLatticeFast : public STriangularLattice, public LatticeFast {
   public:
    STriangularLatticeFast(Hamiltonian h, double t, char m = 'p', int s = -1)
        : Lattice(h, t, m),
          STriangularLattice(h, t, m, s),
          LatticeFast(h, t, m){};

   protected:
    void checkShape() const { STriangularLattice::checkShape(); }
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