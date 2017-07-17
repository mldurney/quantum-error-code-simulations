#ifndef LATTICES_H_
#define LATTICES_H_

#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "hamiltonian.h"

namespace ising {
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

    virtual void checkShape() const {}
    void shapeError() const;
    void setMode(char m) { mode = m; };

    void updateAll();
    void updatePseudo();
    void updateRandom();
    double findProbability(int index);
    int findTotalEnergy();
    int findIndexEnergy(int index);
    double findMagnetism();

   private:
    const char shape;
    const double temp;
    char mode;

    void initSpins();
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
    int findIndexEnergy(int index);
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

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class RectangularLatticeFast : public RectangularLattice, public LatticeFast {
   public:
    RectangularLatticeFast(Hamiltonian h, double t, char m = 'p', int r = -1,
                           int c = -1)
        : RectangularLattice(h, t, m, r, c),
          LatticeFast(h, t, m),
          Lattice(h, t, m){};

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
        : SquareLattice(h, t, m, s), LatticeFast(h, t, m), Lattice(h, t, m){};

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

   private:
    int rows;
    int cols;

    void guessRowsCols();
};

class TriangularLatticeFast : public TriangularLattice, public LatticeFast {
   public:
    TriangularLatticeFast(Hamiltonian h, double t, char m = 'p', int r = -1,
                          int c = -1)
        : TriangularLattice(h, t, m, r, c),
          LatticeFast(h, t, m),
          Lattice(h, t, m){};

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
        : STriangularLattice(h, t, m, s),
          LatticeFast(h, t, m),
          Lattice(h, t, m){};

   protected:
    void checkShape() const { STriangularLattice::checkShape(); }
};
}

#endif /* LATTICES_H_ */
