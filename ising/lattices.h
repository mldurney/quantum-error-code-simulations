#ifndef LATTICES_H_
#define LATTICES_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "hamiltonian.h"

#define E 2.71828182845904523536
enum {ALL = 'a', RANDOM = 'r'};
enum {RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't'};

class Lattice
{
public:
    Lattice(Hamiltonian h, double t, char m = 'r');
    virtual ~Lattice() = default;

    char getShape() const { return shape; }
    double getTemp() const { return temp; }
    char getMode() const { return mode; }
    char getCoupling() const { return coupling; }

    void updateLattice();
    void switchMode();
    int getTotalEnergy() { return (this->*findTotalEnergyPtr)(); }
    double getMagnetism() { return (this->*findMagnetismPtr)(); }
    virtual void printLattice(int cols = -1);

protected:
    Hamiltonian hamiltonian;
    vector< vector<int> > hFunction;
    vector<int> indices;
    int numIndices;
    vector< vector<int> > localTerms;
    vector<int> spins;

    void shapeError() const;
    void setMode(char m) { mode = m; };
    void setCoupling(char c) { coupling = c; }
    void chooseSpeed();
    void setSpeedNormal();
    void setSpeedFast();

    void updateAll();
    void updateAllFast();
    void updateRandom();
    void updateRandomFast();
    double findProbability(int index);
    int findTotalEnergy();
    int findTotalEnergyFast();
    int findIndexEnergy(int index);
    int findIndexEnergyFast(int index);
    double findMagnetism();
    double findMagnetismFast();

    void (Lattice::*updateAllPtr)();
    void (Lattice::*updateRandomPtr)();
    int (Lattice::*findTotalEnergyPtr)();
    int (Lattice::*findIndexEnergyPtr)(int);
    double (Lattice::*findMagnetismPtr)();

private:
    char shape;
    double temp;
    char mode;
    char coupling = '\0';

    void initSpins();
    virtual void checkShape() const {}
};


class RectangularLattice : public Lattice
{
public:
    RectangularLattice(Hamiltonian h, double t, char m = 'r',
            int r = -1, int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice( int empty = -1 ) { Lattice::printLattice(getCols()); }

protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }

private:
    int rows;
    int cols;

    void checkShape() const;
    void guessRowsCols();
};


class SquareLattice : public RectangularLattice
{
public:
    SquareLattice(Hamiltonian h, double t, char m = 'r', int s = -1);
    int getSide() const { return side; }
    void printLattice( int empty = -1 ) { Lattice::printLattice(getSide()); }

protected:
    void setSide(int s) { side = s; }

private:
    int side;

    void checkShape() const;
    void guessSide();
};


class TriangularLattice : public Lattice
{
public:
    TriangularLattice(Hamiltonian h, double t, char m = 'r',
            int r = -1, int c = -1);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice( int empty = -1 ) { Lattice::printLattice(getCols()); }

protected:
    void setRows(int r) { rows = r; }
    void setCols(int c) { cols = c; }

private:
    int rows;
    int cols;

    void checkShape() const;
    void guessRowsCols();
};

#endif /* LATTICES_H_ */
