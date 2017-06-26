#ifndef LATTICES_H_
#define LATTICES_H_

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include "hamiltonian.h"

#define E 2.71828182845904523536
enum {ALL = 'a', RANDOM = 'r'};
enum {RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't'};

class Lattice
{
public:
    Lattice(Hamiltonian h, double t, char m);
    void updateLattice();
    void switchMode(char m);
    double findTotalEnergy();
    double findMagnetism();
    void printLattice(int cols);
    ~Lattice() {}

protected:
    Hamiltonian hamiltonian;
    vector< vector<int> > hFunction;
    vector<int> indices;
    int numIndices;
    vector< vector<int> > localTerms;
    char shape;
    double temp;
    char mode;
    vector<int> spins;

    void shapeError() const;
    void updateAll();
    void updateRandom();
    double findProbability(int index);
    double findIndexEnergy(int index);

    double (Lattice::*findIndexEnergyPtr)(int);

private:
    void initSpins();
    virtual void checkShape() const = 0;
};


class RectangularLattice : public Lattice
{
public:
    RectangularLattice(Hamiltonian h, double t, char m, int r, int c);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice() { Lattice::printLattice(getCols()); }

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
    SquareLattice(Hamiltonian h, double t, char m, int s);
    int getSide() const { return side; }

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
    TriangularLattice(Hamiltonian h, double t, char m, int r, int c);
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    void printLattice() { Lattice::printLattice(getCols()); }

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
