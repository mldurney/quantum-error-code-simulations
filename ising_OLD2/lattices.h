#ifndef LATTICES_H_
#define LATTICES_H_

#include <array>
#include <vector>

using namespace std;

#define E 2.71828182845904523536
enum {ALL = 'a', RANDOM = 'r'};

class Lattice
{
public:
    Lattice(double t, char m, int c);
    Lattice(double t, char m, vector<int> c, bool constCoupling);
    void updateLattice();
    void switchMode(char mode);
    double findTotalEnergy();
    double findMagnetism();
    void printLattice();
    virtual ~Lattice() = default;

protected:
    double temp;
    char mode;
    int numIndices;
    int coupling = 1;
    const bool isCouplingConst;
    vector<int> couplings;
    vector<vector<int>> indices;
    vector<int> spins;
    vector<vector<int>> localTerms;
    vector<vector<array<int, 2>>> hamiltonian;
    void initSpins();
    virtual void initLocalTerms() = 0;
    void initHamiltonian();
    void updateAll();
    void updateRandom();
    double findProbability(int index);
    double findIndexEnergy(int index);
    double findIndexEnergyQuick(int index);
    double (Lattice::*findIndexEnergyPtr)(int);
};


class RectangularLattice : public Lattice
{
public:
    RectangularLattice(int row, int col, double t, char m, int c);
    RectangularLattice(int row, int col, double t, char m, vector<int> c,
            bool constCoupling);

protected:
    int rows;
    int cols;
    void initLocalTerms();
};


class SquareLattice : public RectangularLattice
{
public:
    SquareLattice(int side, double t, char m, int c);
    SquareLattice(int side, double t, char m, vector<int> c);
};


class TriangularLattice : public Lattice
{
public:
    TriangularLattice(int row, int col, double t, char m, int c);
    TriangularLattice(int row, int col, double t, char m, vector<int> c,
            bool constCoupling);

protected:
    int rows;
    int cols;
    void initLocalTerms();
};

#endif /* LATTICES_H_ */
