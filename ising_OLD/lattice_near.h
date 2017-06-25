#ifndef LATTICE_NEAR_H_
#define LATTICE_NEAR_H_

#include <cmath>

double const E = 2.71828182845904523536;
enum {ALL = 'a', RANDOM = 'r'};

class Lattice {
public:
    Lattice(int n, double t, char m);
    Lattice(int r, int c, double t, char m);
    void updateLattice();
    void switchMode(char mode);
    double findTotalEnergy();
    double findIndexEnergy(int x, int y);
    double findMagnetism();
    void printLattice();
    ~Lattice();

private:
    int** board;
    int rows;
    int cols;
    double temp;
    char mode;
    static const int MAX_DISTANCE = 1;
    void updateAll();
    void updateRandom();
    double findProbability(int x, int y);
};

#endif /* LATTICE_NEAR_H_ */
