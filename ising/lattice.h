#ifndef LATTICE_H_
#define LATTICE_H_

#include <math.h>

class Lattice {
public:
    Lattice(int n, float t);
    Lattice(int r, int c, float t);
    void updateLattice();
    float findTotalEnergy();
    float findIndexEnergy(int x, int y);
    float findMagnetism();
    void printLattice();
    ~Lattice();

private:
    int** board;
    float** probabilities;
    int rows;
    int cols;
    float temp;
    static const int MAX_DISTANCE = 1;
    void findProbabilities();
    double distance(int x, int y) {return sqrt(pow(x, 2) + pow(y, 2));}
};

#endif /* LATTICE_H_ */
