#ifndef LATTICE_H_
#define LATTICE_H_

#include <cmath>

double const E = 2.71828182845904523536;
enum {ALL = 'a', RANDOM = 'r'};

class Lattice {
public:
    Lattice(int n, float t, char m);
    Lattice(int r, int c, float t, char m);
    void updateLattice();
    void switchMode(char mode);
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
    char mode;
    static const int MAX_DISTANCE = 1;
    void updateAll();
    void updateRandom();
    void findProbabilities();
    double distance(int x, int y) {return sqrt(pow(x, 2) + pow(y, 2));}
};

#endif /* LATTICE_H_ */
