#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include "lattices.h"

using namespace std;

void help();

int main(int argc, char* argv[])
{
    int rows, cols, n, trials, coupling;
    double t, dt;
    char mode;
    vector<RectangularLattice*> lattices;
    vector<array<double, 2>> magnetisms;

    if (argc == 9)
    {
        rows = atoi(argv[1]);
        cols = atoi(argv[2]);
        t = atof(argv[3]);
        dt = atof(argv[4]);
        n = atoi(argv[5]);
        trials = atoi(argv[6]);
        coupling = atoi(argv[7]);
        mode = (argv[8])[0];
    }
    else if (argc == 1)
    {
        cout << "Enter number of rows: ";
        cin >> rows;
        cout << "Enter number of columns: ";
        cin >> cols;
        cout << "Enter minimum temperature (K): ";
        cin >> t;
        cout << "Enter change in temperature between lattices (K): ";
        cin >> dt;
        cout << "Enter total number of lattices: ";
        cin >> n;
        cout << "Enter number of trials per lattice: ";
        cin >> trials;
        cout << "Enter coupling between lattice pairs: ";
        cin >> coupling;
        cout << "Enter update mode (a - all, r - random): ";
        cin >> mode;
    }
    else
    {
        cout << "Usage: " << argv[0] << " size(int) min_temperature(float) ";
        cout << "change_temperature(float) num_lattices(int) trials(int) ";
        cout << "coupling(int) mode(char)\n\n";
        exit(1);
    }

    for (int i = 0; i < n; i++)
    {
        RectangularLattice *nextL =
                new RectangularLattice(rows, cols, t + i * dt, mode, coupling);
        lattices.push_back(nextL);

        array<double, 2> nextM = {nextL->findMagnetism(), t + i * dt};
        magnetisms.push_back(nextM);
    }

    vector<RectangularLattice*>::iterator itL;
    vector<array<double, 2>>::iterator itM;

    for (itL = lattices.begin(), itM = magnetisms.begin();
            itL != lattices.end() && itM != magnetisms.end(); itL++, itM++)
    {

        for (int i = 1; i < trials; i++)
        {
            (*itL)->updateLattice();
            (*itM)[0] += (*itL)->findMagnetism();
        }

        (*itM)[0] = abs((*itM)[0]) / trials;
    }

    char filename[50];
    sprintf(filename, "mag_%dx%d-%.2fT+%.2fdTx%d-%du-%dc-%c.csv",
            rows, cols, t, dt, n, trials, coupling, mode);
    ofstream file;
    file.open("magnetizations/" + string(filename));

    for (itM = magnetisms.begin(); itM != magnetisms.end(); itM++)
    {
        file << (*itM)[0] << "," << (*itM)[1] << endl;
    }

    printf("Recorded average magnetisms by temperature in %s\n", filename);

    for (itL = lattices.begin(); itL != lattices.end(); itL++)
    {
        delete *itL;
    }

    lattices.clear();
    magnetisms.clear();
}
