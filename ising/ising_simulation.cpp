/*
 * TODO:
 *      Add change tracking
 *      Comment to glory
 *      Make variables more accessible (at least read after start)
 *          Maybe even change...?
 *      Enable percentage tracking
 *      Add more parameters
 *          Like what?
 *      Add safety nets and what not
 *      Create list of boards to see progression over time
 *      Enable output of results
 *      Error checking
 *
 *      Include option to compare changes
 *
 *      Desired behavior:
 *          Plot of average value of absolute magnetization vs temp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include "lattice_near.h"

using namespace std;

void help();

int main() {
    int cols, rows, n, trials;
    double t, dt;
    char mode;
    vector<Lattice*> lattices;
    vector<array<double, 2>> magnetisms;

    cout << "Enter number of rows: ";
    cin >> rows;
    cout << "Enter number of cols: ";
    cin >> cols;
    cout << "Enter minimum temperature (K): ";
    cin >> t;
    cout << "Enter change in temperature between lattices (K): ";
    cin >> dt;
    cout << "Enter total number of lattices: ";
    cin >> n;
    cout << "Enter number of trials per lattice: ";
    cin >> trials;
    cout << "Enter update mode (a - all, r - random): ";
    cin >> mode;

    for (int i = 0; i < n; i++) {
        Lattice *nextL = new Lattice(rows, cols, t + i * dt, mode);
        lattices.push_back(nextL);
        array<double, 2> nextM = {nextL->findMagnetism(), t + i * dt};
        magnetisms.push_back(nextM);
    }
    vector<Lattice*>::iterator itL;
    vector<array<double, 2>>::iterator itM;

    for (itL = lattices.begin(), itM = magnetisms.begin();
            itL != lattices.end() && itM != magnetisms.end(); itL++, itM++) {

        for (int i = 1; i < trials; i++) {
            (*itL)->updateLattice();
            (*itM)[0] += (*itL)->findMagnetism();
        }

        (*itM)[0] /= trials;
    }

    char filename[] = {"average_magnetisms.csv"};
    ofstream file;
    file.open(filename);

    for (itM = magnetisms.begin(); itM != magnetisms.end(); itM++) {
        file << (*itM)[0] << "," << (*itM)[1] << endl;
    }

    printf("\nRecorded average magnetisms by temperature in %s\n\n", filename);

    for (itL = lattices.begin(); itL != lattices.end(); itL++) {
        delete *itL;
    }

    lattices.clear();
    magnetisms.clear();
}
