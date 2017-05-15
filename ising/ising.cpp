/*
 * IDEAS:
 *      Add change tracking
 *      Comment to glory
 *      Make variables more accessible (at least read after start)
 *          Maybe even change...?
 *      Enable percentage tracking
 *      Add more parameters
 *          Like what?
 *      Add safety nets and what not
 *      Optimize in general... Very slow
 *      Explore possibilities of distance
 *      Create list of boards to see progression over time
 *      Enable output of results
 *
 *      Change probabilities:
 *          ef > ei : e^(-1/T)(ef-ei)
 *          ef < ei : 1
 *      Optimize for distance 1
 *      Include option to compare changes
 *      Pick random points L*L times
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "lattice.h"

using namespace std;

void help();

int main() {
    int cols, rows, temp;
    char *input = new char[10];

    cout << "Enter number of rows: ";
    cin >> rows;
    cout << "Enter number of cols: ";
    cin >> cols;
    cout << "Enter temperature (K): ";
    cin >> temp;

    Lattice* lattice = new Lattice(rows, cols, temp);

    help();

    while (true) {
        cin >> input;

        if (*input == 'p') {
            lattice->printLattice();
        } else if (*input == 'u') {
            int num = atoi(&input[1]);

            for (int i = 0; i <= num; i++) {
                lattice->updateLattice();
            }

            printf("\nUpdated lattice %d times.\n", num);
        } else if (*input == 'e') {
            printf("\nTotal lattice energy: %f\n", lattice->findTotalEnergy());
        } else if (*input == 'm') {
            printf("\nLattice magnetization: %f\n", lattice->findMagnetism());
        } else if(*input == 'a') {
            int num = atoi(&input[1]);

            for (int i = 0; i <= num; i++) {
                lattice->updateLattice();
            }

            lattice->printLattice();
            printf("\nUpdated lattice %d times.", num);
            printf("\nTotal lattice energy: %f", lattice->findTotalEnergy());
            printf("\nLattice magnetization: %f\n", lattice->findMagnetism());
        } else if (*input == 'q') {
            cout << "\nExiting program...\n\n";
            return 0;
        } else if (*input == 'h') {
            help();
        } else {
            cout << "Invalid input! Please try again.\n\n";
        }

        cout << "Enter command (h - help): ";
    }
}

void help() {
    cout << "\nCommands: "
         << "\n\tp\t- Print lattice"
         << "\n\tu#\t- Update lattice # times"
         << "\n\te\t- Evaluate and return total lattice energy"
         << "\n\tm\t- Evaluate and return lattice magnetization"
         << "\n\ta#\t- Update lattice # times; "
                        "print lattice, energy, and magnetization"
         << "\n\tq\t- Quit program"
         << "\n\th\t- Access list of commands"
         << "\nEnter a command: ";
}
