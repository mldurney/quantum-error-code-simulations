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
 *      Explore possibilities of distance
 *      Create list of boards to see progression over time
 *      Enable output of results
 *      Error checking
 *
 *      Include option to compare changes
 *
 *      Desired behavior:
 *          Plot of average value of absolute magnetization vs temp

- Obtain a plot of the average (absolute value) of the magnetization as a function of temperature
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "lattice_near.h"

using namespace std;

void help();

int main() {
    int cols, rows, temp;
    char mode;
    char *input = new char[10];

    cout << "Enter number of rows: ";
    cin >> rows;
    cout << "Enter number of cols: ";
    cin >> cols;
    cout << "Enter temperature (K): ";
    cin >> temp;
    cout << "Enter update mode (a - all, r - random): ";
    cin >> mode;

    Lattice* lattice = new Lattice(rows, cols, temp, mode);

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
        } else if (*input == 's') {
            printf("\nSwitching mode to update %s.\n",
                    (mode == ALL) ? "random cells" : "all cells");
            lattice->switchMode(mode);
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
         << "\n\ts\t- Switch modes (update all cells <-> update random cells)"
         << "\n\tq\t- Quit program"
         << "\n\th\t- Access list of commands"
         << "\nEnter a command: ";
}
