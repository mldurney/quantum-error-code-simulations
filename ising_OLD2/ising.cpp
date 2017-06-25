#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "lattices.h"

using namespace std;

void help();

int main() {
    int side;
    double temp;
    char mode;
    char *input = new char[10];

    cout << "Enter number of rows: ";
    cin >> side;
    cout << "Enter temperature (K): ";
    cin >> temp;
    cout << "Enter update mode (a - all, r - random): ";
    cin >> mode;

    SquareLattice* lattice = new SquareLattice(side, temp, mode, 1);

    help();

    while (true) {
        cin >> input;
        int num;

        switch (*input)
        {
            case 'p':
                lattice->printLattice();
                break;

            case 'u':
                num = atoi(&input[1]);

                for (int i = 0; i <= num; i++)
                {
                    lattice->updateLattice();
                }

                printf("\nUpdated lattice %d times.\n", num);

                break;

            case 'e':
                printf("\nTotal lattice energy: %f\n",
                        lattice->findTotalEnergy());
                break;

            case 'm':
                printf("\nLattice magnetization: %f\n",
                        lattice->findMagnetism());
                break;

            case 'a':
                num = atoi(&input[1]);

                for (int i = 0; i <= num; i++) {
                    lattice->updateLattice();
                }

                lattice->printLattice();

                printf("\nUpdated lattice %d times.", num);
                printf("\nTotal lattice energy: %f", lattice->findTotalEnergy());
                printf("\nLattice magnetization: %f\n", lattice->findMagnetism());

                break;

            case 's':
                printf("\nSwitching mode to update %s.\n",
                        (mode == ALL) ? "random cells" : "all cells");
                lattice->switchMode(mode);
                break;

            case 'q':
                cout << "\nExiting program...\n\n";
                return 0;
                break;

            case 'h':
                help();
                break;

            default:
                cout << "Invalid input! Please try again.\n\n";
                break;
        }

        cout << "Enter command (h - help): ";
    }
}

void help()
{
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
