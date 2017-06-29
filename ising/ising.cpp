#include <stdio.h>
#include "lattices.h"

using namespace std;

void help();
void runSimulation(Lattice*);

int main(int argc, char* argv[]) {
    if (argc != 2)
    {
        printf("Usage: %s name_of_hamiltonian_file\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char* filename = argv[1];
    ifstream file(filename);

    if (!file)
    {
        printf("Invalid file name. %s does not exist!\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }


    char shape = '\0';
    int rows = -1;
    int cols = -1;

    if (isalpha(file.peek()))
    {
        string line;
        string field;
        getline(file, line);

        stringstream stream(line);
        stringstream fields[3];

        fields[0] << "\0";
        fields[1] << "-1";
        fields[2] << "-1";

        for (int i = 0; i < 3 && getline(stream, field, ','); ++i)
        {
            fields[i].str(field);
        }

        fields[0] >> shape;
        fields[1] >> rows;
        fields[2] >> cols;
    }


    Hamiltonian hamiltonian =
            Hamiltonian(Hamiltonian::importHamiltonian(file),
                    shape, rows, cols);

    double temp;
    char mode;

    cout << "Enter temperature (K): ";
    cin >> temp;
    cout << "Enter update mode (a - all, r - random): ";
    cin >> mode;

    Lattice* lattice;

    switch (shape)
    {
        case RECTANGLE:
            lattice = new RectangularLattice(hamiltonian, temp, mode);
            break;
        case SQUARE:
            lattice = new SquareLattice(hamiltonian, temp, mode);
            break;
        case TRIANGLE:
            lattice = new TriangularLattice(hamiltonian, temp, mode);
            break;
        default:
            lattice = new Lattice(hamiltonian, temp, mode);
    }

    runSimulation(lattice);

    delete lattice;
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
         << "\n\nEnter a command: ";
}

void runSimulation(Lattice* lattice)
{
    char *input = new char[10];

    help();

    while (true)
    {
        cin >> input;

        int num;
        switch (*input)
        {
            case 'p':
                lattice->printLattice();
                break;

            case 'u':
                num = atoi(&input[1]);

                for (int i = 0; i < num; i++) {
                    lattice->updateLattice();
                }

                printf("\nUpdated lattice %d times.\n", num);
                break;

            case 'e':
                printf("\nTotal lattice energy: %d\n",
                        lattice->getTotalEnergy());
                break;

            case 'm':
                printf("\nLattice magnetization: %f\n",
                        lattice->getMagnetism());
                break;

            case 'a':
                num = atoi(&input[1]);

                for (int i = 0; i < num; i++) {
                    lattice->updateLattice();
                }

                lattice->printLattice();
                printf("\nUpdated lattice %d times.", num);
                printf("\nTotal lattice energy: %d",
                        lattice->getTotalEnergy());
                printf("\nLattice magnetization: %f\n",
                        lattice->getMagnetism());

                break;

            case 's':
                printf("\nSwitching mode to update %s.\n",
                        (lattice->getMode() == ALL) ?
                        "random cells" : "all cells");
                lattice->switchMode();
                break;

            case 'q':
                cout << "\nExiting program...\n\n";
                return;

            case 'h':
                help();
                break;

            default:
                cout << "Invalid input! Please try again.\n\n";
        }

        cout << "\nEnter command (h - help): ";
    }
}
