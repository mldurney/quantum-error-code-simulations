#include "ising.h"

int main(int argc, char* argv[])
{
    string filename = receiveFilename(argc, argv);
    ifstream file(filename);

    if (!file)
    {
        cout<< "Invalid file name. " << filename << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    char shape;
    Hamiltonian hamiltonian = readHamiltonian(file, shape);

    file.close();

    double temp;
    char mode;

    cout << "Enter temperature (K): ";
    cin >> temp;
    cout << "Enter update mode (a - all, p - pseudorandom, r - random): ";
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

    runLattice(lattice);

    delete lattice;
}

string receiveFilename(int argc, char* argv[])
{
    string filename;

    if (argc == 1)
    {
        cout << "Enter Hamiltonian input file: ";
        cin >> filename;
    }
    else if (argc == 2)
    {
        filename = argv[1];
    }
    else
    {
        cout << "Usage: " << argv[0] << " name_of_hamiltonian_file\n\n";
        exit(EXIT_FAILURE);
    }

    return filename;
}

void runLattice(Lattice* lattice)
{
    char *input = new char[10];

    runLatticeHelp();

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

                cout << "\nUpdated lattice " << num << " times.\n";
                break;

            case 'e':
                cout << "\nTotal lattice energy: ";
                cout << lattice->getTotalEnergy() << endl;
                break;

            case 'm':
                cout << "\nLattice magnetization: ";
                cout << lattice->getMagnetism() << endl;;
                break;

            case 'a':
                num = atoi(&input[1]);

                for (int i = 0; i < num; i++) {
                    lattice->updateLattice();
                }

                lattice->printLattice();
                cout << "\nUpdated lattice " << num << " times.";
                cout << "\nTotal lattice energy: ";
                cout << lattice->getTotalEnergy();
                cout << "\nLattice magnetization: ";
                cout << lattice->getMagnetism() << endl;

                break;

            case 's':
                cout<< "\nSwitching mode to " << input[1] << endl;
                lattice->switchMode(input[1]);
                break;

            case 'q':
                cout << "\nExiting program...\n\n";
                return;

            case 'h':
                runLatticeHelp();
                break;

            default:
                cout << "Invalid input! Please try again.\n\n";
        }

        cout << "\nEnter command (h - help): ";
    }
}

void runLatticeHelp() {
    cout << "\nCommands: "
         << "\n\tp\t- Print lattice"
         << "\n\tu#\t- Update lattice # times"
         << "\n\te\t- Evaluate and return total lattice energy"
         << "\n\tm\t- Evaluate and return lattice magnetization"
         << "\n\ta#\t- Update lattice # times; "
                        "print lattice, energy, and magnetization"
         << "\n\ts#\t- Switch modes (a - all, p - pseudorandom, r - random)"
         << "\n\tq\t- Quit program"
         << "\n\th\t- Access list of commands"
         << "\n\nEnter a command: ";
}
