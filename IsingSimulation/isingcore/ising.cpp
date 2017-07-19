#include "ising.h"

using namespace ising;

int main(int argc, char *argv[]) {
    std::string filename = receiveHamiltonianFile(argc, argv);
    std::ifstream file(filename);

    if (!file) {
        std::cout << "Invalid file name. " << filename
                  << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    char shape;
    Hamiltonian hamiltonian = readHamiltonian(file, shape);

    file.close();

    double temp;
    char mode;

    std::cout << "Enter temperature (K): ";
    std::cin >> temp;
    std::cout << "Enter update mode (a - all, p - pseudorandom, r - random): ";
    std::cin >> mode;

    Lattice *lattice = chooseLattice(shape, hamiltonian, temp, mode);
    runLattice(lattice);
    delete lattice;

    return 0;
}

std::string ising::receiveHamiltonianFile(int argc, char *argv[]) {
    std::string filename;

    if (argc == 1) {
        std::cout << "Enter Hamiltonian input file: ";
        std::cin >> filename;
    } else if (argc == 2) {
        filename = argv[1];
    } else {
        std::cout << "Usage: " << argv[0] << " name_of_hamiltonian_file\n\n";
        exit(EXIT_FAILURE);
    }

    return filename;
}

void ising::runLattice(Lattice *lattice) {
    char *input = new char[10];

    runLatticeHelp();

    while (true) {
        std::cin >> input;

        int num;
        switch (*input) {
            case 'p':
                lattice->printLattice();
                break;

            case 'u':
                num = atoi(&input[1]);

                for (int i = 0; i < num; i++) {
                    lattice->updateLattice();
                }

                std::cout << "\nUpdated lattice " << num << " times.\n";
                break;

            case 'e':
                std::cout << "\nTotal lattice energy: ";
                std::cout << lattice->getTotalEnergy() << std::endl;
                break;

            case 'm':
                std::cout << "\nLattice magnetization: ";
                std::cout << lattice->getMagnetism() << std::endl;
                ;
                break;

            case 'a':
                num = atoi(&input[1]);

                for (int i = 0; i < num; i++) {
                    lattice->updateLattice();
                }

                lattice->printLattice();
                std::cout << "\nUpdated lattice " << num << " times.";
                std::cout << "\nTotal lattice energy: ";
                std::cout << lattice->getTotalEnergy();
                std::cout << "\nLattice magnetization: ";
                std::cout << lattice->getMagnetism() << std::endl;

                break;

            case 's':
                std::cout << "\nSwitching mode to " << input[1] << std::endl;
                lattice->switchMode(input[1]);
                break;

            case 'q':
                std::cout << "\nExiting program...\n\n";
                return;

            case 'h':
                runLatticeHelp();
                break;

            default:
                std::cout << "Invalid input! Please try again.\n\n";
        }

        std::cout << "\nEnter command (h - help): ";
    }

    delete[] input;
}

void ising::runLatticeHelp() {
    std::cout
        << "\nCommands: "
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
