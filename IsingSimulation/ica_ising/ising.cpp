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

    double temperature;
    char mode;
    uint preupdates;

    if (argc == 4 || argc == 5) {
        temperature = atof(argv[2]);
        mode = (argv[3])[0];

        if (argc == 5) {
            preupdates = 100;
        }

    } else {
        std::cout << "Enter temperature (K): ";
        std::cin >> temperature;
        std::cout
            << "Enter update mode (a - all, p - pseudorandom, r - random): ";
        std::cin >> mode;
    }

    Lattice *lattice =
        chooseLattice(shape, hamiltonian, temperature, 0, 1, mode);
    Replica replica = lattice->getReplicaCopy(0);
    runReplica(replica, preupdates);
}

std::string ising::receiveHamiltonianFile(int argc, char *argv[]) {
    std::string filename;

    if (argc == 1) {
        std::cout << "Enter Hamiltonian input file: ";
        std::cin >> filename;
    } else if (argc == 2 || argc == 4 || argc == 5) {
        filename = argv[1];
    } else {
        std::cout << "Usage: " << argv[0]
                  << " name_of_hamiltonian_file [temperature mode]\n\n";
        exit(EXIT_FAILURE);
    }

    return filename;
}

void ising::runReplica(Replica replica, uint preupdates) {
    char *input = new char[10];

    if (preupdates == 0) {
        runReplicaHelp();
    } else {
        for (uint i = 0; i < preupdates; i++) {
            replica.update();
        }

        std::cout << std::endl;

        replica.print();
        std::cout << "\nUpdated lattice " << preupdates << " times.";
        std::cout << "\nTotal lattice energy: ";
        std::cout << replica.getTotalEnergy();
        std::cout << "\nLattice magnetization: ";
        std::cout << replica.getMagnetization() << std::endl;

        std::cout << "\nEnter command (h - help): ";
    }

    while (true) {
        std::cin >> input;

        uint num;
        switch (*input) {
            case 'p':
                replica.print();
                break;

            case 'u':
                num = atoi(&input[1]);

                for (uint i = 0; i < num; i++) {
                    replica.update();
                }

                std::cout << "\nUpdated lattice " << num << " times.\n";
                break;

            case 'e':
                std::cout << "\nTotal lattice energy: ";
                std::cout << replica.getTotalEnergy() << std::endl;
                break;

            case 'm':
                std::cout << "\nLattice magnetization: ";
                std::cout << replica.getMagnetization() << std::endl;
                ;
                break;

            case 'a':
                num = atoi(&input[1]);

                for (uint i = 0; i < num; i++) {
                    replica.update();
                }

                replica.print();
                std::cout << "\nUpdated lattice " << num << " times.";
                std::cout << "\nTotal lattice energy: ";
                std::cout << replica.getTotalEnergy();
                std::cout << "\nLattice magnetization: ";
                std::cout << replica.getMagnetization() << std::endl;

                break;

            case 'q':
                std::cout << "\nExiting program...\n\n";
                return;

            case 'h':
                runReplicaHelp();
                break;

            default:
                std::cout << "Invalid input! Please try again.\n\n";
        }

        std::cout << "\nEnter command (h - help): ";
    }

    delete[] input;
}

void ising::runReplicaHelp() {
    std::cout << "\nCommands: "
              << "\n\tp\t- Print lattice"
              << "\n\tu#\t- Update lattice # times"
              << "\n\te\t- Evaluate and return total lattice energy"
              << "\n\tm\t- Evaluate and return lattice magnetization"
              << "\n\ta#\t- Update lattice # times; "
                 "print lattice, energy, and magnetization"
              << "\n\tq\t- Quit program"
              << "\n\th\t- Access list of commands"
              << "\n\nEnter a command: ";
}
