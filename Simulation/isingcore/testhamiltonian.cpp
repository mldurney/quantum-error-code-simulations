#include <iostream>
#include "hamiltonian.h"

using namespace ising;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s name_of_hamiltonian_file\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    std::ifstream file(argv[1]);

    if (!file) {
        printf("Invalid file name. %s does not exist!\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    char shape = '\0';
    int rows = -1;
    int cols = -1;
    char c = (char)file.peek();

    if (isalpha(c)) {
        std::string line;
        getline(file, line);
        std::istringstream lineStream(line);

        lineStream >> shape;
        lineStream.ignore();
        lineStream >> rows;
        lineStream.ignore();
        lineStream >> cols;
    }

    Hamiltonian h =
        Hamiltonian(importHamiltonianVector(file), shape, rows, cols);

    file.close();

    std::cout << std::endl;

    h.printHamiltonian();
    h.printIndices();
    h.printLocations();
    h.printLocalTerms();
    h.printIndInteractions();
}
