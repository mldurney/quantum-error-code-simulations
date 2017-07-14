#include <iostream>
#include "hamiltonian.h"

using namespace ising;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s name_of_hamiltonian_file\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *filename = argv[1];
    std::ifstream file(filename);

    if (!file) {
        printf("Invalid file name. %s does not exist!\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    char shape = '\0';
    char c = file.peek();

    if (isalpha(c)) {
        shape = c;
        file.ignore(256, '\n');
    }

    Hamiltonian h = Hamiltonian(importHamiltonianVector(file), shape);

    file.close();

    std::cout << std::endl;

    h.printHamiltonian();
    h.printIndices();
    h.printLocalTerms();
    h.printIndInteractions();
}
