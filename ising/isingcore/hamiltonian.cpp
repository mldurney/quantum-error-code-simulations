#include "hamiltonian.h"

using namespace ising;

std::vector<std::vector<int>> ising::importHamiltonianVector(
    std::ifstream& file) {
    std::vector<std::vector<int>> hamiltonianVector;
    std::string line;
    int num;

    while (isalpha(file.peek())) {
        file.ignore(256, '\n');
    }

    while (getline(file, line)) {
        std::vector<int> interaction;
        std::istringstream lineStream(line);

        while (lineStream >> num) {
            interaction.push_back(num);

            if (lineStream.peek() == ',') {
                lineStream.ignore();
            }
        }

        hamiltonianVector.push_back(interaction);
    }

    return hamiltonianVector;
}

Hamiltonian::Hamiltonian(std::vector<std::vector<int>> h, char s, int r, int c)
    : hamiltonian(h), shape(s), rows(r), cols(c) {
    generateIndices();
    generateLocalTerms();
    generateIndInteractions();
}

void Hamiltonian::generateIndices() {
    std::vector<std::vector<int>>::iterator it1;
    std::vector<int>::iterator it2;
    std::vector<int>::iterator it3;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2) {
            it3 = find(indices.begin(), indices.end(), *it2);
            if (it3 == indices.end()) {
                indices.push_back(*it2);
            }
        }
    }

    sort(indices.begin(), indices.end());
    numIndices = indices.size();
}

void Hamiltonian::generateLocalTerms() {
    std::vector<std::vector<int>>::iterator it1;
    std::vector<int>::iterator it2;
    std::vector<int>::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2) {
        localTerms[*it2] = std::vector<int>();
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2) {
            for (it3 = it1->begin() + 1; it3 != it1->end(); ++it3) {
                if (it2 == it3) {
                    continue;
                }

                localTerms[*it2].push_back(*it3);
            }
        }
    }
}

void Hamiltonian::generateIndInteractions() {
    std::vector<std::vector<int>>::iterator it1;
    std::vector<int>::iterator it2;
    std::vector<int>::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2) {
        indInteractions[*it2] = std::vector<std::vector<int>>();
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2) {
            std::vector<int> interaction = {*(it1->begin())};

            for (it3 = it1->begin() + 1; it3 != it1->end(); ++it3) {
                if (it2 == it3) {
                    continue;
                }

                interaction.push_back(*it3);
            }

            indInteractions[*it2].push_back(interaction);
        }
    }
}

void Hamiltonian::printHamiltonian() {
    std::vector<std::vector<int>>::iterator it1;
    std::vector<int>::iterator it2;

    std::cout << "Printing Hamiltonian:" << std::endl;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2) {
            std::cout << *it2 << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printIndices() {
    std::vector<int>::iterator it;

    std::cout << "Printing indices:" << std::endl;

    for (it = indices.begin(); it != indices.end(); ++it) {
        std::cout << *it << " ";
    }

    std::cout << std::endl << std::endl;
}

void Hamiltonian::printLocalTerms() {
    std::map<int, std::vector<int>>::iterator it1;
    std::vector<int>::iterator it2;

    std::cout << "Printing local terms:" << std::endl;

    for (it1 = localTerms.begin(); it1 != localTerms.end(); ++it1) {
        std::cout << it1->first << ": ";

        for (it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2) {
            std::cout << *it2 << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printIndInteractions() {
    std::map<int, std::vector<std::vector<int>>>::iterator it1;
    std::vector<std::vector<int>>::iterator it2;
    std::vector<int>::iterator it3;

    std::cout << "Printing index interactions:" << std::endl;

    for (it1 = indInteractions.begin(); it1 != indInteractions.end(); ++it1) {
        std::cout << it1->first << ":\t";

        for (it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2) {
            for (it3 = it2->begin(); it3 != it2->end(); ++it3) {
                std::cout << *it3 << " ";
            }

            std::cout << "\t";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}