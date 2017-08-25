#include "hamiltonian.h"

using namespace ising;

ivector2 ising::importHamiltonianVector(std::ifstream &file) {
    ivector2 hamiltonianVector;
    std::string line;
    int num;

    while (isalpha(file.peek())) {
        file.ignore(256, '\n');
    }

    while (getline(file, line)) {
        ivector interaction;
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

Hamiltonian::Hamiltonian(ivector2 h, char s, int r, int c)
    : hamiltonian(h), shape(s), rows(r), cols(c) {
    generateIndices();
    generateLocations();
    generateLocalTerms();
    generateIndInteractions();
}

void Hamiltonian::generateIndices() {
    ivector::iterator it;

    for (auto &interaction : hamiltonian) {
        for (it = interaction.begin() + 1; it != interaction.end(); ++it) {
            if (std::find(indices.begin(), indices.end(), *it) ==
                indices.end()) {
                indices.push_back(*it);
            }
        }
    }

    std::sort(indices.begin(), indices.end());
    numIndices = (int)indices.size();
}

void Hamiltonian::generateLocations() {
    if (rows == -1 || cols == -1) {
        return;
    }

    for (auto &i : indices) {
        int row = i / cols;
        int col = i % cols;
        locations[i] = {row, col};

        if (row >= rows || col >= cols) {
            locations.clear();
            return;
        }
    }
}

void Hamiltonian::generateLocalTerms() {
    ivector::iterator it1;
    ivector::iterator it2;

    for (auto &i : indices) {
        localTerms[i] = ivector();
    }

    for (auto &h : hamiltonian) {
        for (it1 = h.begin() + 1; it1 != h.end(); ++it1) {
            for (it2 = h.begin() + 1; it2 != h.end(); ++it2) {
                if (it1 == it2) {
                    continue;
                }

                localTerms[*it1].push_back(*it2);
            }
        }
    }
}

void Hamiltonian::generateIndInteractions() {
    ivector::iterator it1, it2;

    for (auto &i : indices) {
        indInteractions[i] = ivector2();
    }

    for (auto &h : hamiltonian) {
        for (it1 = h.begin() + 1; it1 != h.end(); ++it1) {
            ivector interaction = {*(h.begin())};

            for (it2 = h.begin() + 1; it2 != h.end(); ++it2) {
                if (it1 == it2) {
                    continue;
                }

                interaction.push_back(*it2);
            }

            indInteractions[*it1].push_back(interaction);
        }
    }
}

void Hamiltonian::printHamiltonian() const {
    std::cout << "Printing Hamiltonian:" << std::endl;

    for (auto &interaction : hamiltonian) {
        for (auto &i : interaction) {
            std::cout << i << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printIndices() const {
    std::cout << "Printing indices:" << std::endl;

    for (auto &i : indices) {
        std::cout << i << " ";
    }

    std::cout << std::endl << std::endl;
}

void Hamiltonian::printLocations() const {
    std::cout << "Printing locations:" << std::endl;

    for (auto &loc : locations) {
        std::cout << loc.first << ":\t" << loc.second[0] << " " << loc.second[1]
                  << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printLocalTerms() const {
    std::cout << "Printing local terms:" << std::endl;

    for (auto &terms : localTerms) {
        std::cout << terms.first << ":\t";

        for (auto &i : terms.second) {
            std::cout << i << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printIndInteractions() const {
    std::map<int, ivector2>::const_iterator it1;
    ivector2::const_iterator it2;
    ivector::const_iterator it3;

    std::cout << "Printing index interactions:" << std::endl;

    for (auto &interaction : indInteractions) {
        std::cout << interaction.first << ":\t";

        for (auto &ind : interaction.second) {
            for (auto &i : ind) {
                std::cout << i << " ";
            }

            std::cout << "\t";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}
