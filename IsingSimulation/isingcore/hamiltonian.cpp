#include "hamiltonian.h"

using namespace ising;

ivector2 ising::importHamiltonianVector(
    std::ifstream& file) {
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
    generateLocalTerms();
    generateIndInteractions();
    findIsFast();
}

void Hamiltonian::generateIndices() {
    ivector2::iterator it1;
    ivector::iterator it2;
    ivector::iterator it3;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2) {
            it3 = find(indices.begin(), indices.end(), *it2);
            if (it3 == indices.end()) {
                indices.push_back(*it2);
            }
        }
    }

    sort(indices.begin(), indices.end());
    numIndices = (int)indices.size();
}

void Hamiltonian::generateLocalTerms() {
    ivector2::iterator it1;
    ivector::iterator it2;
    ivector::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2) {
        localTerms[*it2] = ivector();
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
    ivector2::iterator it1;
    ivector::iterator it2;
    ivector::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2) {
        indInteractions[*it2] = ivector2();
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2) {
            ivector interaction = {*(it1->begin())};

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

void Hamiltonian::findIsFast() {
    ivector2::iterator it = hamiltonian.begin();
    char c = (char)(*it)[0];

    for (++it; it != hamiltonian.end(); ++it) {
        if ((*it)[0] != c) {
            isFast = false;
            break;
        }
    }

    for (int i = 0; i < numIndices; ++i) {
        if (indices[i] != i) {
            isFast = false;
            break;
        }
    }
}

void Hamiltonian::printHamiltonian() const {
    ivector2::const_iterator it1;
    ivector::const_iterator it2;

    std::cout << "Printing Hamiltonian:" << std::endl;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1) {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2) {
            std::cout << *it2 << " ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void Hamiltonian::printIndices() const {
    ivector::const_iterator it;

    std::cout << "Printing indices:" << std::endl;

    for (it = indices.begin(); it != indices.end(); ++it) {
        std::cout << *it << " ";
    }

    std::cout << std::endl << std::endl;
}

void Hamiltonian::printLocalTerms() const {
    ivectormap::const_iterator it1;
    ivector::const_iterator it2;

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

void Hamiltonian::printIndInteractions() const {
    std::map<int, ivector2>::const_iterator it1;
    ivector2::const_iterator it2;
    ivector::const_iterator it3;

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
