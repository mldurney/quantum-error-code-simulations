#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "common.h"

namespace ising {
ivector2 importHamiltonianVector(std::ifstream &file);

class Hamiltonian {
   public:
    Hamiltonian(ivector2 h, char s = '\0', int r = -1, int c = -1);
    ~Hamiltonian() {}

    const ivector2 getHamiltonian() const { return hamiltonian; }
    uint getNumIndices() const { return numIndices; }
    const ivector getIndices() const { return indices; }
    const i2arraymap getLocations() const { return locations; }
    const ivectormap getLocalTerms() const { return localTerms; }
    const std::map<int, ivector2> getIndInteractions() const {
        return indInteractions;
    }

    char getShape() const { return shape; }
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    void printHamiltonian() const;
    void printIndices() const;
    void printLocations() const;
    void printLocalTerms() const;
    void printIndInteractions() const;

   private:
    void generateIndices();
    void generateLocations();
    void generateLocalTerms();
    void generateIndInteractions();

    ivector2 hamiltonian;
    uint numIndices;
    ivector indices;
    i2arraymap locations;
    ivectormap localTerms;
    std::map<int, ivector2> indInteractions;

    char shape;
    int rows;
    int cols;
};
}

#endif /* HAMILTONIAN_H_ */
