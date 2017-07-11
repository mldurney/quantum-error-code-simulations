#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace ising {
std::vector<std::vector<int>> importHamiltonianVector(std::ifstream &file);

class Hamiltonian {
   public:
    Hamiltonian(std::vector<std::vector<int>> h, char s = '\0', int r = -1,
                int c = -1);
    ~Hamiltonian() {}

    std::vector<std::vector<int>> getHamiltonian() const { return hamiltonian; }
    std::vector<int> getIndices() const { return indices; }
    int getNumIndices() const { return numIndices; }
    std::map<int, std::vector<int>> getLocalTerms() const { return localTerms; }
    std::map<int, std::vector<std::vector<int>>> getIndInteractions() const {
        return indInteractions;
    }

    char getShape() const { return shape; }
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    void printHamiltonian();
    void printIndices();
    void printLocalTerms();
    void printIndInteractions();

   private:
    void generateIndices();
    void generateLocalTerms();
    void generateIndInteractions();

    std::vector<std::vector<int>> hamiltonian;
    std::vector<int> indices;
    int numIndices;
    std::map<int, std::vector<int>> localTerms;
    std::map<int, std::vector<std::vector<int>>> indInteractions;

    char shape;
    int rows;
    int cols;
};
}

#endif /* HAMILTONIAN_H_ */
