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

    const std::vector<std::vector<int>> getHamiltonian() const { return hamiltonian; }
    const std::vector<int> getIndices() const { return indices; }
    const int getNumIndices() const { return numIndices; }
    const std::map<int, std::vector<int>> getLocalTerms() const { return localTerms; }
    const std::map<int, std::vector<std::vector<int>>> getIndInteractions() const {
        return indInteractions;
    }

    const char getShape() const { return shape; }
    const int getRows() const { return rows; }
    const int getCols() const { return cols; }
	const bool getIsFast() const { return isFast; }

    void printHamiltonian() const;
    void printIndices() const;
    void printLocalTerms() const;
    void printIndInteractions() const;

   private:
    void generateIndices();
    void generateLocalTerms();
    void generateIndInteractions();
	void findIsFast();

    std::vector<std::vector<int>> hamiltonian;
    std::vector<int> indices;
    int numIndices;
    std::map<int, std::vector<int>> localTerms;
    std::map<int, std::vector<std::vector<int>>> indInteractions;

    char shape;
    int rows;
    int cols;
	bool isFast;
};
}

#endif /* HAMILTONIAN_H_ */
