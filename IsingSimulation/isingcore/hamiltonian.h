#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <array>

typedef std::vector<int> ivector;
typedef std::vector<ivector> ivector2;
typedef std::map<int, ivector> ivectormap;
typedef std::array<int, 2> i2array;
typedef std::map<int, i2array> i2arraymap;

namespace ising {
ivector2 importHamiltonianVector(std::ifstream &file);

class Hamiltonian {
   public:
    Hamiltonian(ivector2 h, char s = '\0', int r = -1, int c = -1);
    ~Hamiltonian() {}

    const ivector2 getHamiltonian() const { return hamiltonian; }
    const ivector getIndices() const { return indices; }
	const i2arraymap getLocations() const { return locations; }
    const int getNumIndices() const { return numIndices; }
    const ivectormap getLocalTerms() const { return localTerms; }
    const std::map<int, ivector2> getIndInteractions() const {
        return indInteractions;
    }

    const char getShape() const { return shape; }
    const int getRows() const { return rows; }
    const int getCols() const { return cols; }
    const bool getIsFast() const { return isFast; }

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
    void findIsFast();

    ivector2 hamiltonian;
    ivector indices;
	i2arraymap locations;
    int numIndices;
    ivectormap localTerms;
    std::map<int, ivector2> indInteractions;

    char shape;
    int rows;
    int cols;
    bool isFast;
};
}

#endif /* HAMILTONIAN_H_ */
