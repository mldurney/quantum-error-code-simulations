#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Hamiltonian
{
public:
    Hamiltonian(vector< vector<int> > h, char s = '\0',
            int r = -1, int c = -1);
    ~Hamiltonian() {}

    vector< vector<int> > getHamiltonian() const { return hamiltonian; }
    vector<int> getIndices() const { return indices; }
    int getNumIndices() const { return numIndices; }
    vector< vector<int> > getLocalTerms() const { return localTerms; }

    char getShape() const { return shape; }
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    void printHamiltonian();
    void printIndices();
    void printLocalTerms();

    static vector< vector<int> > importHamiltonian(ifstream& file);

private:
    void generateIndices();
    void generateLocalTerms();

    vector< vector<int> > hamiltonian;
    vector<int> indices;
    int numIndices;
    vector< vector<int> > localTerms;

    char shape;
    int rows;
    int cols;
};

#endif /* HAMILTONIAN_H_ */
