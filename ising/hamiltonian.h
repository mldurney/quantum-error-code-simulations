#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <map>
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
    map< int, vector<int> > getLocalTerms() const { return localTerms; }
    map< int, vector< vector<int> > > getIndInteractions() const
            { return indInteractions; }

    char getShape() const { return shape; }
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    void printHamiltonian();
    void printIndices();
    void printLocalTerms();
    void printIndInteractions();

    static vector< vector<int> > importHamiltonianVector(ifstream& file);

private:
    void generateIndices();
    void generateLocalTerms();
    void generateIndInteractions();

    vector< vector<int> > hamiltonian;
    vector<int> indices;
    int numIndices;
    map< int, vector<int> > localTerms;
    map< int, vector< vector<int> > > indInteractions;

    char shape;
    int rows;
    int cols;
};

#endif /* HAMILTONIAN_H_ */
