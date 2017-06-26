#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <vector>
#include <algorithm>

using namespace std;

class Hamiltonian
{
public:
    Hamiltonian(vector< vector<int> > h, char s);
    vector< vector<int> > getHamiltonian() const { return hamiltonian; }
    vector<int> getIndices() const { return indices; }
    int getNumIndices() const { return numIndices; }
    vector< vector<int> > getLocalTerms() const { return localTerms; }
    char getShape() const { return shape; }
    void printHamiltonian();
    void printIndices();
    void printLocalTerms();
    ~Hamiltonian() {}

private:
    void generateIndices();
    void generateLocalTerms();

    vector< vector<int> > hamiltonian;
    vector<int> indices;
    int numIndices;
    vector< vector<int> > localTerms;
    char shape;
};

#endif /* HAMILTONIAN_H_ */
