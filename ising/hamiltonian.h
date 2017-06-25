#ifndef HAMILTONIAN_H_
#define HAMILTONIAN_H_

#include <vector>
#include <algorithm>

using namespace std;

class Hamiltonian
{
public:
    Hamiltonian(vector<vector<int>>, char);
    void printHamiltonian();
    void printIndices();
    void printLocalTerms();

protected:
    vector<vector<int>> hamiltonian;
    vector<int> indices;
    vector<vector<int>> localTerms;

private:
    void generateIndices();
    void generateLocalTerms();
};

#endif /* HAMILTONIAN_H_ */
