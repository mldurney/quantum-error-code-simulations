#include <iostream>
#include "hamiltonian.h"

void Hamiltonian::printHamiltonian()
{
    vector<vector<int>>::iterator it1;
    vector<int>::iterator it2;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            cout << *it2 + " ";
        }

        cout << "\n";
    }
}

void Hamiltonian::printIndices()
{
    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        cout << *it + " ";
    }
}

void Hamiltonian::printLocalTerms()
{
    vector<vector<int>>::iterator it1;
    vector<int>::iterator it2;

    for (it1 = localTerms.begin(); it1 != localTerms.end(); ++it1)
    {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            cout << *it2 + " ";
        }

        cout << "\n";
    }
}

vector<vector<int>> importHamiltonianCSV(char* filename)
{
    hamiltonianVector = vector<vector<int>>;
    return hamiltonianVector;
}

int main()
{
    Hamiltonian h = Hamiltonian(importHamiltonianCSV)
    h.printHamiltonian();
    h.printIndices();
    h.printLocalTerms();
}
