#include "hamiltonian.h"

Hamiltonian::Hamiltonian(vector<vector<int>> h) : hamiltonian(h)
{
    generateIndices();
    generateLocalTerms();
}

void Hamiltonian::generateIndices()
{
    vector<vector<int>>::iterator it1;
    vector<int>::iterator it2;
    vector<int>::iterator it3;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            it3 = find(indices.begin(), indices.end(), *it2);
            if (it3 == indices.end())
            {
                indices.push_back(*it2);
            }
        }
    }

    sort(indices.begin(), indices.end());
}

void Hamiltonian::generateLocalTerms()
{
    vector<vector<int>>::iterator it1;
    vector<int>::iterator it2;
    vector<int>::iterator it3;

    for (it3 = indices.begin(); it3 != indices.end(); ++it3)
    {
        localTerms.push_back(vector<int>());
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = ++(it1->begin()); it2 != it1->end(); ++it2)
        {
            for (it3 = it2 + 1; it3 != it1->end(); ++it3)
            {
                localTerms[*it2].push_back(*it3);
            }
        }
    }
}
