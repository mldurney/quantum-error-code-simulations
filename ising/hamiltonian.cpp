#include "hamiltonian.h"

Hamiltonian::Hamiltonian(vector<vector<int>> h, char s, int r, int c) :
    hamiltonian(h), shape(s), rows(r), cols(c)
{
    generateIndices();
    generateLocalTerms();
    generateIndInteractions();
}

void Hamiltonian::generateIndices()
{
    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;
    vector<int>::iterator it3;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2)
        {
            it3 = find(indices.begin(), indices.end(), *it2);
            if (it3 == indices.end())
            {
                indices.push_back(*it2);
            }
        }
    }

    sort(indices.begin(), indices.end());
    numIndices = indices.size();
}

void Hamiltonian::generateLocalTerms()
{
    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;
    vector<int>::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2)
    {
        localTerms[*it2] = vector<int>();
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2)
        {
            for (it3 = it1->begin() + 1; it3 != it1->end(); ++it3)
            {
                if (it2 == it3)
                {
                    continue;
                }

                localTerms[*it2].push_back(*it3);
            }
        }
    }
}

void Hamiltonian::generateIndInteractions()
{
    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;
    vector<int>::iterator it3;

    for (it2 = indices.begin(); it2 != indices.end(); ++it2)
    {
        indInteractions[*it2] = vector< vector<int> >();
    }

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2)
        {
            vector<int> interaction = {*(it1->begin())};

            for (it3 = it1->begin() + 1; it3 != it1->end(); ++it3)
            {
                if (it2 == it3)
                {
                    continue;
                }

                interaction.push_back(*it3);
            }

            indInteractions[*it2].push_back(interaction);
        }
    }
}

void Hamiltonian::printHamiltonian()
{
    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;

    cout << "Printing Hamiltonian:" << endl;

    for (it1 = hamiltonian.begin(); it1 != hamiltonian.end(); ++it1)
    {
        for (it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            cout << *it2 << " ";
        }

        cout << endl;
    }

    cout << endl;
}

void Hamiltonian::printIndices()
{
    vector<int>::iterator it;

    cout << "Printing indices:" << endl;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        cout << *it << " ";
    }

    cout << endl << endl;
}

void Hamiltonian::printLocalTerms()
{
    map< int, vector<int> >::iterator it1;
    vector<int>::iterator it2;

    cout << "Printing local terms:" << endl;

    for (it1 = localTerms.begin(); it1 != localTerms.end(); ++it1)
    {
        cout << it1->first << ": ";

        for (it2 = (it1->second).begin(); it2 != (it1->second).end(); ++it2)
        {
            cout << *it2 << " ";
        }

        cout << endl;
    }

    cout << endl;
}

void Hamiltonian::printIndInteractions()
{
    map< int, vector< vector<int> > >::iterator it1;
    vector< vector<int> >::iterator it2;
    vector<int>::iterator it3;

    cout << "Printing index interactions:" << endl;

    for (it1 = indInteractions.begin(); it1 != indInteractions.end(); ++it1)
    {
        cout << it1->first << ":\t";

        for (it2 = (it1->second).begin(); it2 != (it1->second).begin(); ++it2)
        {
            for (it3 = it2->begin(); it3 != it2->end(); ++it3)
            {
                cout << *it3 << " ";
            }

            cout << "\t";
        }

        cout << endl;
    }

    cout << endl;
}

vector< vector<int> > Hamiltonian::importHamiltonianVector(ifstream& file)
{
    vector< vector<int> > hamiltonianVector;
    string line;
    int num;

    while (isalpha(file.peek()))
    {
        file.ignore(256, '\n');
    }

    while (getline(file, line))
    {
        vector<int> interaction;
        istringstream lineStream(line);

        while (lineStream >> num)
        {
            interaction.push_back(num);

            if (lineStream.peek() == ',')
            {
                lineStream.ignore();
            }
        }

        hamiltonianVector.push_back(interaction);
    }

    return hamiltonianVector;
}
