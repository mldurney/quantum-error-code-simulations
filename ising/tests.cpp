#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "hamiltonian.h"

vector< vector<int> > importHamiltonianCSV(ifstream&);

void Hamiltonian::printHamiltonian()
{
    vector< vector<int> >::iterator it1;
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
    vector< vector<int> >::iterator it1;
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

vector< vector<int> > importHamiltonianCSV(ifstream& file)
{
    vector< vector<int> > hamiltonianVector;
    string line;
    int num;

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

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s name_of_hamiltonian_file\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char* filename = argv[1];
    ifstream file(filename);

    char shape = '\0';
    char c = file.peek();

    if (isalpha(c))
    {
        shape = c;
        file.ignore(256, '\n');
    }

    Hamiltonian h = Hamiltonian(importHamiltonianCSV(file), shape);
    h.printHamiltonian();
    h.printIndices();
    h.printLocalTerms();
}
