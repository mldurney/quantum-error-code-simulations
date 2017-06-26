#include "lattices.h"


/////////////
// Lattice //
/////////////

Lattice::Lattice(Hamiltonian h, double t, char m = 'r') :
        hamiltonian(h), temp(t), mode(m)
{
    hFunction = hamiltonian.getHamiltonian();
    indices = hamiltonian.getIndices();
    numIndices = hamiltonian.getNumIndices();
    localTerms = hamiltonian.getLocalTerms();
    shape = hamiltonian.getShape();

    initSpins();
}

void Lattice::initSpins()
{
    srand(time(NULL));

    vector<int>::iterator it;

    spins.resize(*(indices.end()) + 1);

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        spins.insert(it, (rand() % 2 == 1) ? 1 : -1);
    }
}

void Lattice::updateLattice() {
    switch (mode)
    {
        case ALL: updateAll(); break;
        case RANDOM: updateRandom(); break;
    }
}

void Lattice::updateAll()
{
    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        if (findProbability(*it) > (double) rand() / (double) RAND_MAX)
        {
            spins[*it] *= -1;
        }
    }
}

void Lattice::updateRandom()
{
    int index;

    for (int i = 0; i < numIndices; ++i)
    {
        index = indices[rand() % numIndices];

        if (findProbability(index)
                > (double) rand() / (double) RAND_MAX)
        {
            spins[index] *= -1;
        }
    }
}

void Lattice::switchMode(char mode)
{
    switch (mode)
    {
        case ALL:
            mode = ALL;
            break;
        case RANDOM:
            mode = RANDOM;
            break;
        default:
            cout << "INVALID MODE. Exiting...\n\n";
            exit(EXIT_FAILURE);
    }
}

double Lattice::findProbability(int index)
{
    double initEnergy = findIndexEnergy(index);
    spins[index] *= -1;

    double finalEnergy = findIndexEnergy(index);
    spins[index] *= -1;

    if (finalEnergy > initEnergy)
    {
        return pow(E, (-1 / temp) * (finalEnergy - initEnergy));
    }
    else
    {
        return 1;
    }
}

double Lattice::findTotalEnergy()
{
    double energy = 0;
    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        energy += findIndexEnergy(*it);
    }

    return energy;
}

double Lattice::findIndexEnergy(int index)
{
    double energy = -1;

    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;

    for (it1 = hFunction.begin(); it1 != hFunction.end(); ++it1)
    {
        for (it2 = it1->begin() + 1; it2 != it1->end(); ++it2)
        {
            if (index == *it2)
            {
                for (it2 = it1->begin(); it2 != it1->end(); ++it2)
                {
                    energy *= *it2;
                }

                break;
            }
        }
    }

    return energy;
}

double Lattice::findMagnetism()
{
    double magnetism = 0;
    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        magnetism += spins[*it];
    }

    return magnetism / numIndices;
}

void Lattice::shapeError() const
{
    cout << "Wrong lattice type. Aborting...\n\n";
    exit(EXIT_FAILURE);
}

void Lattice::printLattice(int cols = -1)
{
    if (cols == -1)
    {
        cols = (int) sqrt(numIndices);
    }

    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        for (int col = 0; col < cols && it != indices.end(); ++col, ++it)
        {
            (spins[*it] == 1) ? cout << "+ " : cout << "- ";
        }

        cout << "\n";
    }
}


////////////////////////
// RectangularLattice //
////////////////////////

RectangularLattice::RectangularLattice(Hamiltonian h, double t, char m = 'r',
        int r = -1, int c = -1) : Lattice(h, t, m), rows(r), cols(c)
{
    checkShape();

    if (rows == -1 || cols == -1)
    {
        guessRowsCols();
    }
}

void RectangularLattice::checkShape() const
{
    if (shape != RECTANGLE || shape != SQUARE)
    {
        cout << "Invalid shape parameter -- not a rectangle ('r')!\n";
        shapeError();
    }
}

void RectangularLattice::guessRowsCols()
{
    for (int guess = (int) sqrt(numIndices); guess > 0; --guess)
    {
        if (numIndices % guess == 0)
        {
            setRows(guess);
            setCols(numIndices / guess);
            return;
        }
    }
}


///////////////////
// SquareLattice //
///////////////////


SquareLattice::SquareLattice(Hamiltonian h, double t, char m = 'r',
        int s = -1) : RectangularLattice(h, t, m, s, s)
{
    checkShape();

    if (side == -1)
    {
        guessSide();
    }
}

void SquareLattice::checkShape() const
{
    if (shape != SQUARE)
    {
        cout << "Invalid shape parameter -- not a square ('s')!\n";
        shapeError();
    }
}

void SquareLattice::guessSide()
{
    int guess = (int) sqrt(numIndices);

    if (numIndices / guess != guess)
    {
        cout << "Lattice is not a square (invalid number of indices)!";
        shapeError();
    }

    setSide(guess);
    setRows(guess);
    setCols(guess);
}


///////////////////////
// TriangularLattice //
///////////////////////


TriangularLattice::TriangularLattice(Hamiltonian h, double t, char m = 'r',
        int r = -1, int c = -1) : Lattice(h, t, m), rows(r), cols(c)
{
    checkShape();

    if (rows == -1 || cols == -1)
    {
        guessRowsCols();
    }
}

void TriangularLattice::checkShape() const
{
    if (shape != TRIANGLE)
    {
        cout << "Invalid shape parameter -- not a triangle ('t')!\n";
        shapeError();
    }
}

void TriangularLattice::guessRowsCols()
{
    for (int guess = (int) sqrt(numIndices); guess > 0; --guess)
    {
        if (numIndices % guess == 0)
        {
            setRows(guess);
            setCols(numIndices / guess);
            return;
        }
    }
}


int main()
{
    return 0;
}
