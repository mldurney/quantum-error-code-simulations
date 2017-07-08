#include "lattices.h"


/////////////
// Lattice //
/////////////

Lattice::Lattice(Hamiltonian h, double t, char m) :
        hamiltonian(h), temp(t), mode(m)
{
    hFunction = hamiltonian.getHamiltonian();
    indices = hamiltonian.getIndices();
    numIndices = hamiltonian.getNumIndices();
    localTerms = hamiltonian.getLocalTerms();
    indInteractions = hamiltonian.getIndInteractions();
    shape = hamiltonian.getShape();

    initSpins();
    chooseSpeed();
}

void Lattice::initSpins()
{
    srand(time(NULL));

    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        spins[*it] = (rand() % 2 == 1) ? 1 : -1;
    }
}

void Lattice::chooseSpeed()
{
    bool isFast = true;

    vector< vector<int> >::iterator it = hFunction.begin();
    char c = (*it)[0];

    for (++it; it != hFunction.end(); ++it)
    {
        if ((*it)[0] != c)
        {
            isFast = false;
            break;
        }
    }

    for (int i = 0; i < numIndices; ++i)
    {
        if (indices[i] != i)
        {
            isFast = false;
            break;
        }
    }

    if (isFast)
    {
        setCoupling(c);
        setSpeedFast();
    }
    else
    {
        setSpeedNormal();
    }
}

void Lattice::setSpeedNormal()
{
    updateAllPtr = &Lattice::updateAll;
    updatePseudoPtr = &Lattice::updatePseudo;
    updateRandomPtr = &Lattice::updateRandom;
    findTotalEnergyPtr = &Lattice::findTotalEnergy;
    findIndexEnergyPtr = &Lattice::findIndexEnergy;
    findMagnetismPtr = &Lattice::findMagnetism;
}

void Lattice::setSpeedFast()
{
    updateAllPtr = &Lattice::updateAllFast;
    updatePseudoPtr = &Lattice::updatePseudo;
    updateRandomPtr = &Lattice::updateRandomFast;
    findTotalEnergyPtr = &Lattice::findTotalEnergyFast;
    findIndexEnergyPtr = &Lattice::findIndexEnergyFast;
    findMagnetismPtr = &Lattice::findMagnetismFast;
}

void Lattice::updateLattice() {
    switch (mode)
    {
        case ALL: (this->*updateAllPtr)(); break;
        case PSEUDO: (this->*updatePseudoPtr)(); break;
        case RANDOM: (this->*updateRandomPtr)(); break;
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

void Lattice::updateAllFast()
{
    for (int i = 0; i < numIndices; ++i)
    {
        if (findProbability(i) > (double) rand() / (double) RAND_MAX)
        {
            spins[i] *= -1;
        }
    }
}

void Lattice::updatePseudo()
{
    random_shuffle(indices.begin(), indices.end());

    for (int i = 0; i < numIndices; ++i)
    {
        if (findProbability(indices[i]) > (double) rand() / (double) RAND_MAX)
        {
            spins[indices[i]] *= -1;
        }
    }
}

void Lattice::updateRandom()
{
    int index;

    for (int i = 0; i < numIndices; ++i)
    {
        index = indices[rand() % numIndices];

        if (findProbability(index) > (double) rand() / (double) RAND_MAX)
        {
            spins[index] *= -1;
        }
    }
}

void Lattice::updateRandomFast()
{
    int index;

    for (int i = 0; i < numIndices; ++i)
    {
        index = rand() % numIndices;

        if (findProbability(index) > (double) rand() / (double) RAND_MAX)
        {
            spins[index] *= -1;
        }
    }
}

void Lattice::switchMode(char m)
{
    switch (m)
    {
        case ALL:
            setMode(ALL);
            break;
        case PSEUDO:
            setMode(PSEUDO);
            break;
        case RANDOM:
            setMode(RANDOM);
            break;
        default:
            cout << "INVALID MODE. Exiting...\n\n";
            exit(EXIT_FAILURE);
    }
}

double Lattice::findProbability(int index)
{
    int initEnergy = (this->*findIndexEnergyPtr)(index);
    spins[index] *= -1;

    int finalEnergy = (this->*findIndexEnergyPtr)(index);
    spins[index] *= -1;

    if (finalEnergy > initEnergy)
    {
        return pow(E, (-1 / getTemp()) * (finalEnergy - initEnergy));
    }
    else
    {
        return 1;
    }
}

int Lattice::findTotalEnergy()
{
    int energy = 0;
    vector<int>::iterator it;

    for (it = indices.begin(); it != indices.end(); ++it)
    {
        energy += (this->*findIndexEnergyPtr)(*it);
    }

    return energy;
}

int Lattice::findTotalEnergyFast()
{
    int energy = 0;

    for (int i = 0; i < numIndices; ++i)
    {
        energy += (this->*findIndexEnergyPtr)(i);
    }

    return energy;
}

int Lattice::findIndexEnergy(int index)
{
    int energy = 0;

    vector< vector<int> >::iterator it1;
    vector<int>::iterator it2;

    for (it1 = indInteractions[index].begin();
            it1 != indInteractions[index].end(); ++it1)
    {
        it2 = it1->begin();
        int couplingEnergy = *it2;

        for (++it2; it2 != it1->end(); ++it2)
        {
            couplingEnergy *= spins[*it2];
        }

        energy -= couplingEnergy;
    }

    return spins[index] * energy;
}

int Lattice::findIndexEnergyFast(int index)
{
    int energy = 0;

    for (int i = localTerms[index].size() - 1; i >= 0; --i)
    {
        energy -= spins[localTerms[index][i]];
    }

    return coupling * spins[index] * energy;
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

double Lattice::findMagnetismFast()
{
    double magnetism = 0;

    for (int i = 0; i < numIndices; ++i)
    {
        magnetism += spins[i];
    }

    return magnetism / numIndices;
}

void Lattice::shapeError() const
{
    cout << "Wrong lattice type. Aborting...\n\n";
    exit(EXIT_FAILURE);
}

void Lattice::printLattice(int cols)
{
    if (cols == -1)
    {
        cols = (int) sqrt(numIndices);
    }

    vector<int>::iterator it = indices.begin();

    while (it != indices.end())
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

RectangularLattice::RectangularLattice(Hamiltonian h, double t, char m,
        int r, int c) : Lattice(h, t, m), rows(r), cols(c)
{
    checkShape();

    if (rows == -1 || cols == -1)
    {
        guessRowsCols();
    }
}

void RectangularLattice::checkShape() const
{
    if (getShape() != RECTANGLE && getShape() != SQUARE)
    {
        cout << "Invalid shape parameter -- not a rectangle ('r')!\n";
        shapeError();
    }
}

void RectangularLattice::guessRowsCols()
{
    if (hamiltonian.getRows() != -1 && hamiltonian.getCols() != -1)
    {
        setRows(hamiltonian.getRows());
        setCols(hamiltonian.getCols());
        return;
    }

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


SquareLattice::SquareLattice(Hamiltonian h, double t, char m, int s) :
        RectangularLattice(h, t, m, s, s), side(s)
{
    checkShape();

    if (getSide() == -1)
    {
        guessSide();
    }
}

void SquareLattice::checkShape() const
{
    if (getShape() != SQUARE)
    {
        cout << "Invalid shape parameter -- not a square ('s')!\n";
        shapeError();
    }
}

void SquareLattice::guessSide()
{
    int guess = hamiltonian.getRows();

    if (guess == -1)
    {
        int guess = (int) sqrt(numIndices);

        if (guess * guess != numIndices)
        {
            cout << "Lattice is not a square (invalid number of indices)!";
            shapeError();
        }
    }

    setSide(guess);
    setRows(guess);
    setCols(guess);
}


///////////////////////
// TriangularLattice //
///////////////////////


TriangularLattice::TriangularLattice(Hamiltonian h, double t, char m,
        int r, int c) : Lattice(h, t, m), rows(r), cols(c)
{
    checkShape();

    if (rows == -1 || cols == -1)
    {
        guessRowsCols();
    }
}

void TriangularLattice::checkShape() const
{
    if (getShape() != TRIANGLE && getShape() != STRIANGLE)
    {
        cout << "Invalid shape parameter -- not a triangle ('t')!\n";
        shapeError();
    }
}

void TriangularLattice::guessRowsCols()
{
    if (hamiltonian.getRows() != -1 && hamiltonian.getCols() != -1)
    {
        setRows(hamiltonian.getRows());
        setCols(hamiltonian.getCols());
        return;
    }

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


////////////////////////
// STriangularLattice //
////////////////////////


STriangularLattice::STriangularLattice(Hamiltonian h, double t, char m,
        int s) : TriangularLattice(h, t, m, s, s), side(s)
{
    checkShape();

    if (getSide() == -1)
    {
        guessSide();
    }
}

void STriangularLattice::checkShape() const
{
    if (getShape() != STRIANGLE)
    {
        cout << "Invalid shape parameter -- not a square triangle ('v')!\n";
        shapeError();
    }
}

void STriangularLattice::guessSide()
{
    int guess = hamiltonian.getRows();

    if (guess == -1)
    {
        int guess = (int) sqrt(numIndices);

        if (guess * guess != numIndices)
        {
            cout << "Lattice is not a square triangle ";
            cout << "(invalid number of indices)!";
            shapeError();
        }
    }

    setSide(guess);
    setRows(guess);
    setCols(guess);
}
