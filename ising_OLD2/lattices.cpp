#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <cassert>
#include "lattices.h"


//////////////////////////////////
// Generic Lattice constructors //
//////////////////////////////////

/**
 * [Lattice::Lattice   Construct generic lattice with same coupling throughout]
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {int}       Coupling between all pairs of indices
 */
Lattice::Lattice(double t, char m, int c) :
        Lattice(t, m, vector<int>(1, c), true) {}

/**
 * [Lattice::Lattice   Construct generic lattice where coupling between pairs
 *                     cycles through an int vector of options]
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {vector<int>}   List of couplings cycled through and assigned
 *                              to each index pair
 */
Lattice::Lattice(double t, char m, vector<int> c = vector<int>(1, 1),
        bool constCoupling = false) :
        temp(t), mode(m), isCouplingConst(constCoupling)
{
    couplings = c;

    if (isCouplingConst)
    {
        coupling = couplings[0];
        findIndexEnergyPtr = &Lattice::findIndexEnergyQuick;
    }
    else
    {
        findIndexEnergyPtr = &Lattice::findIndexEnergy;
    }
}


//////////////////////////////////////////
// Define protected vectors of Lattice  //
//////////////////////////////////////////

/**
 * [Lattice::initSpins  Define int vector of spins, randomly assigning state
 *                      1 or -1 separately to n spins, where n is the number
 *                      of indices in the lattice]
 *
 *                      Give starting spins for Lattice
 */
void Lattice::initSpins()
{
    srand(time(NULL));

    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    for (itRow = indices.begin(); itRow != indices.end(); ++itRow)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol)
        {
            spins.push_back((rand() % 2 == 1) ? 1 : -1);
        }
    }
}

/**
 * [Lattice::initHamiltonian    Define vector of vectors of size 2 arrays
 *                              The outer vector represents each index,
 *                              the inner vectors represent each interaction
 *                              with that index, and each int array gives
 *                              a specific interaction {other_index, coupping}]
 *
 *                              Give hamiltonian for Lattice
 */
void Lattice::initHamiltonian()
{
    vector<vector<int>>::iterator itLocal;
    vector<int>::iterator itTerms;
    vector<int>::iterator itCoup = couplings.begin();

    for (itLocal = localTerms.begin(); itLocal != localTerms.end(); ++itLocal)
    {
        vector<array<int, 2>> indexHamiltonian;

        for (itTerms = itLocal->begin(); itTerms != itLocal->end();
                ++itTerms, ++itCoup)
        {
            if (itCoup == couplings.end())
            {
                itCoup = couplings.begin();
            }

            int currCoupling = *itCoup;
            int currTerm = *itTerms;

            indexHamiltonian.push_back({currCoupling, currTerm});
        }

        hamiltonian.push_back(indexHamiltonian);
    }
}


//////////////////////////////////
// Functions to update Lattice  //
//////////////////////////////////

void Lattice::updateLattice() {
    switch (mode)
    {
        case ALL: updateAll(); break;
        case RANDOM: updateRandom(); break;
    }
}

/**
 * [Lattice::updateAll  If ALL mode active, iterate through each index and
 *                      update it probabilitistically based on initial and
 *                      final energy of local region]
 */
void Lattice::updateAll()
{
    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    for (itRow = indices.begin(); itRow != indices.end(); ++itRow)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol)
        {
            if (findProbability(*itCol) > (double) rand() / (double) RAND_MAX)
            {
                spins[*itCol] *= -1;
            }
        }
    }
}

/**
 * [Lattice::updateRandom   If RANDOM mode active, randomly select n indices,
 *                          where n is the total number of indices, and update
 *                          each selected indices probabilitistically based on
 *                          initial and final energy of local region]
 */
void Lattice::updateRandom()
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

/**
 * [Lattice::switchMode     Switch update mode]
 * @param mode  {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 */
void Lattice::switchMode(char mode)
{
    switch (mode)
    {
        case ALL: mode = ALL; break;
        case RANDOM: mode = RANDOM; break;
        default: printf("\nINVALID MODE. Exiting...\n"); exit(EXIT_FAILURE);
    }
}

/**
 * [Lattice::findProbability    Determine chance that index spin will flip in
 *                              update]
 * @param  index    {int}       Integer representation of index location
 * @return          {double}    Probability that index will flip, from 0.0 to
 *                              1.0
 */
double Lattice::findProbability(int index)
{
    double initEnergy = (this->*findIndexEnergyPtr)(index);
    spins[index] *= -1;

    double finalEnergy = (this->*findIndexEnergyPtr)(index);
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


////////////////////////////////////////////
// Functions to find energy and magnetism //
////////////////////////////////////////////

/**
 * [Lattice::findTotalEnergy    Find total energy for current Lattice]
 * @return [description]
 */
double Lattice::findTotalEnergy()
{
    double energy = 0;
    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    for (itRow = indices.begin(); itRow != indices.end(); ++itRow)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol)
        {
            energy += (this->*findIndexEnergyPtr)(*itCol);
        }
    }

    return energy;
}

double Lattice::findIndexEnergy(int index)
{
    double energy = 0;
    vector<array<int, 2>>::iterator it;

    for (it = hamiltonian[index].begin(); it != hamiltonian[index].end(); ++it)
    {
        energy -= (*it)[0] * spins[(*it)[1]];
    }

    return spins[index] * energy;
}

double Lattice::findIndexEnergyQuick(int index)
{
    double energy = 0;

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

    for (it = spins.begin(); it != spins.end(); ++it)
    {
        magnetism += *it;
    }

    return magnetism / numIndices;
}


/////////////////////////////////////
// Miscellaneous Lattice functions //
/////////////////////////////////////

/**
 * [Lattice::printLattice   Print Lattice to terminal
 *                          1 spins represented by + and -1 spins by - in grid]
 */
void Lattice::printLattice()
{
    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    for (itRow = indices.begin(); itRow != indices.end(); ++itRow)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol)
        {
            (spins[*itCol] == 1) ? printf("+ ") : printf("- ");
        }

        printf("\n");
    }
}


/////////////////////////////////////
// RectangularLattice constructors //
/////////////////////////////////////

/**
 * [RectangularLattice::RectangularLattice  Construct rectangular lattice
 *                                          Derived from generic Lattice
 *                                          Coupling same throughout]
 * @param row   {int}       Number of rows in rectangle
 * @param col   {int}       Number of columns in rectangle
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {int}       Coupling between all pairs of indices
 */
RectangularLattice::RectangularLattice(int row, int col, double t, char m,
        int c) : RectangularLattice(row, col, t, m, vector<int>(1, c), true) {}

/**
 * [RectangularLattice::RectangularLattice  Construct rectangular lattice
 *                                          Derived from generic Lattice
 *                                          Coupling between pairs cycles
 *                                          through int vector of options]
 * @param row   {int}       Number of rows in rectangle
 * @param col   {int}       Number of columns in rectangle
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {vector<int>}   List of couplings cycled through and assigned
 *                              to each index pair
 */
RectangularLattice::RectangularLattice(int row, int col, double t, char m,
        vector<int> c, bool constCoupling = false) :
        Lattice(t, m, c, constCoupling), rows(row), cols(col)
{
    int index = 0;

    for (int i = 0; i < rows; ++i)
    {
        vector<int> currRow;

        for (int j = 0; j < cols; ++j, ++index)
        {
            currRow.push_back(index);
        }

        indices.push_back(currRow);
    }

    initSpins();
    initLocalTerms();
    initHamiltonian();

    numIndices = spins.size();
}


//////////////////////////////////////////////////////////////////////////
// Specialized functions to define vector members of RectangularLattice //
//////////////////////////////////////////////////////////////////////////

/**
 * [RectangularLattice::initLocalTerms  Define vector of int vectors,
 *                                      where outer vector corresponds
 *                                      to a specific index and inner vector
 *                                      lists other int indices that
 *                                      interact with index
 *                                      For RectangularLattice, each index on
 *                                      a vertex with all other indices on
 *                                      vertices connected to original vertex
 *                                      by edge
 *
 *                                      Give localTerms for RectangularLattice]
 */
void RectangularLattice::initLocalTerms()
{
    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    int row = 0, col = 0;
    for (itRow = indices.begin(); itRow != indices.end(); ++itRow, ++row)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol, ++col)
        {
            int left = row * cols + ((col != 0) ? col - 1 : cols - 1);
            int right = row * cols + ((col != cols - 1) ? col + 1 : 0);
            int top = col + cols * ((row != 0) ? row - 1 : rows - 1);
            int bottom = col + cols * ((row != rows - 1) ? row + 1 : 0);

            localTerms.push_back({left, right, top, bottom});
        }
    }
}


////////////////////////////////
// SquareLattice constructors //
////////////////////////////////

/**
 * [SquareLattice::SquareLattice    Construct square lattice
 *                                  Derived from generic Lattice
 *                                  Coupling same throughout]
 * @param size   {int}      Side length
 *                              Number of rows and columns in square
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {int}       Coupling between all pairs of indices
 */
SquareLattice::SquareLattice(int side, double t, char m, int c) :
        RectangularLattice(side, side, t, m, vector<int>(1, c), true) {}

/**
 * [SquareLattice::SquareLattice    Construct square lattice
 *                                  Derived from generic Lattice
 *                                  Coupling between pairs cycles through
 *                                  int vector of options]
 * @param size   {int}      Side length
 *                              Number of rows and columns in square
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {vector<int>}   List of couplings cycled through and assigned
 *                              to each index pair
 */
SquareLattice::SquareLattice(int side, double t, char m, vector<int> c) :
        RectangularLattice(side, side, t, m, c)
{
    initSpins();
    initLocalTerms();
    initHamiltonian();
}


////////////////////////////////////
// TriangularLattice constructors //
////////////////////////////////////

/**
 * [TriangularLattice::TriangularLattice    Construct triangular lattice
 *                                          Derived from generic Lattice
 *                                          Coupling same throughout]
 * @param row   {int}       Number of rows in rectangle
 * @param col   {int}       Number of columns in rectangle
 *                              Two traingles fit in each cell where row and
 *                              column intersect
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {int}       Coupling between all pairs of indices
 */
TriangularLattice::TriangularLattice(int row, int col, double t, char m,
        int c) : TriangularLattice(row, col, t, m, vector<int>(1, c), true) {}

/**
 * [TriangularLattice::TriangularLattice    Construct triangular lattice
 *                                          Derived from generic Lattice
 *                                          Coupling between pairs cycles
 *                                          through int vector of options]
 * @param row   {int}       Number of rows in rectangle
 * @param col   {int}       Number of columns in rectangle
 *                              Two traingles fit in each cell where row and
 *                              column intersect
 * @param t     {double}    Temperature of lattice
 * @param m     {char}      Mode of updating lattice
 *                              'a' (ALL) - Every index updated sequentially
 *                              'r' (RANDOM) - n indices randomly selected and
 *                                             updated, where n is the total
 *                                             number of indices
 * @param c     {vector<int>}   List of couplings cycled through and assigned
 *                              to each index pair
 */
TriangularLattice::TriangularLattice(int row, int col, double t, char m,
        vector<int> c = vector<int>(1, 1), bool constCoupling = false) :
        Lattice(t, m, c, constCoupling), rows(row), cols(col)
{
    int index = 0;

    for (int i = 0; i < rows; ++i)
    {
        vector<int> currRow;

        for (int j = 0; j < cols * 2; ++j, ++index)
        {
            currRow.push_back(index);
        }

        indices.push_back(currRow);
    }

    initSpins();
    initLocalTerms();
    initHamiltonian();

    numIndices = spins.size();
}


/////////////////////////////////////////////////////////////////////////
// Specialized functions to define vector members of TriangularLattice //
/////////////////////////////////////////////////////////////////////////

/**
 * [TriangularLattice::initLocalTerms   Define vector of int vectors,
 *                                      where outer vector corresponds
 *                                      to a specific index and inner vector
 *                                      lists other int indices that
 *                                      interact with index
 *                                      For TriangularLattice, each index on
 *                                      a vertex with all other indices on
 *                                      vertices connected to original vertex
 *                                      by edge
 *
 *                                      Give localTerms for TriangularLattice]
 */
void TriangularLattice::initLocalTerms()
{
    vector<vector<int>>::iterator itRow;
    vector<int>::iterator itCol;

    int row, col= 0;
    for (itRow = indices.begin(); itRow != indices.end(); ++itRow, ++row)
    {
        for (itCol = itRow->begin(); itCol != itRow->end(); ++itCol, ++col)
        {
            int left = row * cols + ((col != 0) ? col - 1 : cols - 1);
            int right = row * cols + ((col != cols - 1) ? col + 1 : 0);
            int top = col + cols * ((row != 0) ? row - 1 : rows - 1);
            int bottom = col + cols * ((row != rows - 1) ? row + 1 : 0);
            int topleft = ((col != 0) ? col - 1 : cols - 1)
                    + cols * ((row != 0) ? row - 1 : rows - 1);
            int bottomright = row * cols + ((col != cols - 1) ? col + 1 : 0)
                    + col + cols * ((row != rows - 1) ? row + 1 : 0);

            localTerms.push_back({left, right, top, bottom,
                    topleft, bottomright});
        }
    }
}
