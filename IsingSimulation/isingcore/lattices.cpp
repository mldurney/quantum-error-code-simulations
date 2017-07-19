#include "lattices.h"

using namespace ising;

/////////////
// Lattice //
/////////////

Lattice::Lattice(Hamiltonian h, double t, char m)
    : hamiltonian(h),
      hFunction(h.getHamiltonian()),
      indices(h.getIndices()),
      numIndices(h.getNumIndices()),
      localTerms(h.getLocalTerms()),
      indInteractions(h.getIndInteractions()),
      shape(h.getShape()),
      temp(t),
      mode(m) {
    srand((unsigned int)time(NULL));
    zSeed = rand();
    wSeed = rand();

    randomizedIndices = indices;
    initSpins();
	const std::map<int, std::map<int, double>> distances = generateDistances();
}

void Lattice::initSpins() {
    spins = std::vector<int>(numIndices, 0);

    std::vector<int>::const_iterator it;

    for (it = indices.begin(); it != indices.end(); ++it) {
        spins[*it] = (MWC() % 2 == 0) ? 1 : -1;
    }
}

std::map<int, std::map<int, double>> Lattice::generateDistances() {
	std::map<int, std::map<int, double>> dists;
	for (auto i : indices) {
		for (auto j : indices) {
			dists[i][j] = findDistance(i, j);
		}
	}

	return dists;
}

void Lattice::updateLattice() {
    switch (mode) {
        case ALL:
            updateAll();
            break;
        case PSEUDO:
            updatePseudo();
            break;
        case RANDOM:
            updateRandom();
            break;
    }
}

void Lattice::updateAll() {
    std::vector<int>::const_iterator it;

    for (it = indices.begin(); it != indices.end(); ++it) {
        if (findProbability(*it) > randFloatCO()) {
            spins[*it] *= -1;
        }
    }
}

void Lattice::updatePseudo() {
    for (int i = 0; i < numIndices; ++i) {
        int j = MWC() % numIndices;
        std::swap(randomizedIndices[i], randomizedIndices[j]);
    }

    for (int i = 0; i < numIndices; ++i) {
        if (findProbability(randomizedIndices[i]) > randFloatCO()) {
            spins[randomizedIndices[i]] *= -1;
        }
    }
}

void Lattice::updateRandom() {
    int index;

    for (int i = 0; i < numIndices; ++i) {
        index = indices[MWC() % numIndices];

        if (findProbability(index) > randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

void Lattice::switchMode(char m) {
    switch (m) {
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
            std::cout << "INVALID MODE. Exiting...\n\n";
            exit(EXIT_FAILURE);
    }
}

double Lattice::findProbability(int index) {
    int initEnergy = findIndexEnergy(index);

    if (initEnergy < 0) {
        return pow(E, (-1 / getTemp()) * (-2 * initEnergy));
    } else {
        return 1;
    }
}

int Lattice::findTotalEnergy() {
    int energy = 0;
    std::vector<int>::const_iterator it;

    for (it = indices.begin(); it != indices.end(); ++it) {
        energy += findIndexEnergy(*it);
    }

    return energy;
}

int Lattice::findIndexEnergy(int index) {
    /*
    return -spins[index] * std::accumulate(indInteractions.at(index).begin(),
            indInteractions.at(index).end(), 0,
            [&](int lhs, const auto& rhs) { return lhs +
            std::accumulate(rhs.begin() + 1, rhs.end(), rhs[0],
            [&](int a, int b) { return a * spins[b]; }); });
    */

    int energy = 0;

    std::vector<std::vector<int>>::const_iterator it1;
    std::vector<int>::const_iterator it2;

    auto end = indInteractions.at(index).end();
    for (it1 = indInteractions.at(index).begin(); it1 != end; ++it1) {
        it2 = it1->begin();
        int couplingEnergy = *it2;

        for (++it2; it2 != it1->end(); ++it2) {
            couplingEnergy *= spins[*it2];
        }

        energy -= couplingEnergy;
    }

    return spins[index] * energy;
}

double Lattice::findMagnetism() {
    double magnetism = 0;
    std::vector<int>::const_iterator it;

    for (it = indices.begin(); it != indices.end(); ++it) {
        magnetism += spins[*it];
    }

    return magnetism / numIndices;
}

void Lattice::shapeError() const {
    std::cout << "Wrong lattice type. Aborting...\n\n";
    exit(EXIT_FAILURE);
}

void Lattice::printLattice(int cols) const {
    if (cols == -1) {
        cols = (int)sqrt(numIndices);
    }

    std::vector<int>::const_iterator it = indices.begin();

    while (it != indices.end()) {
        for (int col = 0; col < cols && it != indices.end(); ++col, ++it) {
            (spins[*it] == 1) ? std::cout << "+ " : std::cout << "- ";
        }

        std::cout << "\n";
    }
}

float Lattice::asFloat(uint32_t i) {
    union {
        uint32_t i;
        float f;
    } pun = {i};
    return pun.f;
}

float Lattice::randFloatCO() {
    return asFloat(0x3F800000U | (MWC() >> 9)) - 1.0f;
}

/////////////////
// LatticeFast //
/////////////////

void LatticeFast::updateAll() {
    for (int i = 0; i < numIndices; ++i) {
        if (findProbability(i) > randFloatCO()) {
            spins[i] *= -1;
        }
    }
}

void LatticeFast::updateRandom() {
    int index;

    for (int i = 0; i < numIndices; ++i) {
        index = MWC() % numIndices;

        if (findProbability(index) > randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

int LatticeFast::findTotalEnergy() {
    int energy = 0;

    for (int i = 0; i < numIndices; ++i) {
        energy += findIndexEnergy(i);
    }

    return energy;
}

int LatticeFast::findIndexEnergy(int index) {
    int energy = 0;

    for (int i = (int)localTerms.at(index).size() - 1; i >= 0; --i) {
        energy -= spins[localTerms.at(index)[i]];
    }

    return coupling * spins[index] * energy;
}

double LatticeFast::findMagnetism() {
    double magnetism = 0;

    for (int i = 0; i < numIndices; ++i) {
        magnetism += spins[i];
    }

    return magnetism / numIndices;
}

////////////////////////
// RectangularLattice //
////////////////////////

RectangularLattice::RectangularLattice(Hamiltonian h, double t, char m, int r,
                                       int c)
    : Lattice(h, t, m), rows(r), cols(c) {
    checkShape();

    if (rows == -1 || cols == -1) {
        guessRowsCols();
    }
}

void RectangularLattice::checkShape() const {
    if (getShape() != RECTANGLE && getShape() != SQUARE) {
        std::cout << "Invalid shape parameter -- not a rectangle ('r')!\n";
        shapeError();
    }
}

double RectangularLattice::findDistance(int i, int j) {
	int iRow = i % getCols();
	int iCol = i / getCols();
	int jRow = j % getCols();
	int jCol = j / getCols();

	int distRow = iRow - jRow;
	int distCol = iCol - jCol;
	return sqrt(pow(distRow, 2) + pow(distCol, 2));
}

void RectangularLattice::guessRowsCols() {
    if (hamiltonian.getRows() != -1 && hamiltonian.getCols() != -1) {
        setRows(hamiltonian.getRows());
        setCols(hamiltonian.getCols());
        return;
    }

    for (int guess = (int)sqrt(numIndices); guess > 0; --guess) {
        if (numIndices % guess == 0) {
            setRows(guess);
            setCols(numIndices / guess);
            return;
        }
    }
}

///////////////////
// SquareLattice //
///////////////////

SquareLattice::SquareLattice(Hamiltonian h, double t, char m, int s)
    : Lattice(h, t, m), RectangularLattice(h, t, m, s, s), side(s) {
    checkShape();

    if (getSide() == -1) {
        guessSide();
    }
}

void SquareLattice::checkShape() const {
    if (getShape() != SQUARE) {
        std::cout << "Invalid shape parameter -- not a square ('s')!\n";
        shapeError();
    }
}

void SquareLattice::guessSide() {
    int guess = hamiltonian.getRows();

    if (guess == -1) {
        guess = (int)sqrt(numIndices);

        if (guess * guess != numIndices) {
            std::cout << "Lattice is not a square (invalid number of indices)!";
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

TriangularLattice::TriangularLattice(Hamiltonian h, double t, char m, int r,
                                     int c)
    : Lattice(h, t, m), rows(r), cols(c) {
    checkShape();

    if (rows == -1 || cols == -1) {
        guessRowsCols();
    }
}

void TriangularLattice::checkShape() const {
    if (getShape() != TRIANGLE && getShape() != STRIANGLE) {
        std::cout << "Invalid shape parameter -- not a triangle ('t')!\n";
        shapeError();
    }
}

double TriangularLattice::findDistance(int i, int j) {
	int iRow = i % getCols();
	int iCol = i / getCols();
	int jRow = j % getCols();
	int jCol = j / getCols();

	int distRow = iRow - jRow;
	int distCol = iCol - jCol;
	return sqrt(pow(distRow, 2) + pow(distCol, 2));
}

void TriangularLattice::guessRowsCols() {
    if (hamiltonian.getRows() != -1 && hamiltonian.getCols() != -1) {
        setRows(hamiltonian.getRows());
        setCols(hamiltonian.getCols());
        return;
    }

    for (int guess = (int)sqrt(numIndices); guess > 0; --guess) {
        if (numIndices % guess == 0) {
            setRows(guess);
            setCols(numIndices / guess);
            return;
        }
    }
}

////////////////////////
// STriangularLattice //
////////////////////////

STriangularLattice::STriangularLattice(Hamiltonian h, double t, char m, int s)
    : Lattice(h, t, m), TriangularLattice(h, t, m, s, s), side(s) {
    checkShape();

    if (getSide() == -1) {
        guessSide();
    }
}

void STriangularLattice::checkShape() const {
    if (getShape() != STRIANGLE) {
        std::cout
            << "Invalid shape parameter -- not a square triangle ('v')!\n";
        shapeError();
    }
}

void STriangularLattice::guessSide() {
    int guess = hamiltonian.getRows();

    if (guess == -1) {
        guess = (int)sqrt(numIndices);

        if (guess * guess != numIndices) {
            std::cout << "Lattice is not a square triangle ";
            std::cout << "(invalid number of indices)!";
            shapeError();
        }
    }

    setSide(guess);
    setRows(guess);
    setCols(guess);
}
