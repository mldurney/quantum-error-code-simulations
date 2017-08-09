#include "lattices.h"

using namespace ising;

/////////////
// Lattice //
/////////////

Lattice::Lattice(Hamiltonian h, double t, char m, bool init)
    : hamiltonian(h),
      numIndices(h.getNumIndices()),
      shape(h.getShape()),
      temp(t),
      mode(m) {
    mapsToSequences();
    setType("default");
    setSize((int)sqrt(getNumIndices()));
    srand((unsigned int)time(NULL));
    zSeed = rand();
    wSeed = rand();
    randomizedIndices = indices;
    generateDistances();

    if (init) {
        initSpins();
    }
}

void Lattice::mapsToSequences() {
    // Create map from original indices to new sequential locations
    auto origIndices = hamiltonian.getIndices();
    std::map<int, int> indMap;
    for (int i = 0; i < numIndices; ++i) {
        indMap[origIndices[i]] = i;
    }

    // Initialize sequential indices
    indices.resize(numIndices);

    // Initialize sequential Hamiltonian function
    auto origHFunction = hamiltonian.getHamiltonian();
    for (auto& origInteraction : origHFunction) {
        ivector interaction;
        bool isFirst = true;
        for (auto& i : origInteraction) {
            if (isFirst) {
                interaction.push_back(i);
                isFirst = false;
                continue;
            }

            interaction.push_back(indMap[i]);
        }

        hFunction.push_back(interaction);
    }

    // Initialize sequential locations
    locations.resize(numIndices);
    auto origLocations = hamiltonian.getLocations();
    for (auto& loc : origLocations) {
        locations[indMap[loc.first]] = loc.second;
    }

    // Initialize sequential index interactions
    indInteractions.resize(numIndices);
    auto origIndInteractions = hamiltonian.getIndInteractions();
    for (auto& interaction : origIndInteractions) {
        indInteractions[indMap[interaction.first]] = interaction.second;
    }
}

void Lattice::initSpins() {
    spinReplicas.resize(REPLICAS);

    for (auto& replica : spinReplicas) {
        for (auto& i : replica) {
            i = (MWC() % 2 == 0) ? 1 : -1;
        }
    }
}

void Lattice::generateDistances() {
    for (auto& i : indices) {
        for (auto& j : indices) {
            xDisplacements[i][j] = findXDisplacement(i, j);
            yDisplacements[i][j] = findYDisplacement(i, j);
            distances[i][j] = findDistance(i, j);
        }
    }
}

void Lattice::metropolisUpdate(int replicaNum) {
    switch (mode) {
        case ALL:
            updateAll(replicaNum);
            break;
        case PSEUDO:
            updatePseudo(replicaNum);
            break;
        case RANDOM:
            updateRandom(replicaNum);
            break;
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

void Lattice::updateAll(int replicaNum) {
    for (int i = 0; i < numIndices; ++i) {
        if (findProbability(replicaNum, i) > randFloatCO()) {
            spinReplicas[replicaNum][i] *= -1;
        }
    }
}

void Lattice::updatePseudo(int replicaNum) {
    for (int i = 0; i < numIndices; ++i) {
        int j = MWC() % numIndices;
        std::swap(randomizedIndices[i], randomizedIndices[j]);
    }

    for (int i = 0; i < numIndices; ++i) {
        int index = randomizedIndices[i];
        if (findProbability(replicaNum, index) > randFloatCO()) {
            spinReplicas[replicaNum][index] *= -1;
        }
    }
}

void Lattice::updateRandom(int replicaNum) {
    int index;

    for (int i = 0; i < numIndices; ++i) {
        index = indices[MWC() % numIndices];

        if (findProbability(replicaNum, index) > randFloatCO()) {
            spinReplicas[replicaNum][index] *= -1;
        }
    }
}

double Lattice::findProbability(int replicaNum, int index) {
    int initEnergy = findIndexEnergy(replicaNum, index);

    if (initEnergy < 0) {
        return pow(E, (-1 / getTemp()) * (-2 * initEnergy));
    } else {
        return 1;
    }
}

int Lattice::findTotalEnergy(int replicaNum) {
    int energy = 0;

    for (int i = 0; i < numIndices; ++i) {
        energy += findIndexEnergy(replicaNum, i);
    }

    return energy;
}

int Lattice::findIndexEnergy(int index, int replicaNum) {
    int energy = 0;

    ivector2::const_iterator it1;
    ivector::const_iterator it2;

    auto end = indInteractions[index].end();
    for (it1 = indInteractions[index].begin(); it1 != end; ++it1) {
        it2 = it1->begin();
        int couplingEnergy = *it2;

        for (++it2; it2 != it1->end(); ++it2) {
            couplingEnergy *= spinReplicas[replicaNum][*it2];
        }

        energy -= couplingEnergy;
    }

    return spinReplicas[replicaNum][index] * energy;
}

double Lattice::findMagnetization(int replicaNum) {
    double magnetism = 0;

    for (int i = 0; i < numIndices; ++i) {
        magnetism += spinReplicas[replicaNum][i];
    }

    return magnetism / numIndices;
}

void Lattice::shapeError() const {
    std::cout << "Wrong lattice type. Aborting...\n\n";
    exit(EXIT_FAILURE);
}

void Lattice::flipSpins(int replicaNum) {
    for (auto& s : spinReplicas[replicaNum]) {
        s *= -1;
    }
}

void Lattice::print(int replicaNum, int cols) const {
    if (cols == -1) {
        cols = (int)sqrt(numIndices);
    }

    for (int col = 0, i = 0; col < cols && i < numIndices; ++col, ++i) {
        (spinReplicas[replicaNum][i] == 1) ? std::cout << "+ "
                                           : std::cout << "- ";
    }

    std::cout << "\n";
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

////////////////////////
// RectangularLattice //
////////////////////////

RectangularLattice::RectangularLattice(Hamiltonian h, double t, char m,
                                       bool init, int r, int c)
    : Lattice(h, t, m, init), rows(r), cols(c) {
    setType("rectangle");
    checkShape();

    if (rows == -1 || cols == -1) {
        guessRowsCols();
    }

    generateDistances();
}

void RectangularLattice::checkShape() const {
    if (getShape() != RECTANGLE && getShape() != SQUARE) {
        std::cout << "Invalid shape parameter -- not a rectangle ('r')!\n";
        shapeError();
    }
}

int RectangularLattice::findXDisplacement(int i, int j) {
    int iRow = locations[i][0];
    int jRow = locations[j][0];
    return iRow - jRow;
}

int RectangularLattice::findYDisplacement(int i, int j) {
    int iCol = locations[i][1];
    int jCol = locations[i][1];
    return iCol - jCol;
}

double RectangularLattice::findDistance(int i, int j) {
    return sqrt(pow(xDisplacements[i][j], 2) + pow(yDisplacements[i][j], 2));
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

SquareLattice::SquareLattice(Hamiltonian h, double t, char m, bool init, int s)
    : Lattice(h, t, m, init), RectangularLattice(h, t, m, init, s, s), side(s) {
    setType("square");
    checkShape();

    if (getSide() == -1) {
        guessSide();
    }

    setSize(getSide());
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

TriangularLattice::TriangularLattice(Hamiltonian h, double t, char m, bool init,
                                     int r, int c)
    : Lattice(h, t, m, init), rows(r), cols(c) {
    setType("triangle");
    checkShape();

    if (rows == -1 || cols == -1) {
        guessRowsCols();
    }

    generateDistances();
}

void TriangularLattice::checkShape() const {
    if (getShape() != TRIANGLE && getShape() != STRIANGLE) {
        std::cout << "Invalid shape parameter -- not a triangle ('t')!\n";
        shapeError();
    }
}

int TriangularLattice::findXDisplacement(int i, int j) {
    int iCol = locations[i][0];
    int jCol = locations[j][0];
    return iCol - jCol;
}

int TriangularLattice::findYDisplacement(int i, int j) {
    int iCol = locations[i][1];
    int jCol = locations[j][1];
    return iCol - jCol;
}

double TriangularLattice::findDistance(int i, int j) {
    return sqrt(pow(xDisplacements[i][j], 2) + pow(yDisplacements[i][j], 2));
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

STriangularLattice::STriangularLattice(Hamiltonian h, double t, char m,
                                       bool init, int s)
    : Lattice(h, t, m, init), TriangularLattice(h, t, m, init, s, s), side(s) {
    setType("square triangle");
    checkShape();

    if (getSide() == -1) {
        guessSide();
    }

    setSize(getSide());
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
