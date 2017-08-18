#include "lattices.h"

using namespace ising;

/////////////
// Lattice //
/////////////

Lattice::Lattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m)
    : prop(h, t, dt, n, m) {
    mapsToSequences();
    setType("default");
    setSize((int)sqrt(getNumIndices()));
    setJTemperature(t + n / 2 * dt);
    generateDistances();
    gen = RandomGenerator();

    for (uint i = 0; i < prop.numT; ++i) {
        replicavector replicas;
        for (uint j = 0; j < REPLICAS; ++j) {
            replicas.emplace_back(std::make_unique<Replica>(prop, i));
        }
        configs.push_back(std::move(replicas));
        replicaIndices.push_back(i);
    }

    std::sort(replicaIndices.begin(), replicaIndices.end());
}

void Lattice::mapsToSequences() {
    // Create map from original indices to new sequential locations
    auto origIndices = prop.hamiltonian.getIndices();
    std::map<int, int> indMap;
    for (uint i = 0; i < prop.numIndices; ++i) {
        indMap[origIndices[i]] = i;
    }

    // Initialize sequential indices
    for (unsigned i = 0; i < prop.numIndices; ++i) {
        prop.indices.push_back(i);
    }

    // Initialize sequential Hamiltonian function
    auto origHFunction = prop.hamiltonian.getHamiltonian();
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

        prop.hFunction.push_back(interaction);
    }

    // Initialize sequential locations
    prop.locations.resize(prop.numIndices);
    auto origLocations = prop.hamiltonian.getLocations();
    for (auto& loc : origLocations) {
        prop.locations[indMap[loc.first]] = loc.second;
    }

    // Initialize sequential index local terms
    prop.localTerms.resize(prop.numIndices);
    auto origLocalTerms = prop.hamiltonian.getLocalTerms();
    for (auto& terms : origLocalTerms) {
        prop.localTerms[indMap[terms.first]] = terms.second;
    }

    // Initialize sequential index interactions
    prop.indInteractions.resize(prop.numIndices);
    auto origIndInteractions = prop.hamiltonian.getIndInteractions();
    for (auto& interaction : origIndInteractions) {
        prop.indInteractions[indMap[interaction.first]] = interaction.second;
    }
}

void Lattice::generateDistances() {
    prop.xDisplacements.resize(prop.numIndices);
    prop.yDisplacements.resize(prop.numIndices);
    prop.distances.resize(prop.numIndices);

    for (uint i = 0; i < prop.numIndices; ++i) {
        prop.xDisplacements[i].resize(prop.numIndices);
        prop.yDisplacements[i].resize(prop.numIndices);
        prop.distances[i].resize(prop.numIndices);
    }

    for (auto& i : prop.indices) {
        for (auto& j : prop.indices) {
            prop.xDisplacements[i][j] = findXDisplacement(i, j);
            prop.yDisplacements[i][j] = findYDisplacement(i, j);
            prop.distances[i][j] = findDistance(i, j);
        }
    }
}

Replica Lattice::getReplicaCopy(uint i, uint j) {
    if (i >= configs.size() && j > REPLICAS) {
        std::cout << "INVALID CONFIGURATION/REPLICA INDICES! Exiting...\n\n";
        exit(EXIT_FAILURE);
    }

    return *configs[i][j].get();
}

void Lattice::monteCarloSweep() {
    for (auto& replicas : configs) {
        for (auto& replica : replicas) {
            replica->update();
        }
    }
}

void Lattice::houdayerClusterMove() {
    for (uint i = 0; i < configs.size(); ++i) {
        houdayerClusterMove(i);
    }
}

void Lattice::houdayerClusterMove(uint index) {
    auto spins0 = configs[index][0]->getSpins();
    auto spins1 = configs[index][1]->getSpins();

    bool done = false;
    while (!done) {
        ivector qIndices;
        for (uint i = 0; i < prop.numIndices; ++i) {
            spins0[i] *= spins1[i];
            if (spins0[i] == 0) {
                qIndices.push_back(i);
            }
        }

        if (qIndices.size() == 0) {
            done = true;
        } else if (2 * qIndices.size() > prop.numIndices) {
            configs[index][0]->flipSpins();
        } else {
            ivector cluster = {qIndices[gen.MWC() % qIndices.size()]};
            ivector unbranched = cluster;

            while (!unbranched.empty()) {
                int index = unbranched.back();
                unbranched.pop_back();

                for (auto& i : prop.localTerms[index]) {
                    if (std::find(cluster.begin(), cluster.end(), i) ==
                            cluster.end() &&
                        std::binary_search(qIndices.begin(), qIndices.end(),
                                           i)) {
                        cluster.push_back(i);
                        unbranched.push_back(i);
                    }
                }
            }

            for (auto& i : cluster) {
                configs[index][0]->flipSpin(i);
                configs[index][1]->flipSpin(i);
            }

            done = true;
        }
    }
}

void Lattice::parallelTemperingUpdate() {
    for (uint i = 0; i < prop.numT - 1; ++i) {
        ldouble dEnergy = configs[i][0]->getTotalEnergy() -
                          configs[i + 1][0]->getTotalEnergy();
        ldouble dBoltzmann = 1 / (KB * configs[i][0]->getTemperature()) -
                             1 / (KB * configs[i + 1][0]->getTemperature());
        ldouble probability = dEnergy * dBoltzmann;

        if (probability > 1 || probability > gen.randFloatCO()) {
            for (uint j = 0; j < REPLICAS; ++j) {
                swapConfigs(i, i + 1);
            }
        }
    }
}

void Lattice::swapConfigs(uint i, uint j) {
    if (i == j || i >= configs.size() || j >= configs.size()) {
        std::cout << "INVALID CONFIGURATION INDEX! Exiting...\n\n";
        exit(EXIT_FAILURE);
    }

    auto t1 = configs[i][0]->getTemperature();
    auto t2 = configs[j][0]->getTemperature();

    for (uint k = 0; k < REPLICAS; ++k) {
        configs[i][k]->setTemperature(t2);
        configs[j][k]->setTemperature(t1);
        configs[i][k].swap(configs[j][k]);
    }
}

void Lattice::HCA() {
    monteCarloSweep();
    houdayerClusterMove();
    parallelTemperingUpdate();
}

void Lattice::ICA() {
    monteCarloSweep();

    for (uint i = 0; i < configs.size(); ++i) {
        ldouble t = configs[i][0]->getTemperature();
        if (t < jTemperature) {
            houdayerClusterMove(i);
        }
    }

    parallelTemperingUpdate();
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

void Lattice::shapeError() const {
    std::cout << "Wrong lattice type. Aborting...\n\n";
    exit(EXIT_FAILURE);
}

////////////////////////
// RectangularLattice //
////////////////////////

RectangularLattice::RectangularLattice(Hamiltonian h, ldouble t, ldouble dt,
                                       int n, char m, int r, int c)
    : Lattice(h, t, dt, n, m) {
    setType("rectangle");
    setRows(r);
    setCols(c);
    checkShape();

    if (getRows() == -1 || getRows() == -1) {
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
    int iRow = prop.locations[i][0];
    int jRow = prop.locations[j][0];
    return iRow - jRow;
}

int RectangularLattice::findYDisplacement(int i, int j) {
    int iCol = prop.locations[i][1];
    int jCol = prop.locations[j][1];
    return iCol - jCol;
}

ldouble RectangularLattice::findDistance(int i, int j) {
    return sqrt(pow(prop.xDisplacements[i][j], 2) +
                pow(prop.yDisplacements[i][j], 2));
}

void RectangularLattice::guessRowsCols() {
    if (prop.hamiltonian.getRows() != -1 && prop.hamiltonian.getCols() != -1) {
        setRows(prop.hamiltonian.getRows());
        setCols(prop.hamiltonian.getCols());
        return;
    }

    for (int guess = (int)sqrt(prop.numIndices); guess > 0; --guess) {
        if (prop.numIndices % guess == 0) {
            setRows(guess);
            setCols(prop.numIndices / guess);
            return;
        }
    }
}

///////////////////
// SquareLattice //
///////////////////

SquareLattice::SquareLattice(Hamiltonian h, ldouble t, ldouble dt, int n, char m,
                             int s)
    : Lattice(h, t, dt, n, m), RectangularLattice(h, t, dt, n, m, s, s) {
    setType("square");
    checkShape();

    if (getRows() == -1 || getRows() != getCols()) {
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
    int guess = prop.hamiltonian.getRows();

    if (guess == -1) {
        guess = (int)sqrt(prop.numIndices);

        if (guess * guess != (int)prop.numIndices) {
            std::cout << "Lattice is not a square (invalid number of indices)!";
            shapeError();
        }
    }

    setRows(guess);
    setCols(guess);
}

///////////////////////
// TriangularLattice //
///////////////////////

TriangularLattice::TriangularLattice(Hamiltonian h, ldouble t, ldouble dt, int n,
                                     char m, int r, int c)
    : Lattice(h, t, dt, n, m) {
    setType("triangle");
    checkShape();
    setRows(r);
    setCols(c);

    if (getRows() == -1 || getRows() == -1) {
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
    int iCol = prop.locations[i][0];
    int jCol = prop.locations[j][0];
    return iCol - jCol;
}

int TriangularLattice::findYDisplacement(int i, int j) {
    int iCol = prop.locations[i][1];
    int jCol = prop.locations[j][1];
    return iCol - jCol;
}

ldouble TriangularLattice::findDistance(int i, int j) {
    return sqrt(pow(prop.xDisplacements[i][j], 2) +
                pow(prop.yDisplacements[i][j], 2));
}

void TriangularLattice::guessRowsCols() {
    if (prop.hamiltonian.getRows() != -1 && prop.hamiltonian.getCols() != -1) {
        setRows(prop.hamiltonian.getRows());
        setCols(prop.hamiltonian.getCols());
        return;
    }

    for (int guess = (int)sqrt(prop.numIndices); guess > 0; --guess) {
        if (prop.numIndices % guess == 0) {
            setRows(guess);
            setCols(prop.numIndices / guess);
            return;
        }
    }
}

////////////////////////
// STriangularLattice //
////////////////////////

STriangularLattice::STriangularLattice(Hamiltonian h, ldouble t, ldouble dt,
                                       int n, char m, int s)
    : Lattice(h, t, dt, n, m), TriangularLattice(h, t, dt, n, m, s, s) {
    setType("square triangle");
    checkShape();

    if (getRows() == -1 || getRows() != getCols()) {
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
    int guess = prop.hamiltonian.getRows();

    if (guess == -1) {
        guess = (int)sqrt(prop.numIndices);

        if (guess * guess != (int)prop.numIndices) {
            std::cout << "Lattice is not a square triangle ";
            std::cout << "(invalid number of indices)!";
            shapeError();
        }
    }

    setRows(guess);
    setCols(guess);
}
