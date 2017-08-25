#include "replica.h"

using namespace ising;

Replica::Replica(const LatticeProperties& properties, uint n)
    : prop(properties), replicaIndex(n) {
    temperature = prop.minT + replicaIndex * prop.dT;
    randomizedIndices = prop.indices;
    spins.resize(prop.numIndices);
    initSpins();
}

void Replica::initSpins() {
    for (auto& s : spins) {
        s = (gen.MWC() % 2 == 0) ? 1 : -1;
    }
}

void Replica::setTemperature(double t) {
    if (t < prop.minT || t > prop.minT + prop.numT * prop.dT) {
        std::cout << "TEMPERATURE OUTSIDE VALID RANGE! Exiting...\n\n";
        exit(EXIT_FAILURE);
    }

    temperature = t;
}

void Replica::update() {
    switch (prop.mode) {
        case ALL:
            updateAll();
            break;
        case PSEUDO:
            updatePseudo();
            break;
        case RANDOM:
            updateRandom();
            break;
        default:
            std::cout << "INVALID MODE! Exiting...\n\n";
            exit(EXIT_FAILURE);
    }
}

void Replica::updateAll() {
    for (uint i = 0; i < prop.numIndices; ++i) {
        if (findProbability(i) > gen.randFloatCO()) {
            spins[i] *= -1;
        }
    }
}

void Replica::updatePseudo() {
    for (uint i = 0; i < prop.numIndices; ++i) {
        int j = gen.MWC() % prop.numIndices;
        std::swap(randomizedIndices[i], randomizedIndices[j]);
    }

    for (uint i = 0; i < prop.numIndices; ++i) {
        int index = randomizedIndices[i];
        if (findProbability(index) > gen.randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

void Replica::updateRandom() {
    int index;

    for (uint i = 0; i < prop.numIndices; ++i) {
        index = prop.indices[gen.MWC() % prop.numIndices];

        if (findProbability(index) > gen.randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

double Replica::findProbability(int index) {
    int initEnergy = findIndexEnergy(index);

    if (initEnergy < 0) {
        return pow(E, (-1 / getTemperature()) * (-2 * initEnergy));
    } else {
        return 1;
    }
}

int Replica::findTotalEnergy() {
    int energy = 0;

    for (uint i = 0; i < prop.numIndices; ++i) {
        energy += findIndexEnergy(i);
    }

    return energy;
}

int Replica::findIndexEnergy(int index) {
    int energy = 0;

    ivector2::const_iterator it1;
    ivector::const_iterator it2;

    auto end = prop.indInteractions[index].end();
    for (it1 = prop.indInteractions[index].begin(); it1 != end; ++it1) {
        it2 = it1->begin();
        int couplingEnergy = *it2;

        for (++it2; it2 != it1->end(); ++it2) {
            couplingEnergy *= spins[*it2];
        }

        energy -= couplingEnergy;
    }

    return spins[index] * energy;
}

double Replica::findMagnetization() {
    double magnetism = 0;

    for (uint i = 0; i < prop.numIndices; ++i) {
        magnetism += spins[i];
    }

    return magnetism / prop.numIndices;
}

void Replica::flipSpins() {
    for (auto& s : spins) {
        s *= -1;
    }
}

void Replica::print() const {
    uint side = prop.cols;
    if (prop.cols == -1) {
        side = (int)sqrt(prop.numIndices);
    }

    for (uint i = 0; i < prop.numIndices;) {
        for (uint col = 0; col < side && i < prop.numIndices; ++col, ++i) {
            (spins.at(i) == 1) ? std::cout << "+ " : std::cout << "- ";
        }

        std::cout << "\n";
    }
}