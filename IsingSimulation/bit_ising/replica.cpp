#include "replica.h"

using namespace ising;

void Replica::initSpins() {
    for (auto& i : spins) {
        i = (MWC() % 2 == 0) ? 1 : -1;
    }
}

void Replica::metropolisUpdate() {
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

void Replica::updateAll() {
    for (int i = 0; i < numIndices; ++i) {
        if (findProbability(i) > randFloatCO()) {
            spins[i] *= -1;
        }
    }
}

void Replica::updatePseudo() {
    for (int i = 0; i < numIndices; ++i) {
        int j = MWC() % numIndices;
        std::swap(randomizedIndices[i], randomizedIndices[j]);
    }

    for (int i = 0; i < numIndices; ++i) {
        int index = randomizedIndices[i];
        if (findProbability(index) > randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

void Replica::updateRandom() {
    int index;

    for (int i = 0; i < numIndices; ++i) {
        index = indices[MWC() % numIndices];

        if (findProbability(index) > randFloatCO()) {
            spins[index] *= -1;
        }
    }
}

double Replica::findProbability(int index) {
    int initEnergy = findIndexEnergy(index);

    if (initEnergy < 0) {
        return pow(E, (-1 / getTemp()) * (-2 * initEnergy));
    } else {
        return 1;
    }
}

int Replica::findTotalEnergy() {
    int energy = 0;

    for (int i = 0; i < numIndices; ++i) {
        energy += findIndexEnergy(i);
    }

    return energy;
}

int Replica::findIndexEnergy(int index) {
    int energy = 0;

    ivector2::const_iterator it1;
    ivector::const_iterator it2;

    auto end = indInteractions[index].end();
    for (it1 = indInteractions[index].begin(); it1 != end; ++it1) {
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

    for (int i = 0; i < numIndices; ++i) {
        magnetism += spins[i];
    }

    return magnetism / numIndices;
}

void Replica::flipSpins() {
    for (auto& i : spins) {
        i *= -1;
    }
}