#include "spins.h"

void Lattice::updateAll() {
    for (int n = spins.size() - 1; n >= 0; --n) {
        for (int i = 0; i < numIndices; ++i) {
            if (findProbability(i, n) > randFloatCO()) {
                spins[n][i] *= -1;
            }
        }
    }
}

void Lattice::updatePseudo() {
    for (int i = 0; i < numIndices; ++i) {
        int j = MWC() % numIndices;
        std::swap(randomizedIndices[i], randomizedIndices[j]);
    }

    for (int n = spins.size() - 1; n >= 0; --n) {
        for (int i = 0; i < numIndices; ++i) {
            int index = randomizedIndices[i];
            if (findProbability(index, n) > randFloatCO()) {
                spins[n][index] *= -1;
            }
        }
    }
}

void Lattice::updateRandom() {
    int index;

    for (int n = spins.size() - 1; n >= 0; --n) {
        for (int i = 0; i < numIndices; ++i) {
            index = indices[MWC() % numIndices];

            if (findProbability(index, n) > randFloatCO()) {
                spins[n][index] *= -1;
            }
        }
    }
}