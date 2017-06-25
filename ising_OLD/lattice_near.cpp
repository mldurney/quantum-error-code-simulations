#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lattice_near.h"

Lattice::Lattice(int n, double t, char m) : Lattice(n, n, t, m) {}

Lattice::Lattice(int r, int c, double t, char m) :
        rows(r), cols(c), temp(t), mode(m) {
    board = new int*[rows];

    for (int i = 0; i < rows; i++) {
        board[i] = new int[cols];
    }

    srand(time(NULL));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            board[i][j] = (rand() % 2 == 1) ? 1 : -1;
        }
    }
}

void Lattice::updateLattice() {
    switch (mode) {
        case ALL: updateAll(); break;
        case RANDOM: updateRandom(); break;
    }
}

void Lattice::updateAll() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (findProbability(i, j) > (double) rand() / (double) RAND_MAX) {
                board[i][j] *= -1;
            }
        }
    }
}

void Lattice::updateRandom() {
    for (int i = 0; i < rows * cols; i++) {
        int x = rand() % rows;
        int y = rand() % cols;

        if (findProbability(x, y) > (double) rand() / (double) RAND_MAX) {
            board[x][y] *= -1;
        }
    }
}

void Lattice::switchMode(char mode) {
    switch (mode) {
        case ALL: mode = ALL; break;
        case RANDOM: mode = RANDOM; break;
        default: printf("\nINVALID MODE. Exiting...\n"); exit(EXIT_FAILURE);
    }
}

double Lattice::findProbability(int x, int y) {
    double energyInit = findIndexEnergy(x, y);
    board[x][y] *= -1;
    double energyFinal = findIndexEnergy(x, y);
    board[x][y] *= -1;

    if (energyFinal > energyInit) {
        return pow(E, (-1 / temp) * (energyFinal - energyInit));
    } else {
        return 1;
    }
}

double Lattice::findTotalEnergy() {
    double energy = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            energy += findIndexEnergy(i, j);
        }
    }

    return energy;
}

double Lattice::findIndexEnergy(int x, int y) {
    return -(board[x][y] *
            ( board[x][(y + 1) % cols]
            + board[x][(((y - 1) % cols) + cols) % cols]
            + board[(x + 1) % rows][y]
            + board[(((x - 1) % rows) + rows) % rows][y]));
}

double Lattice::findMagnetism() {
    double magnetism = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            magnetism += board[i][j];
        }
    }

    return magnetism / (rows * cols);
}

void Lattice::printLattice() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            (board[i][j] == 1) ? printf("+ ") : printf("- ");
        }

        printf("\n");
    }
}

Lattice::~Lattice() {
    for (int i = 0; i < rows; i++) {
        delete[] board[i];
    }

    delete[] board;
}
