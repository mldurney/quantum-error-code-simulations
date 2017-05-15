#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lattice.h"

Lattice::Lattice(int n, float t) : Lattice(n, n, t) {}

Lattice::Lattice(int r, int c, float t) : rows(r), cols(c), temp(t) {
    board = new int*[rows];
    probabilities = new float*[rows];

    for (int i = 0; i < rows; i++) {
        board[i] = new int[cols];
        probabilities[i] = new float[cols];
    }

    srand(time(NULL));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            board[i][j] = (rand() % 2 == 1) ? 1 : -1;
            probabilities[i][j] = 0;
        }
    }
}

void Lattice::updateLattice() {
    findProbabilities();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (probabilities[i][j] > (float) rand() / (float) RAND_MAX) {
                board[i][j] *= -1;
                //printf("%f", probabilities[i][j]);
            }
        }
    }
}

void Lattice::findProbabilities() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float energy1 = findTotalEnergy();
            board[i][j] *= -1;
            float energy2 = findTotalEnergy();
            board[i][j] *= -1;

            probabilities[i][j] = pow(M_E, (-1 / temp) * (energy2 - energy1));
        }
    }
}

float Lattice::findTotalEnergy() {
    float energy = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            energy += findIndexEnergy(i, j);
        }
    }

    return energy;
}

float Lattice::findIndexEnergy(int x, int y) {
    float energy = 0;
    for (int i = 0; i <= MAX_DISTANCE; i++) {
        for (int j = 0; distance(i, j) <= MAX_DISTANCE; j++) {
            if (distance(i, j) < 1) {
                continue;
            } else {
                if (i == 0) {
                    energy -= board[x][y] / pow(distance(i, j), 2)
                        * (board[x][(y + j) % cols]
                        +  board[x][(((y - j) % cols) + cols) % cols]);
                } else if (j == 0) {
                    energy -= board[x][y] / pow(distance(i, j), 2)
                        * (board[(x + i) % rows][y]
                        +  board[(((x - i) % rows) + rows) % rows][y]);
                } else {
                    energy -= board[x][y] / pow(distance(i, j), 2)
                        * (board[(x + i) % rows][(y + j) % cols]
                        +  board[(x + i) % rows][(((y - j) % cols) + cols) % cols]
                        +  board[(((x - i) % rows) + rows) % rows][(y + j) % cols]
                        +  board[(((x - i) % rows) + rows) % rows][(((y - j) % cols) + cols) % cols]);
                }
            }
        }
    }

    return energy;
}

float Lattice::findMagnetism() {
    float magnetism = 0;
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
        delete[] probabilities[i];
    }

    delete[] board;
    delete[] probabilities;
}
