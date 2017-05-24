/*
 * TODO:
 *      Add change tracking
 *      Comment to glory
 *      Make variables more accessible (at least read after start)
 *          Maybe even change...?
 *      Enable percentage tracking
 *      Add more parameters
 *          Like what?
 *      Add safety nets and what not
 *      Optimize in general... Very slow
 *      Explore possibilities of distance
 *      Create list of boards to see progression over time
 *      Enable output of results
 *      Error checking
 *
 *      Change probabilities:
 *          ef > ei : e^(-1/T)(ef-ei)
 *          ef < ei : 1
 *      Optimize for distance 1
 *      Include option to compare changes
 *      Pick random points L*L times
 *
 *      Desired behavior:
 *          Plot of average value of absolute magnetization vs temp
 *
 - Focus on nearest-neighbour. We can discuss next nearest neighbours at a later time
- Change the probability function based off energy difference
- Optimize the code for calculating energy differences at a given time
- Obtain a plot of the average (absolute value) of the magnetization as a function of temperature
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lattice.h"

Lattice::Lattice(int n, float t, char m) : Lattice(n, n, t, m) {}

Lattice::Lattice(int r, int c, float t, char m) :
        rows(r), cols(c), temp(t), mode(m) {
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
    switch (mode) {
        case ALL: updateAll(); break;
        case RANDOM: updateRandom(); break;
    }
}

void Lattice::updateAll() {
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

void Lattice::updateRandom() {  // incompatible probabilities
    findProbabilities();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int n1 = rand() % i;
            int n2 = rand() % j;
            if (probabilities[n1][n2] > (float) rand() / (float) RAND_MAX) {
                board[i][j] *= -1;
            }
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

void Lattice::findProbabilities() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float energy1 = findTotalEnergy();
            board[i][j] *= -1;
            float energy2 = findTotalEnergy();
            board[i][j] *= -1;

            probabilities[i][j] = pow(E, (-1 / temp) * (energy2 - energy1));
        }
    }
}

float Lattice::findTotalEnergy() {
    float energy = 0;

    if (MAX_DISTANCE == 1) {
        energy -= board[x][y] *
                ( board[x][(y + 1) % cols]
                + board[x][(((y - 1) % cols) + cols) % cols]
                + board[(x + 1) % rows][y]
                + board[(((x - 1) % rows) + rows) % rows][y]);
    } else {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                energy += findIndexEnergy(i, j);
            }
        }
    }

    return energy;
}

float Lattice::findIndexEnergy(int x, int y) {
    float energy = 0;

    if (MAX_DISTANCE == 1) {
        energy -= board[x][y] *
                ( board[x][(y + 1) % cols]
                + board[x][(((y - 1) % cols) + cols) % cols]
                + board[(x + 1) % rows][y]
                + board[(((x - 1) % rows) + rows) % rows][y]);
    } else {
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
