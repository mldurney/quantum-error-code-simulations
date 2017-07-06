#ifndef ISINGHELPERS_H_
#define ISINGHELPERS_H_

#include <stdio.h>
#include "lattices.h"

const int MAX_FILENAME_SIZE = 255;

string receiveFilename(int argc, char* argv[]);
Hamiltonian readHamiltonian(ifstream& file, char& shape);
string getOutFilename(string inFilename, string oldDir, string newDir);
void writeOutput(string filename, vector<double> temp, vector<double> results);

#endif /* ISINGHELPERS_H_ */
