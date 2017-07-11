#include "isinghelpers.h"

using namespace ising;

Hamiltonian ising::readHamiltonian(std::ifstream &file, char &shape) {
    shape = '\0';
    int rows = -1;
    int cols = -1;

    if (isalpha(file.peek())) {
        std::string line;
        std::string field;
        getline(file, line);

        std::stringstream stream(line);
        std::stringstream fields[3];

        fields[0] << "\0";
        fields[1] << "-1";
        fields[2] << "-1";

        for (int i = 0; i < 3 && getline(stream, field, ','); ++i) {
            fields[i].str(field);
        }

        fields[0] >> shape;
        fields[1] >> rows;
        fields[2] >> cols;
    }

    return Hamiltonian(importHamiltonianVector(file), shape, rows, cols);
}

std::string ising::getOutFilename(std::string inFilename, std::string oldDir,
                                  std::string newDir) {
    std::string outFilename;

    if (inFilename.find('/') == inFilename.rfind('/')) {
        std::string objectName = inFilename.substr(0, inFilename.rfind('.'));

        char buffer[MAX_FILENAME_SIZE + 1];
        sprintf(buffer, "%s/%s.csv", newDir.c_str(), objectName.c_str());

        outFilename = buffer;
    } else {
        outFilename = inFilename;
        outFilename.replace(outFilename.find(oldDir), oldDir.length(), newDir);
    }

    return outFilename;
}

void ising::writeOutput(std::string filename, std::vector<double> temp,
                        std::vector<double> results) {
    bool isNewFile = (std::ifstream(filename)) ? false : true;
    std::ofstream file(filename.c_str(),
                       std::ofstream::out | std::ofstream::app);

    std::vector<double>::iterator it;

    if (isNewFile) {
        for (it = temp.begin(); it != temp.end(); ++it) {
            file << *it;

            if (it + 1 != temp.end()) {
                file << ",";
            }
        }
    }

    file << std::endl;

    for (it = results.begin(); it != results.end(); ++it) {
        file << *it;

        if (it + 1 != results.end()) {
            file << ",";
        }
    }

    file.close();

    std::cout << "Recorded results by temperature in " << filename;
    std::cout << std::endl;
}
