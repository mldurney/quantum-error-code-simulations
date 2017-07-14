#include "isinghelpers.h"
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
static const std::string SLASH = "\\";
#else
static const std::string SLASH = "/";
#endif

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

Lattice* ising::chooseLattice(const char shape, const Hamiltonian& hamiltonian, const double temp, 
	const char mode) {
	Lattice* lattice;

	if (hamiltonian.getIsFast()) {
		switch (shape) {
		case RECTANGLE:
			lattice = new RectangularLatticeFast(hamiltonian, temp, mode);
			break;
		case SQUARE:
			lattice = new SquareLatticeFast(hamiltonian, temp, mode);
			break;
		case TRIANGLE:
			lattice = new TriangularLatticeFast(hamiltonian, temp, mode);
			break;
		case STRIANGLE:
			lattice = new STriangularLatticeFast(hamiltonian, temp, mode);
		default:
			lattice = new LatticeFast(hamiltonian, temp, mode);
		}
	}

	else {
		switch (shape) {
		case RECTANGLE:
			lattice = new RectangularLattice(hamiltonian, temp, mode);
			break;
		case SQUARE:
			lattice = new SquareLattice(hamiltonian, temp, mode);
			break;
		case TRIANGLE:
			lattice = new TriangularLattice(hamiltonian, temp, mode);
			break;
		case STRIANGLE:
			lattice = new STriangularLattice(hamiltonian, temp, mode);
		default:
			lattice = new Lattice(hamiltonian, temp, mode);
		}
	}

	return lattice;
}

std::string ising::getOutFilename(const std::string& inFilename, const std::string& oldDir,
                                  const std::string& newDir) {
    std::string outFilename;

    if (inFilename.find(SLASH) == inFilename.rfind(SLASH)) {
        std::string objectName = inFilename.substr(0, inFilename.rfind('.'));

		std::ostringstream outFilenameBuilder;
		outFilenameBuilder << newDir << SLASH << objectName << ".csv";
		outFilename = outFilenameBuilder.str();
    } else {
        outFilename = inFilename;
        outFilename.replace(outFilename.find(oldDir), oldDir.length(), newDir);
    }

    return outFilename;
}

void ising::writeOutput(const std::string& filename, const std::vector<double>& temp,
                        const std::vector<double>& results) {
    bool isNewFile = (std::ifstream(filename)) ? false : true;
    std::ofstream file(filename.c_str(),
                       std::ofstream::out | std::ofstream::app);

    std::vector<double>::const_iterator it;

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
