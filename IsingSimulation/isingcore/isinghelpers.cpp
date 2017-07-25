#include "isinghelpers.h"
#if defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32) && !defined(__CYGWIN__)
static const std::string SLASH = "\\";
#else
static const std::string SLASH = "/";
#endif

using namespace ising;

Hamiltonian ising::readHamiltonian(std::ifstream& file, char& shape) {
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

Lattice* ising::chooseLattice(const char shape, const Hamiltonian& hamiltonian,
                              const double temp, const char mode, const bool init) {
    Lattice* lattice;

    if (hamiltonian.getIsFast()) {
        switch (shape) {
            case RECTANGLE:
                lattice = new RectangularLatticeFast(hamiltonian, temp, mode, init);
                break;
            case SQUARE:
                lattice = new SquareLatticeFast(hamiltonian, temp, mode, init);
                break;
            case TRIANGLE:
                lattice = new TriangularLatticeFast(hamiltonian, temp, mode, init);
                break;
            case STRIANGLE:
                lattice = new STriangularLatticeFast(hamiltonian, temp, mode, init);
                break;
            default:
                lattice = new LatticeFast(hamiltonian, temp, mode, init);
        }
    }

    else {
        switch (shape) {
            case RECTANGLE:
                lattice = new RectangularLattice(hamiltonian, temp, mode, init);
                break;
            case SQUARE:
                lattice = new SquareLattice(hamiltonian, temp, mode, init);
                break;
            case TRIANGLE:
                lattice = new TriangularLattice(hamiltonian, temp, mode, init);
                break;
            case STRIANGLE:
                lattice = new STriangularLattice(hamiltonian, temp, mode, init);
                break;
            default:
                lattice = new Lattice(hamiltonian, temp, mode, init);
        }
    }

    return lattice;
}

std::string ising::getOutFilename(const std::string& inFilename,
								  const std::string& newDir) {
	fs::path outFilename(inFilename);
	std::ostringstream latticeName;
	latticeName << outFilename.filename();

	outFilename.remove_filename();
	outFilename.replace_filename(newDir + "/");
	outFilename /= latticeName.str();

	return outFilename.string();
}

std::string ising::getOutFilename(const std::string& inFilename,
								  const std::string& oldDir,
								  const std::string& newDir) {
	fs::path outFilename(inFilename);
	std::ostringstream latticeName;
	latticeName << outFilename.filename();

	if (outFilename.string().find(oldDir) != std::string::npos) {
		outFilename.remove_filename();
	}

	outFilename.replace_filename(newDir + "/");
	outFilename /= latticeName.str();

	return outFilename.string();
}

void ising::writeOutput(const std::string& filename,
                        const dvector& temp,
                        const dvector& results) {
    bool isNewFile = (std::ifstream(filename)) ? false : true;
    std::ofstream file(filename.c_str(),
                       std::ofstream::out | std::ofstream::app);

    dvector::const_iterator it;

    if (isNewFile) {
        std::ostringstream header;
        for (it = temp.begin(); it != temp.end(); ++it) {
            header << *it;

            if (it + 1 != temp.end()) {
                header << ",";
            }
        }
        header << "\n";
        file << header.str();
    }

    std::ostringstream row;
    for (it = results.begin(); it != results.end(); ++it) {
        row << *it;

        if (it + 1 != results.end()) {
            row << ",";
        }
    }
    row << "\n";
    file << row.str();

    file.close();

    std::cout << "Recorded results by temperature in " << filename;
    std::cout << std::endl;
}
