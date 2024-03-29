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

Lattice* ising::chooseLattice(char shape, const Hamiltonian& h, double t,
                              double dt, uint n, char m) {
    Lattice* lattice;

    switch (shape) {
        case RECTANGLE:
            lattice = new RectangularLattice(h, t, dt, n, m);
            break;
        case SQUARE:
            lattice = new SquareLattice(h, t, dt, n, m);
            break;
        case TRIANGLE:
            lattice = new TriangularLattice(h, t, dt, n, m);
            break;
        case STRIANGLE:
            lattice = new STriangularLattice(h, t, dt, n, m);
            break;
        default:
            lattice = new Lattice(h, t, dt, n, m);
    }

    return lattice;
}

std::string ising::getOutFilename(const std::string& inFilename,
                                  const std::string& newDir) {
    fs::path outFilename(inFilename);
    std::ostringstream latticeName;
    latticeName << outFilename.filename().string();

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
    latticeName << outFilename.filename().string();

    if (outFilename.string().find(oldDir) != std::string::npos) {
        outFilename.remove_filename();
    }

    outFilename.replace_filename(newDir + "/");
    outFilename /= latticeName.str();

    return outFilename.string();
}

void ising::writeOutput(const std::string& filename, const dmap& temperatures,
                        const dmap& results) {
    fs::path directory(filename);
    directory.remove_filename();
    fs::create_directory(directory);

    std::ofstream file(filename.c_str());
    file << "temperature,result\n";

    ivector indices;
    for (auto& t : temperatures) {
        indices.push_back(t.first);
    }

    for (auto& i : indices) {
        std::ostringstream row;
        row << temperatures.at(i) << "," << results.at(i) << "\n";
        file << row.str();
    }

    file.close();

    std::cout << "Recorded results by temperature in " << filename;
    std::cout << std::endl;
}
