#include "simulatedlattice.h"

using namespace ising;

SimulatedLattice::SimulatedLattice(Lattice *latt, int preupdates)
    : lattice(latt) {
    for (int i = 0; i < preupdates; ++i) {
        latt->updateLattice();
    }
}

SimulatedLattice::~SimulatedLattice() { delete lattice; }

void SimulatedLattice::setAveMag(double mag) {
    if (mag >= 0 && mag <= 1) {
        aveMag = mag;
    } else {
        std::cout << "\nInvalid average magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::setAveMag2(double mag2) {
    if (mag2 >= 0 && mag2 <= 1) {
        aveMag2 = mag2;
    } else {
        std::cout << "\nInvalid average squared magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::setAveMag4(double mag4) {
    if (mag4 >= 0 && mag4 <= 1) {
        aveMag4 = mag4;
    } else {
        std::cout << "\nInvalid average fourth-power magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

double SimulatedLattice::getBinderCumulant() {
    return 1 - getAveMag4() / (3 * pow(getAveMag2(), 2));
}
