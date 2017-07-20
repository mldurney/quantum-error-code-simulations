#include "simulatedlattice.h"

using namespace ising;

SimulatedLattice::SimulatedLattice(Lattice *latt, const unsigned int updates,
                                   const unsigned int preupdates)
    : lattice(latt), updates(updates), preupdates(preupdates) {
    indLattice = numLattices++;
    temperatures.push_back(-1);
    magnetizations.push_back(-1);
    binderCumulants.push_back(-1);
    runPreupdates();
}

void SimulatedLattice::runPreupdates() {
    if (preupdates > 0) {
        for (unsigned int i = 0; i < preupdates; ++i) {
            lattice->updateLattice();
        }
    } else {
        while (preupdates < INITUPDATES) {
            lattice->updateLattice();
            ++preupdates;
        }

        bool continueUpdating = true;
        double preMagnetization = lattice->getMagnetism();
        unsigned int numUpdates = INITUPDATES * 2;
        unsigned int maxUpdates = 110000;
        double threshold = .05;

        while (continueUpdating) {
            while (preupdates < numUpdates) {
                lattice->updateLattice();
                ++preupdates;
            }

            double updatedMagnetization = lattice->getMagnetism();
            double change = abs((updatedMagnetization - preMagnetization) /
                                preMagnetization);

            if (change < threshold && numUpdates < maxUpdates) {
                continueUpdating = false;
            }

            numUpdates *= 2;
        }
    }
}

void SimulatedLattice::runLatticeSimulation() {
    unsigned int skip = 10;
    double runningMag = 0;
    double runningMag2 = 0;
    double runningMag4 = 0;

    for (unsigned int i = 0; i < updates; ++i) {
        for (unsigned int j = 0; j < skip; ++j) {
            lattice->updateLattice();
        }

        double magnetization = lattice->getMagnetism();
        runningMag += magnetization;
        runningMag2 += pow(magnetization, 2);
        runningMag4 += pow(magnetization, 4);
    }

    setAveMag(fabs(runningMag) / updates);
    setAveMag2(fabs(runningMag2) / updates);
    setAveMag4(fabs(runningMag4) / updates);

    double currT = lattice->getTemp();
    double currM = getAveMag();
    double currBC = getBinderCumulant();

    // std::lock_guard<std::mutex> guard(data_mutex);
    temperatures[indLattice] = currT;
    magnetizations[indLattice] = currM;
    binderCumulants[indLattice] = currBC;
}

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
