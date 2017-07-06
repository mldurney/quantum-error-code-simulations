#include "simulatedlattice.h"

SimulatedLattice::SimulatedLattice(Lattice* latt, int preupdates) :
        lattice(latt)
{
    for (int i = 0; i < preupdates; ++i)
    {
        latt->updateLattice();
    }
}

SimulatedLattice::~SimulatedLattice()
{
    delete lattice;
}

void SimulatedLattice::setAveMag(double mag)
{
    if (mag >= 0 && mag <= 1)
    {
        aveMag = mag;
    }
    else
    {
        cout << "\nInvalid average magnetization!";
        cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::setAveMag2(double mag2)
{
    if (mag2 >= 0 && mag2 <= 1)
    {
        aveMag2 = mag2;
    }
    else
    {
        cout << "\nInvalid average squared magnetization!";
        cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::setAveMag4(double mag4)
{
    if (mag4 >= 0 && mag4 <= 1)
    {
        aveMag4 = mag4;
    }
    else
    {
        cout << "\nInvalid average fourth-power magnetization!";
        cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

double SimulatedLattice::getBinderCumulant()
{
    return 1 - getAveMag4() / (3 * pow(getAveMag2(), 2));
}
