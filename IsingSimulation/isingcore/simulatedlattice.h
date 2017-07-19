#ifndef SIMULATEDLATTICE_H_
#define SIMULATEDLATTICE_H_

#include <mutex>
#include "lattices.h"

namespace ising {
const int INITUPDATES = 100;

class SimulatedLattice {
   public:
    SimulatedLattice(Lattice *lattice, const unsigned int updates,
                     const unsigned int preupdates = 0);
    ~SimulatedLattice() {}
    void runPreupdates();
    void runLatticeSimulation();
    unsigned int getIndLattice() const { return indLattice; }
    Lattice *getLattice() const { return lattice; }
    unsigned int getUpdates() const { return updates; }
    double getAveMag() const { return aveMag; }
    double getAveMag2() const { return aveMag2; }
    double getAveMag4() const { return aveMag4; }
    double getBinderCumulant();

    static int getNumLattices() { return numLattices; }
    static std::vector<double> getTemperatures() { return temperatures; }
    static std::vector<double> getMagnetizations() { return magnetizations; }
    static std::vector<double> getBinderCumulants() { return binderCumulants; }

   protected:
    void setAveMag(double mag);
    void setAveMag2(double mag2);
    void setAveMag4(double mag4);

   private:
    unsigned int indLattice;
    Lattice *lattice;
    unsigned int updates;
    unsigned int preupdates;
    double aveMag;
    double aveMag2;
    double aveMag4;

    static int numLattices;
    static std::vector<double> temperatures;
    static std::vector<double> magnetizations;
    static std::vector<double> binderCumulants;
    static std::mutex data_mutex;
};
}

#endif /* SIMULATEDLATTICE_H_ */