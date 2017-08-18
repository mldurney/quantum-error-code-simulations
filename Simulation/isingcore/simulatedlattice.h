#ifndef SIMULATEDLATTICE_H_
#define SIMULATEDLATTICE_H_

#include <experimental/filesystem>
#include <mutex>
#include <numeric>
#include "lattices.h"

namespace fs = std::experimental::filesystem;

namespace ising {
const uint PREUPDATES = 500;
const uint BASEUPDATES = 5;
const uint SKIP = 10;
enum { MAG = 'm', CHI0 = 'x', ENERGY = 'e' };

class SimulatedLattice {
   public:
    SimulatedLattice(Lattice* lattice, const std::string& filename, uint index,
                     uint updates, uint preupdates = PREUPDATES,
                     bool suppress = false);
    ~SimulatedLattice() { delete lattice; }
    void runLatticeSimulation();
    Lattice* getLattice() const { return lattice; }
    uint getIndLattice() const { return indLattice; }
    uint getUpdates() const { return updates; }
    uint getPreupdates() const { return preupdates; }
    uint getSize() const { return getLattice()->getSize(); }
    ldouble getQ() const { return q; }
    char getStabilityMode() const { return stabilityMode; }

    const dmap& getTemperatures() const { return temperatures; }
    const dmap& getAvgMag() const { return avgMag; }
    const dmap& getAvgMag2() const { return avgMag2; }
    const dmap& getAvgMag4() const { return avgMag4; }
    const cdmap& getChi0() const { return chi0; }
    const cdmap& getChiq() const { return chiq; }

   protected:
    void setQ(ldouble qNew) { q = qNew; }
    void setStabilityMode(char m);

    void addAvgMag(uint index, ldouble mag);
    void addAvgMag2(uint index, ldouble mag2);
    void addAvgMag4(uint index, ldouble mag4);
    void addChi0(uint index, cdouble chi) { chi0[index] = chi; }
    void addChiq(uint index, cdouble chi) { chiq[index] = chi; }

   private:
    Lattice* lattice;
    uint indLattice;
    uint updates;
    uint preupdates;
    ldouble q;
    char stabilityMode = MAG;

    dmap temperatures;
    dmap avgMag;
    dmap avgMag2;
    dmap avgMag4;
    cdmap chi0;
    cdmap chiq;

    fs::path tempDirectory;
    fs::path tempFile;
    static std::mutex file_mutex;

    void initTempFile(const std::string& filename);
    void updateTempFile();
    void runPreupdates();
    void runUpdates();
    void runUpdatesStable();
    uint reachStability();
    uint reachStabilityMag();
    uint reachStabilityChi0();
    uint reachStabilityEnergy();
};
}

#endif /* SIMULATEDLATTICE_H_ */