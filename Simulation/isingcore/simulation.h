#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "isinghelpers.h"
#include "simulatedlattice.h"
#include "threadpool.h"
#include "threadpoolhelpers.h"

typedef std::unique_ptr<ising::SimulatedLattice> latticeptr;
typedef std::map<int, latticeptr> latticemap;

namespace ising {
const ldouble MIN_PERCENTILE = .5;
const ldouble MAX_PERCENTILE = .9;

class Simulation {
   public:
    Simulation(const std::string &filename, ldouble t, ldouble dt, uint n,
               uint updates = 0, uint trials = 1, char mode = 'p');
    Simulation(const std::string &filename, ldouble t, ldouble dt, uint n,
               uint updates, uint preupdates, uint trials, char mode);
    ~Simulation() = default;
    void runSimulation();

    const std::string getFilename() const { return inFilename; };
    ldouble getMinT() const { return minT; }
    ldouble getDT() const { return dT; }
    uint getNumT() const { return numT; }
    uint getUpdates() const { return updates; }
    uint getPreupdates() const { return preupdates; }
    uint getTrials() const { return trials; }
    char getMode() const { return mode; }

    ldouble getAvgMag(uint n) { return findAverage(avgMag[n]); }
    ldouble getAvgMag2(uint n) { return findAverage(avgMag2[n]); }
    ldouble getAvgMag4(uint n) { return findAverage(avgMag4[n]); }
    cdouble getChi0(uint n) { return findAverage(chi0[n]); }
    cdouble getChiq(uint n) { return findAverage(chiq[n]); }

    ldouble getBinderCumulant(uint n);
    cdouble getCorrelationFunction(uint n);

    const dmap &getTemperatures() const { return temperatures; }
    const dmap &getMagnetizations() const { return magnetizations; }
    const dmap &getBinderCumulants() const { return binderCumulants; }
    const cdmap &getCorrelationFunctions() const {
        return correlationFunctions;
    }

    dmap getRealCorrelationFunctions();

   protected:
    void addAvgMag(uint n, ldouble mag);
    void addAvgMag2(uint n, ldouble mag2);
    void addAvgMag4(uint n, ldouble mag4);
    void addChi0(uint n, cdouble chi) { chi0[n].push_back(chi); }
    void addChiq(uint n, cdouble chi) { chiq[n].push_back(chi); }
    ldouble findAverage(dvector &v);
    cdouble findAverage(cdvector &v);
    ldouble findAverageNoOutliers(dvector &v,
                                  ldouble percentile1 = MIN_PERCENTILE,
                                  ldouble percentile2 = MAX_PERCENTILE);
    cdouble findAverageNoOutliers(cdvector &v,
                                  ldouble percentile1 = MIN_PERCENTILE,
                                  ldouble percentile2 = MAX_PERCENTILE);

   private:
    void checkInputFile();
    void initTempDirectory();
    void loadTempData();
    uint getLeadingInt(const fs::path &filename);
    uint temperatureToIndex(ldouble t);
    void initLattices();
    void initPersonalLattice();
    void addLattice(uint trial);
    void runTrials();
    void runTrial(uint trial);
    void initRunTrials();
    void initRunTrial(uint trial);

    const std::string &inFilename;
    ldouble minT;
    ldouble dT;
    uint numT;
    uint updates;
    uint preupdates;
    uint trials;
    char mode;

    latticemap lattices;
    latticeptr personalLattice;
    fs::path tempDirectory;
    ivector remainingTrials;

    dvectormap avgMag;
    dvectormap avgMag2;
    dvectormap avgMag4;
    cdvectormap chi0;
    cdvectormap chiq;

    dmap temperatures;
    dmap magnetizations;
    dmap binderCumulants;
    cdmap correlationFunctions;

    static std::mutex file_mutex;
    static std::mutex trial_mutex;
    static std::mutex avgMag_mutex;
    static std::mutex avgMag2_mutex;
    static std::mutex avgMag4_mutex;
    static std::mutex chi0_mutex;
    static std::mutex chiq_mutex;
};
}

#endif /* SIMULATION_H_ */