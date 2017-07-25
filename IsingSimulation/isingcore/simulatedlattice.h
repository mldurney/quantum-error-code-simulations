#ifndef SIMULATEDLATTICE_H_
#define SIMULATEDLATTICE_H_

#include <complex>
#include <experimental/filesystem>
#include <mutex>
#include "lattices.h"

namespace fs = std::experimental::filesystem;
typedef std::complex<double> cdouble;
typedef std::vector<cdouble> cvector;

namespace ising {
const int INITUPDATES = 100;
const double MIN_PERCENTILE = .4;
const double MAX_PERCENTILE = 1;

#if defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32) && !defined(__CYGWIN__)
static const std::string SLASH = "\\";
#else
static const std::string SLASH = "/";
#endif

class SimulatedLattice {
   public:
    SimulatedLattice(Lattice *lattice, const std::string &filename,
                     const unsigned int updates,
                     const unsigned int preupdates = 0,
                     const unsigned int trials = 1);
    ~SimulatedLattice() {}
    void runLatticeSimulation();
    unsigned int getIndLattice() const { return indLattice; }
    Lattice *getLattice() const { return lattice; }
    unsigned int getUpdates() const { return updates; }
    double getAvgMag() { return findAverage(avgMag); }
    double getAvgMag2() { return findAverage(avgMag2); }
    double getAvgMag4() { return findAverage(avgMag4); }
    cdouble getChi0() { return findAverage(chi0); }
    cdouble getChiq() { return findAverage(chiq); }
    double getBinderCumulant();
    cdouble getCorrelationLength();

    static int getNumLattices() { return numLattices; }
    static dvector getTemperatures() { return temperatures; }
    static dvector getMagnetizations() { return magnetizations; }
    static dvector getBinderCumulants() { return binderCumulants; }
    static cvector getCorrelationLengths() { return correlationLengths; }
    static dvector getRealCorrelationLengths();

   protected:
    void addAvgMag(double mag);
    void addAvgMag2(double mag2);
    void addAvgMag4(double mag4);
    void addChi0(cdouble chi) { chi0.push_back(chi); }
    void addChiq(cdouble chi) { chiq.push_back(chi); }
    double findAverage(dvector &v, double percentile1 = MIN_PERCENTILE,
                       double percentile2 = MAX_PERCENTILE);
    cdouble findAverage(cvector &v, double percentile1 = MIN_PERCENTILE,
                        double percentile2 = MAX_PERCENTILE);

   private:
    unsigned int indLattice;
    Lattice *lattice;
    fs::path tempDirectory;
    fs::path tempData;
    unsigned int updates;
    unsigned int preupdates;
    unsigned int trials;
    double q;
    dvector avgMag;
    dvector avgMag2;
    dvector avgMag4;
    cvector chi0;
    cvector chiq;

    static int numLattices;
    static dvector temperatures;
    static dvector magnetizations;
    static dvector binderCumulants;
    static cvector correlationLengths;
    static std::mutex data_mutex;
    static std::mutex file_mutex;

    void initTempFiles(const std::string &filename);
    void updateTempFile();
    void runPreupdates();
    void runTrials();
};
}

#endif /* SIMULATEDLATTICE_H_ */