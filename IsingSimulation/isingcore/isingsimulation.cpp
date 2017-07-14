#include "isingsimulation.h"

using namespace ising;

int main(int argc, char *argv[]) {
    std::string inFilename;
    double t, dt;
    int n, updates;
    char mode;

    receiveInputCMD(argc, argv, inFilename, t, dt, n, updates, mode);
    manageSimulations(inFilename, t, dt, n, updates, mode);
}

void ising::manageSimulations(const std::string& inFilename, const double t, const double dt,
                              const int n, const int updates, const char mode) {
    std::ifstream inFile(inFilename);

    if (!inFile) {
        std::cout << "Invalid file name. " << inFilename
                  << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    char shape;
    Hamiltonian h = readHamiltonian(inFile, shape);

    std::vector<double> temperatures;
    std::vector<double> magnetizations(n, 0);
    std::vector<double> binderCumulants(n, 0);

    int remainingTrials = n;
    int trial = 0;

    Lattice *lattice;
    while (remainingTrials > 0) {
        std::vector<std::thread> threads;
        int numThreads = getNumThreads(remainingTrials);
		// Thread pools, queue instead of vector
		// Mutex to protect
		// Use lock_guard

        for (int i = 0; i < numThreads; ++i, ++trial, --remainingTrials) {
            double currT = t + trial * dt;
            temperatures.push_back(currT);
			lattice = chooseLattice(shape, h, currT, mode);

            threads.emplace_back(runLatticeSimulation, lattice, updates,
                                 std::ref(magnetizations[trial]),
                                 std::ref(binderCumulants[trial]));
        }

        for (auto &th : threads) {
            th.join();
        }
    }

    std::string oldDir = "hamiltonians";
    std::string outMag = getOutFilename(inFilename, oldDir, "magnetizations");
    std::string outBinder =
        getOutFilename(inFilename, oldDir, "binder_cumulants");

    writeOutput(outMag, temperatures, magnetizations);
    writeOutput(outBinder, temperatures, binderCumulants);
}

void ising::receiveInputCMD(int argc, char *argv[], std::string &filename,
                            double &t, double &dt, int &n, int &updates,
                            char &mode) {
    if (argc == 1) {
        std::cout << "Enter Hamiltonian input file: ";
        std::cin >> filename;
        std::cout << "Enter minimum temperature (K): ";
        std::cin >> t;
        std::cout << "Enter change in temperature between lattices (K): ";
        std::cin >> dt;
        std::cout << "Enter total number of lattices: ";
        std::cin >> n;
        std::cout << "Enter number of updates per lattice: ";
        std::cin >> updates;
        std::cout
            << "Enter update mode (a - all, p - pseudorandom, r - random): ";
        std::cin >> mode;
    } else if (argc == 7) {
        filename = argv[1];
        t = atof(argv[2]);
        dt = atof(argv[3]);
        n = atoi(argv[4]);
        updates = atoi(argv[5]);
        mode = (argv[6])[0];
    } else {
        std::cout << "Usage: " << argv[0] << " filename(std::string) ";
        std::cout << "min_temperature(float) change_temperature(float) ";
        std::cout << "num_lattices(int) updates(int) mode(char)\n\n";
        exit(EXIT_FAILURE);
    }
}

int ising::getNumThreads(const int remaining) {
    int numThreads = getMaxThreads();

    if (remaining < numThreads) {
        numThreads = remaining;
    }

    return numThreads;
}

int ising::getMaxThreads() {
    int maxThreads = std::thread::hardware_concurrency();

    if (maxThreads != 0) {
        return maxThreads;
    } else {
        return 4;
    }
}

void ising::runLatticeSimulation(Lattice *latt, const int updates,
                                 double &indMagnetization,
                                 double &indBinderCumulant) {
    SimulatedLattice simulation = SimulatedLattice(latt, PREUPDATES);
    Lattice *lattice = simulation.getLattice();

    double runningMag = 0;
    double runningMag2 = 0;
    double runningMag4 = 0;

    for (int i = 0; i < updates; ++i) {
        lattice->updateLattice();
        double magnetization = lattice->getMagnetism();

        runningMag += magnetization;
        runningMag2 += pow(magnetization, 2);
        runningMag4 += pow(magnetization, 4);
    }

    simulation.setAveMag(fabs(runningMag) / updates);
    simulation.setAveMag2(fabs(runningMag2) / updates);
    simulation.setAveMag4(fabs(runningMag4) / updates);

    indMagnetization = simulation.getAveMag();
    indBinderCumulant = simulation.getBinderCumulant();
}
