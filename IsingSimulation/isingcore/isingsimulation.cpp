#include "isingsimulation.h"

using namespace ising;

int SimulatedLattice::numLattices = 0;
dvector SimulatedLattice::temperatures;
dvector SimulatedLattice::magnetizations;
dvector SimulatedLattice::binderCumulants;
cvector SimulatedLattice::correlationLengths;
std::mutex SimulatedLattice::data_mutex;
std::mutex SimulatedLattice::file_mutex;

int main(int argc, char *argv[]) {
    std::string inFilename;
    double t, dt;
    int n, updates, trials;
    char mode;

    receiveInputCMD(argc, argv, inFilename, t, dt, n, updates, trials, mode);
    manageSimulations(inFilename, t, dt, n, updates, trials, mode);
}

void ising::manageSimulations(const std::string &inFilename, const double t,
                              const double dt, const int n, const int updates,
                              const int trials, const char mode) {
    std::ifstream inFile(inFilename);

    if (!inFile) {
        std::cout << "Invalid file name. " << inFilename
                  << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    char shape;
    Hamiltonian h = readHamiltonian(inFile, shape);

    std::vector<SimulatedLattice *> simulations;
    std::vector<std::function<void(void)>> simFunctions;

    for (int i = 0; i < n; ++i) {
        double currT = t + i * dt;
        Lattice *lattice = chooseLattice(shape, h, currT, mode, false);
        SimulatedLattice *sim = new SimulatedLattice(
            lattice, inFilename, updates, PREUPDATES, trials);
        simulations.push_back(sim);
    }

    for (auto it = simulations.begin(); it != simulations.end(); ++it) {
        // (*it)->runLatticeSimulation();
        SimulatedLattice &latticeRef = **it;
        simFunctions.push_back([&] { latticeRef.runLatticeSimulation(); });
    }

    runInPool(simFunctions.begin(), simFunctions.end(), getMaxThreads());

    for (auto it = simulations.begin(); it != simulations.end(); ++it) {
        delete *it;
    }

    dvector temperatures = SimulatedLattice::getTemperatures();
    dvector magnetizations = SimulatedLattice::getMagnetizations();
    dvector binderCumulants = SimulatedLattice::getBinderCumulants();
    dvector correlationLengths = SimulatedLattice::getRealCorrelationLengths();

    std::string outMag = getOutFilename(inFilename, "magnetizations");
    std::string outBC = getOutFilename(inFilename, "binder_cumulants");
    std::string outCL = getOutFilename(inFilename, "correlation_lengths");

    writeOutput(outMag, temperatures, magnetizations);
    writeOutput(outBC, temperatures, binderCumulants);
    writeOutput(outCL, temperatures, correlationLengths);
}

void ising::receiveInputCMD(int argc, char *argv[], std::string &filename,
                            double &t, double &dt, int &n, int &updates,
                            int &trials, char &mode) {
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
        std::cout << "Enter number of separate trials per lattice: ";
        std::cin >> trials;
        std::cout
            << "Enter update mode (a - all, p - pseudorandom, r - random): ";
        std::cin >> mode;
    } else if (argc == 8) {
        filename = argv[1];
        t = atof(argv[2]);
        dt = atof(argv[3]);
        n = atoi(argv[4]);
        updates = atoi(argv[5]);
        trials = atoi(argv[6]);
        mode = (argv[7])[0];
    } else {
        std::cout << "Usage: " << argv[0] << " filename(std::string) ";
        std::cout << "min_temperature(float) change_temperature(float) ";
        std::cout << "num_lattices(int) updates(int) trials(int) mode(char)";
        std::cout << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
}

int ising::getMaxThreads() {
    int maxThreads = std::thread::hardware_concurrency();

    if (maxThreads != 0) {
        return maxThreads;
    } else {
        return 4;
    }
}

int ising::getNumThreads(const int remaining) {
    int numThreads = getMaxThreads();

    if (remaining < numThreads) {
        numThreads = remaining;
    }

    return numThreads;
}

template <typename Iter>
void ising::runInPool(Iter begin, Iter end, int threadCount) {
    ThreadPool pool(threadCount);
    for (; begin != end; begin = std::next(begin)) pool.enqueue(*begin);
}