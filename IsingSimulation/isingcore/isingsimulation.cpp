#include "isingsimulation.h"

using namespace ising;

int SimulatedLattice::numLattices = 0;
std::vector<double> SimulatedLattice::temperatures;
std::vector<double> SimulatedLattice::magnetizations;
std::vector<double> SimulatedLattice::binderCumulants;
std::mutex SimulatedLattice::data_mutex;

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

	std::vector<SimulatedLattice*> simulations;
	std::vector<std::function<void(void)>> simFunctions;

	for (int trial = 0; trial < n; ++trial) {
		double currT = t + trial * dt;
		Lattice *lattice = chooseLattice(shape, h, currT, mode);
		SimulatedLattice* sim = new SimulatedLattice(lattice, updates, PREUPDATES);
		simulations.push_back(sim);
	}

	for (auto it = simulations.begin(); it != simulations.end(); ++it) {
		(*it)->runLatticeSimulation();
		//SimulatedLattice& latticeRef = **it;
		//simFunctions.push_back([&] { latticeRef.runLatticeSimulation(); });
	}

	//runInPool(simFunctions.begin(), simFunctions.end(), getMaxThreads());

	for (auto it = simulations.begin(); it != simulations.end(); ++it) {
		delete *it;
	}

	std::vector<double> temperatures = SimulatedLattice::getTemperatures();
	std::vector<double> magnetizations = SimulatedLattice::getMagnetizations();
	std::vector<double> binderCumulants = SimulatedLattice::getBinderCumulants();

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
	for (; begin != end; begin = std::next(begin))
		pool.enqueue(*begin);
}
