#include "isingsimulation.h"

using namespace ising;

std::mutex SimulatedLattice::file_mutex;
std::mutex Simulation::file_mutex;
std::mutex Simulation::trial_mutex;
std::mutex Simulation::avgMag_mutex;
std::mutex Simulation::avgMag2_mutex;
std::mutex Simulation::avgMag4_mutex;
std::mutex Simulation::chi0_mutex;
std::mutex Simulation::chiq_mutex;

int main(int argc, char *argv[]) {
    std::string inFilename;
    double t, dt;
    int n, updates, trials;
    char mode;

    receiveSimulationInput(argc, argv, inFilename, t, dt, n, updates, trials,
                           mode);
    manageSimulation(inFilename, t, dt, n, updates, trials, mode);
}

void ising::manageSimulation(const std::string &inFilename, const double t,
                             const double dt, const int n, const int updates,
                             const int trials, const char mode) {
    Simulation simulation(inFilename, t, dt, n, updates, trials, mode);
    simulation.runSimulation();

    dmap temperatures = simulation.getTemperatures();
    dmap magnetizations = simulation.getMagnetizations();
    dmap binderCumulants = simulation.getBinderCumulants();
    dmap correlationFunctions = simulation.getRealCorrelationFunctions();

    std::string outMag = getOutFilename(inFilename, "magnetizations");
    std::string outBC = getOutFilename(inFilename, "binder_cumulants");
    std::string outCL = getOutFilename(inFilename, "correlation_functions");

    writeOutput(outMag, temperatures, magnetizations);
    writeOutput(outBC, temperatures, binderCumulants);
    writeOutput(outCL, temperatures, correlationFunctions);
}

void ising::receiveSimulationInput(int argc, char *argv[],
                                   std::string &filename, double &t, double &dt,
                                   int &n, int &updates, int &trials,
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
        std::cout << "Enter number of trials: ";
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