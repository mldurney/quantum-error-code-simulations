#include "isingsimulation.h"

int main(int argc, char* argv[])
{
    string inFilename;
    double t, dt;
    int n, updates;
    char mode, shape;

    receiveInput(argc, argv, inFilename, t, dt, n, updates, mode);

    ifstream inFile(inFilename);

    if (!inFile)
    {
        cout << "Invalid file name. " << argv[1] << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    Hamiltonian h = readHamiltonian(inFile, shape);

    vector<thread> threads;
    vector<double> temperatures;
    vector<double> magnetizations(n, 0);
    vector<double> binderCumulants(n, 0);

    Lattice *lattice;
    for (int i = 0; i < n; ++i)
    {
        double currT = t + i * dt;
        temperatures.push_back(currT);

        switch (shape)
        {
            case RECTANGLE:
                lattice = new RectangularLattice(h, currT, mode);
                break;
            case SQUARE:
                lattice = new SquareLattice(h, currT, mode);
                break;
            case TRIANGLE:
                lattice = new TriangularLattice(h, currT, mode);
                break;
            default:
                lattice = new Lattice(h, currT, mode);
        }

        threads.emplace_back(runLatticeSimulation, lattice, updates,
                &magnetizations[i], &binderCumulants[i]);
    }

    for (auto &t : threads)
    {
        t.join();
    }

    string oldDir = "hamiltonians";
    string outMag = getOutFilename(inFilename, oldDir, "magnetizations");
    string outBinder = getOutFilename(inFilename, oldDir, "binder_cumulants");

    writeOutput(outMag, temperatures, magnetizations);
    writeOutput(outBinder, temperatures, binderCumulants);
}

void receiveInput(int argc, char* argv[], string& filename, double& t,
        double& dt, int& n, int& updates, char& mode)
{
    if (argc == 1)
    {
        cout << "Enter Hamiltonian input file: ";
        cin >> filename;
        cout << "Enter minimum temperature (K): ";
        cin >> t;
        cout << "Enter change in temperature between lattices (K): ";
        cin >> dt;
        cout << "Enter total number of lattices: ";
        cin >> n;
        cout << "Enter number of updates per lattice: ";
        cin >> updates;
        cout << "Enter update mode (a - all, p - pseudorandom, r - random): ";
        cin >> mode;
    }
    else if (argc == 7)
    {
        filename = argv[1];
        t = atof(argv[2]);
        dt = atof(argv[3]);
        n = atoi(argv[4]);
        updates = atoi(argv[5]);
        mode = (argv[6])[0];
    }
    else
    {
        cout << "Usage: " << argv[0] << " filename(string) ";
        cout << "min_temperature(float) change_temperature(float) ";
        cout << "num_lattices(int) updates(int) mode(char)\n\n";
        exit(EXIT_FAILURE);
    }
}

void runLatticeSimulation(Lattice* latt, int updates,
        double* magnetizationPtr, double* binderCumulantPtr)
{
    SimulatedLattice simulation = SimulatedLattice(latt, PREUPDATES);
    Lattice *lattice = simulation.getLattice();

    double runningMag = 0;
    double runningMag2 = 0;
    double runningMag4 = 0;

    for (int i = 0; i < updates; ++i)
    {
        lattice->updateLattice();
        double magnetization = lattice->getMagnetism();

        runningMag += magnetization;
        runningMag2 += pow(magnetization, 2);
        runningMag4 += pow(magnetization, 4);
    }

    simulation.setAveMag(abs(runningMag) / updates);
    simulation.setAveMag2(abs(runningMag2) / updates);
    simulation.setAveMag4(abs(runningMag4) / updates);

    *magnetizationPtr = simulation.getAveMag();
    *binderCumulantPtr = simulation.getBinderCumulant();
}
