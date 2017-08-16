#include "simulation.h"

using namespace ising;

Simulation::Simulation(const std::string &filename, double t, double dt, uint n,
                       uint updates, uint trials, char mode)
    : Simulation(filename, t, dt, n, updates, PREUPDATES, trials, mode) {}

Simulation::Simulation(const std::string &filename, double t, double dt, uint n,
                       uint updates, uint preupdates, uint trials, char mode)
    : inFilename(filename),
      minT(t),
      dT(dt),
      numT(n),
      updates(updates),
      preupdates(preupdates),
      trials(trials),
      mode(mode) {
    for (unsigned i = 0; i < trials; ++i) {
        remainingTrials.push_back(i);
    }

    for (unsigned i = 0; i < n; ++i) {
        double currT = minT + dt * n;
        temperatures[temperatureToIndex(currT)] = currT;
    }

    checkInputFile();
    initTempDirectory();
    loadTempData();
    initLattices();
}

void Simulation::runSimulation() {
    runTrials();

    for (int i = 0; i < numT; ++i) {
        magnetizations[i] = getAvgMag();
        binderCumulants[i] = getBinderCumulant();
        correlationFunctions[i] = getCorrelationFunctions();
    }
}

void Simulation::checkInputFile() {
    std::ifstream file(inFilename);

    if (!file) {
        std::cout << "Invalid file name. " << inFilename
                  << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }

    file.close();
}

void Simulation::initTempDirectory() {
    tempDirectory = fs::path(inFilename);
    tempDirectory.remove_filename();
    tempDirectory.replace_filename("temp/");
    fs::create_directory(tempDirectory);
}

void Simulation::loadTempData() {
    std::vector<fs::path> tempFiles;
    for (auto &p : fs::directory_iterator(tempDirectory)) {
        if (p.path().extension == ".csv" &&
            p.path().filename.c_str()[0].isdigit()) {
            tempFiles.push_back(p);
        }
    }

    for (auto &p : tempFiles) {
        uint trial = getLeadingInt(p);

        auto it =
            std::find(remainingTrials.begin(), remainingTrials.end(), trial);
        if (it == remainingTrials.end()) {
            std::cout << "Key not in remaining trials for file " << p.c_str()
                      << "! Exiting...\n\n";
            exit(EXIT_FAILURE);
        } else {
            remainingTrials.erase(it);
        }

        std::ifstream file(p);
        std::string line;
        double num;
        uint unreadT = numT;

        while (isalpha(file.peek())) {
            file.ignore(256, '\n');
        }

        while (getline(file, line)) {
            dvector data;
            std::istringstream lineStream(line);

            while (lineStream >> num) {
                data.push_back(num);

                if (lineStream.peek() == ',') {
                    lineStream.ignore();
                }
            }

            if (data.size() != 8) {
                std::cout << "Insufficient number of entries in row of "
                          << p.c_str() << "! Exiting...\n\n";
                exit(EXIT_FAILURE);
            }

            uint index = temperatureToIndex(data[0]);
            --unreadT;

            avgMag[index].push_back(data[1]);
            avgMag2[index].push_back(data[2]);
            avgMag4[index].push_back(data[3]);
            chi0[index].push_back((data[4], data[5]));
            chiq[index].push_back((data[6], data[7]));
        }

        if (unreadT != 0) {
            std::cout << "Not all temperatures read from file " << p.c_str()
                      << "! Exiting...\n\n";
            exit(EXIT_FAILURE);
        }
    }
}

uint Simulation::getLeadingInt(const fs::path &filename) {
    std::string name = filename.filename;
    std::string trialString;

    for (char &c : name) {
        if (isdigit(c)) {
            trialString += c;
        } else {
            break;
        }
    }

    return stoi(trialString);
}

uint Simulation::temperatureToIndex(double t) {
    if (t < (minT - .5 * dT)) {
        std::cout << "INVALID TEMPERATURE -- below minimum. Exiting... ";
        exit(EXIT_FAILURE);
    }

    return (uint)((t - minT) / dT + .5);
}

void Simulation::initLattices() {
    std::vector<std::function<void(void)>> initializers;

    for (auto &trial : remainingTrials) {
        initializers.push_back([&] { addLattice(trial); });
    }

    runInPool(initializers.begin(), initializers.end(), getMaxThreads());
}

void Simulation::addLattice(uint trial) {
    std::unique_lock<std::mutex> lock(file_mutex);
    std::ifstream file(inFilename);

    char shape;
    Hamiltonian hamiltonian = readHamiltonian(file, shape);

    file.close();
    lock.unlock();

    Lattice *lattice = chooseLattice(shape, hamiltonian, minT, dT, numT, mode);
    auto simLattice = std::make_unique<SimulatedLattice>(
        lattice, inFilename, trial, updates, preupdates);

    std::lock_guard<std::mutex> guard(trial_mutex);
    lattices[trial] = std::move(simLattice);
}

void Simulation::runTrials() {
    std::vector<std::function<void(void)>> runs;

    while (!remainingTrials.empty()) {
        auto &trial = remainingTrials.back();
        remainingTrials.pop_back();
        runs.push_back([&] { runTrial(trial); });
    }

    runInPool(runs.begin(), runs.end(), getMaxThreads());
}

void Simulation::runTrial(uint trial) {
    lattices[trial]->runLatticeSimulation();

    auto lAvgMag = lattices[trial]->getAvgMag();
    auto lAvgMag2 = lattices[trial]->getAvgMag2();
    auto lAvgMag4 = lattices[trial]->getAvgMag4();
    auto lChi0 = lattices[trial]->getChi0();
    auto lChiq = lattices[trial]->getChiq();

    std::unique_lock<std::mutex> lock;

    lock = std::unique_lock<std::mutex>(avgMag_mutex);
    for (auto &mag : lAvgMag) {
        addAvgMag(mag.first, mag.second);
    }
    lock.unlock;

    lock = std::unique_lock<std::mutex>(avgMag2_mutex);
    for (auto &mag : lAvgMag2) {
        addAvgMag2(mag.first, mag.second);
    }
    lock.unlock;

    lock = std::unique_lock<std::mutex>(avgMag4_mutex);
    for (auto &mag : lAvgMag4) {
        addAvgMag4(mag.first, mag.second);
    }
    lock.unlock;

    lock = std::unique_lock<std::mutex>(chi0_mutex);
    for (auto &x : lChi0) {
        addChi0(x.first, x.second);
    }
    lock.unlock;

    lock = std::unique_lock<std::mutex>(chiq_mutex);
    for (auto &x : lChiq) {
        addChiq(x.first, x.second);
    }
    lock.unlock;
}

void Simulation::addAvgMag(uint index, double mag) {
    if (mag >= 0 && mag <= 1) {
        avgMag[index].push_back(mag);
    } else {
        std::cout << "\nInvalid average magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void Simulation::addAvgMag2(uint index, double mag2) {
    if (mag2 >= 0 && mag2 <= 1) {
        avgMag2[index].push_back(mag2);
    } else {
        std::cout << "\nInvalid average squared magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void Simulation::addAvgMag4(uint index, double mag4) {
    if (mag4 >= 0 && mag4 <= 1) {
        avgMag4[index].push_back(mag4);
    } else {
        std::cout << "\nInvalid average fourth-power magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

double Simulation::findAverage(dvector &v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

cdouble Simulation::findAverage(cdvector &v) {
    dvector vReal;
    dvector vImag;

    for (cdouble &c : v) {
        vReal.push_back(c.real());
        vImag.push_back(c.imag());
    }

    double avgReal = findAverage(vReal);
    double avgImag = findAverage(vImag);

    return cdouble(avgReal, avgImag);
}

double Simulation::findAverageNoOutliers(dvector &v, double percentile1,
                                         double percentile2) {
    if (percentile1 < 0 || percentile2 < 0 || percentile1 > 1 ||
        percentile2 > 1) {
        std::cout << "\nInvalid percentile! Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }

    dvector mid;
    double minPercentile = std::min(percentile1, percentile2);
    double maxPercentile = std::max(percentile1, percentile2);

    std::sort(v.begin(), v.end());

    for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i) {
        double percentile = static_cast<double>(i) / v.size();
        if (percentile <= maxPercentile && percentile >= minPercentile) {
            mid.push_back(v[i]);
        }
    }

    double mean = std::accumulate(mid.begin(), mid.end(), 0.0) / mid.size();
    double std = std::sqrt(
        std::abs(std::accumulate(mid.begin(), mid.end(), 0.0,
                                 [mean](double lhs, double rhs) {
                                     return rhs + std::pow(lhs - mean, 2);
                                 })) /
        mid.size());

    double minStd = 1e-6;
    std = std::max(std, minStd);

    dvector noOutliers;
    for (auto &d : v) {
        if (std::abs(d - mean) <= std) {
            noOutliers.push_back(d);
        }
    }

    return std::accumulate(noOutliers.begin(), noOutliers.end(), 0.0) /
           noOutliers.size();
}

cdouble Simulation::findAverageNoOutliers(cdvector &v, double percentile1,
                                          double percentile2) {
    dvector vReal;
    dvector vImag;

    for (cdouble &c : v) {
        vReal.push_back(c.real());
        vImag.push_back(c.imag());
    }

    double avgReal = findAverageNoOutliers(vReal, percentile1, percentile2);
    double avgImag = findAverageNoOutliers(vImag, percentile1, percentile2);
    return cdouble(avgReal, avgImag);
}

double SimulatedLattice::getBinderCumulant(uint n) {
    return 1 - getAvgMag4(n) / (3 * pow(getAvgMag2(n), 2));
}

cdouble SimulatedLattice::getCorrelationFunction(uint n) {
    return cdouble(1 / (2 * lattices[n]->getSize() * sin(q)) *
                   sqrt((getChi0(n) / getChiq(n)) - cdouble(1)));
}

dmap SimulatedLattice::getRealCorrelationFunctions() {
    dmap realLengths;
    for (auto &value : correlationFunctions) {
        realLengths[value.first] = value.second.real();
    }
    return realLengths;
}

#include "simulatedlattice.h"

/*
using namespace ising;

void SimulatedLattice::runLatticeSimulation() {
    runTrials();

    auto currT = lattice->getTemp();
    auto currM = getAvgMag();
    auto currBC = getBinderCumulant();
    auto currCL = getCorrelationFunction();

    std::lock_guard<std::mutex> guard(data_mutex);
    temperatures[indLattice] = currT;
    magnetizations[indLattice] = currM;
    binderCumulants[indLattice] = currBC;
    correlationFunctions[indLattice] = currCL;
}

void SimulatedLattice::runTrials() {
    for (uint trial = 0; trial < trials; ++trial) {
        lattice->reinit();
        runPreupdates();

        if (updates == 0) {
            runUpdatesStable();
        } else {
            runUpdates();
        }

        updateTempFile();
    }
}

void SimulatedLattice::runUpdatesStable() {
    double runningMag = 0;
    double runningMag2 = 0;
    double runningMag4 = 0;
    std::map<int, imap> runningCorr;

    uint powerMax = 5;
    uint power = reachStability();
    if (power > powerMax) {
        power = powerMax;
    }
    uint cycleUpdates =
        BASEUPDATES * static_cast<uint>(std::pow(2, power));

    auto displacements = lattice->getXDisplacements();
    auto indices = lattice->getIndices();

    for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
        for (uint num2 = 0; num2 < SKIP; ++num2) {
            lattice->updateLattice();
        }

        double magnetization = lattice->getMagnetism();
        runningMag += magnetization;
        runningMag2 += pow(magnetization, 2);
        runningMag4 += pow(magnetization, 4);

        imap spins = lattice->getSpins();
        for (auto &i : indices) {
            for (auto &j : indices) {
                runningCorr[i][j] += spins[i] * spins[j];
            }
        }
    }

    cdouble sumCorrK0 = 0;
    cdouble sumCorrKq = 0;
    for (auto &i : indices) {
        for (auto &j : indices) {
            runningCorr[i][j] /= (int)cycleUpdates;
            sumCorrK0 += runningCorr[i][j];
            sumCorrKq += cdouble(runningCorr[i][j]) *
                         std::exp(cdouble(0, q * displacements[i][j]));
        }
    }

    addAvgMag(fabs(runningMag) / cycleUpdates);
    addAvgMag2(fabs(runningMag2) / cycleUpdates);
    addAvgMag4(fabs(runningMag4) / cycleUpdates);
    addChi0(sumCorrK0 / cdouble(lattice->getNumIndices()));
    addChiq(sumCorrKq / cdouble(lattice->getNumIndices()));
}
*/