#include "simulatedlattice.h"

using namespace ising;

SimulatedLattice::SimulatedLattice(Lattice *latt, const std::string &filename,
                                   uint index, uint updates, uint preupdates,
                                   bool suppress)
    : lattice(latt),
      indLattice(index),
      updates(updates),
      preupdates(preupdates) {
    setQ(2 * ising::PI / lattice->getSize());

    if (!suppress) {
        initTempFile(filename);
    }

    auto replicaIndices = lattice->getReplicaIndices();
    for (auto &i : replicaIndices) {
        temperatures[i] = getLattice()->getMinTemperature() +
                          getLattice()->getChangeTemperature() * i;
    }
}

void SimulatedLattice::initTempFile(const std::string &filename) {
    tempDirectory = fs::path(filename);
    std::ostringstream latticeName;
    latticeName << indLattice << tempDirectory.filename().string();

    tempDirectory.remove_filename();
    tempDirectory.replace_filename("temp/");
    tempFile = tempDirectory / latticeName.str();

    std::unique_lock<std::mutex> lock(file_mutex);
    fs::create_directory(tempDirectory);
    lock.unlock();

    std::ofstream file(tempFile);
    std::ostringstream header;
    header << "temperature,avg_mag,avg_mag2,avg_mag4,chi0_re,chi0_im,chiq_re,"
           << "chiq_im\n";
    file << header.str();
    file.close();
}

void SimulatedLattice::updateTempFile() {
    auto replicaIndices = lattice->getReplicaIndices();
    for (auto &i : replicaIndices) {
        std::ofstream file(tempFile, std::ofstream::out | std::ofstream::app);
        std::ostringstream row;
        row << temperatures[i] << "," << avgMag[i] << "," << avgMag2[i] << ","
            << avgMag4[i] << "," << chi0[i].real() << "," << chi0[i].imag()
            << "," << chiq[i].real() << "," << chiq[i].imag() << "\n";
        file << row.str();
        file.close();
    }
}

void SimulatedLattice::runLatticeSimulation() {
    runPreupdates();

    if (updates == 0) {
        runUpdatesStable();
    } else {
        runUpdates();
    }

    updateTempFile();
}

void SimulatedLattice::runPreupdates() {
    for (uint i = 0; i < preupdates; ++i) {
        lattice->ICA();
    }
}

void SimulatedLattice::runUpdates() {
    dmap runningMag, runningMag2, runningMag4;
    imap3 runningCorr;

    auto configs = lattice->getConfigs();
    auto displacements = lattice->getXDisplacements();
    auto indices = lattice->getIndices();

    for (uint num1 = 0; num1 < updates; ++num1) {
        for (uint num2 = 0; num2 < SKIP; ++num2) {
            lattice->ICA();
        }

        for (auto &replicas : configs) {
            uint index = replicas[0]->getReplicaIndex();
            double magnetization = replicas[0]->getMagnetization();
            runningMag[index] += magnetization;
            runningMag2[index] += pow(magnetization, 2);
            runningMag4[index] += pow(magnetization, 4);

            auto spins = replicas[0]->getSpins();
            for (auto &i : indices) {
                for (auto &j : indices) {
                    runningCorr[index][i][j] += spins[i] * spins[j];
                }
            }
        }
    }

    cdmap sumCorrK0, sumCorrKq;
    for (auto &replicas : configs) {
        for (auto &i : indices) {
            for (auto &j : indices) {
                uint index = replicas[0]->getReplicaIndex();
                runningCorr[index][i][j] /= (int)updates;
                sumCorrK0[index] += runningCorr[index][i][j];
                sumCorrKq[index] +=
                    cdouble(runningCorr[index][i][j]) *
                    std::exp(cdouble(0, q * displacements[i][j]));
            }
        }
    }

    auto replicaIndices = lattice->getReplicaIndices();
    for (auto &i : replicaIndices) {
        addAvgMag(i, fabs(runningMag[i]) / updates);
        addAvgMag2(i, fabs(runningMag2[i]) / updates);
        addAvgMag4(i, fabs(runningMag4[i]) / updates);
        addChi0(i, sumCorrK0[i] / cdouble(lattice->getNumIndices()));
        addChiq(i, sumCorrKq[i] / cdouble(lattice->getNumIndices()));
    }
}

void SimulatedLattice::runUpdatesStable() {
    dmap runningMag, runningMag2, runningMag4;
    imap3 runningCorr;

    auto configs = lattice->getConfigs();
    auto displacements = lattice->getXDisplacements();
    auto indices = lattice->getIndices();

    uint powerMax = 5;
    uint power = reachStability();
    if (power > powerMax) {
        power = powerMax;
    }
    uint cycleUpdates = BASEUPDATES * static_cast<uint>(std::pow(2, power));

    for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
        for (uint num2 = 0; num2 < SKIP; ++num2) {
            lattice->ICA();
        }

        for (auto &replicas : configs) {
            uint index = replicas[0]->getReplicaIndex();
            double magnetization = replicas[0]->getMagnetization();
            runningMag[index] += magnetization;
            runningMag2[index] += pow(magnetization, 2);
            runningMag4[index] += pow(magnetization, 4);

            auto spins = replicas[0]->getSpins();
            for (auto &i : indices) {
                for (auto &j : indices) {
                    runningCorr[index][i][j] += spins[i] * spins[j];
                }
            }
        }
    }

    cdmap sumCorrK0, sumCorrKq;
    for (auto &replicas : configs) {
        for (auto &i : indices) {
            for (auto &j : indices) {
                uint index = replicas[0]->getReplicaIndex();
                runningCorr[index][i][j] /= cycleUpdates;
                sumCorrK0[index] += runningCorr[index][i][j];
                sumCorrKq[index] +=
                    cdouble(runningCorr[index][i][j]) *
                    std::exp(cdouble(0, q * displacements[i][j]));
            }
        }
    }

    auto replicaIndices = lattice->getReplicaIndices();
    for (auto &i : replicaIndices) {
        addAvgMag(i, fabs(runningMag[i]) / cycleUpdates);
        addAvgMag2(i, fabs(runningMag2[i]) / cycleUpdates);
        addAvgMag4(i, fabs(runningMag4[i]) / cycleUpdates);
        addChi0(i, sumCorrK0[i] / cdouble(lattice->getNumIndices()));
        addChiq(i, sumCorrKq[i] / cdouble(lattice->getNumIndices()));
    }
}

int SimulatedLattice::reachStability() {
    auto configs = lattice->getConfigs();
    auto replicaIndices = lattice->getReplicaIndices();

    uint cycleUpdates;
    uint cycle;
    uint maxCycles = 10;
    dmapvector bins;
    uint binsToCompare = 3;
    bool continueUpdating = true;

    for (cycle = 1; cycle < binsToCompare; ++cycle) {
        cycleUpdates = BASEUPDATES * static_cast<uint>(std::pow(2, cycle));
        dmap mags;

        for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
            for (uint num2 = 0; num2 < SKIP; ++num2) {
                lattice->ICA();
            }

            for (auto &i : replicaIndices) {
                mags[i] += configs[i][0]->getMagnetization();
            }
        }

        for (auto &m : mags) {
            m.second = fabs(m.second) / cycleUpdates;
        }

        bins.push_back(mags);
    }

    while (continueUpdating) {
        cycleUpdates = BASEUPDATES * static_cast<uint>(std::pow(2, cycle));
        dvectormap cycleMags;

        for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
            for (uint num2 = 0; num2 < SKIP; ++num2) {
                lattice->ICA();
            }

            for (auto &i : replicaIndices) {
                double magnetization = configs[i][0]->getMagnetization();
                cycleMags[i].push_back(magnetization);
            }
        }

        dmap means, stds;

        for (auto &i : replicaIndices) {
            means[i] =
                std::accumulate(cycleMags[i].begin(), cycleMags[i].end(), 0.0) /
                cycleUpdates;
            stds[i] = std::sqrt(
                std::accumulate(cycleMags[i].begin(), cycleMags[i].end(), 0.0,
                                [&](double lhs, double rhs) {
                                    return rhs + std::pow(lhs - means[i], 2);
                                }) /
                cycleUpdates);
        }

        continueUpdating = false;

        for (auto &i : replicaIndices) {
            if (cycle >= maxCycles || continueUpdating == true) {
                break;
            }

            for (int j = 1; j < static_cast<int>(binsToCompare); ++j) {
                double difference =
                    abs((means[i] - bins.end()[-j][i]) / means[i]);

                if (difference >= stds[i]) {
                    continueUpdating = true;
                }
            }
        }

        if (continueUpdating) {
            bins.push_back(means);
            ++cycle;
        }
    }

    return cycle;
}

int SimulatedLattice::reachStabilityChi0() {
    auto configs = lattice->getConfigs();
    auto displacements = lattice->getXDisplacements();
    auto indices = lattice->getIndices();
    auto numIndices = lattice->getNumIndices();
    auto replicaIndices = lattice->getReplicaIndices();

    uint cycleUpdates;
    uint cycle;
    uint maxCycles = 10;
    dmapvector bins;
    uint binsToCompare = 3;
    bool continueUpdating = true;

    for (cycle = 1; cycle < binsToCompare; ++cycle) {
        cycleUpdates = BASEUPDATES * static_cast<uint>(std::pow(2, cycle));
        dmap chi0s;
        imap3 runningCorr;

        for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
            for (uint num2 = 0; num2 < SKIP; ++num2) {
                lattice->ICA();
            }

            for (auto &i : replicaIndices) {
                auto spins = configs[i][0]->getSpins();
                for (auto &j : indices) {
                    for (auto &k : indices) {
                        runningCorr[i][j][k] += spins[j] * spins[k];
                    }
                }
            }
        }

        for (auto &i : replicaIndices) {
            double sumCorrK0 = 0;

            for (auto &j : indices) {
                for (auto &k : indices) {
                    sumCorrK0 += runningCorr[i][j][k] / (int)cycleUpdates;
                }
            }

            chi0s[i] = sumCorrK0 / numIndices;
        }

        bins.push_back(chi0s);
    }

    while (continueUpdating) {
        cycleUpdates = BASEUPDATES * static_cast<uint>(std::pow(2, cycle));
        dvectormap cycleChi0s;

        for (uint num1 = 0; num1 < cycleUpdates; ++num1) {
            for (uint num2 = 0; num2 < SKIP; ++num2) {
                lattice->ICA();
            }

            for (auto &i : replicaIndices) {
                auto spins = configs[i][0]->getSpins();
                double sumCorrK0 = 0;

                for (auto &j : indices) {
                    for (auto &k : indices) {
                        sumCorrK0 += spins[j] * spins[k];
                    }
                }

                cycleChi0s[i].push_back(sumCorrK0 / numIndices);
            }
        }

        dmap means, stds;

        for (auto &i : replicaIndices) {
            means[i] = std::accumulate(cycleChi0s[i].begin(),
                                       cycleChi0s[i].end(), 0.0) /
                       cycleUpdates;
            stds[i] = std::sqrt(
                std::accumulate(cycleChi0s[i].begin(), cycleChi0s[i].end(), 0.0,
                                [&](double lhs, double rhs) {
                                    return rhs + std::pow(lhs - means[i], 2);
                                }) /
                cycleUpdates);
        }

        continueUpdating = false;

        for (auto &i : replicaIndices) {
            if (cycle >= maxCycles || continueUpdating == true) {
                break;
            }

            for (int j = 1; j < static_cast<int>(binsToCompare); ++j) {
                double difference =
                    abs((means[i] - bins.end()[-j][i]) / means[i]);

                if (difference >= stds[i]) {
                    continueUpdating = true;
                }
            }
        }

        if (continueUpdating) {
            bins.push_back(means);
            ++cycle;
        }
    }

    return cycle;
}

void SimulatedLattice::addAvgMag(uint index, double mag) {
    if (mag >= 0 && mag <= 1) {
        avgMag[index] = mag;
    } else {
        std::cout << "\nInvalid average magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::addAvgMag2(uint index, double mag2) {
    if (mag2 >= 0 && mag2 <= 1) {
        avgMag2[index] = mag2;
    } else {
        std::cout << "\nInvalid average squared magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::addAvgMag4(uint index, double mag4) {
    if (mag4 >= 0 && mag4 <= 1) {
        avgMag4[index] = mag4;
    } else {
        std::cout << "\nInvalid average fourth-power magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}