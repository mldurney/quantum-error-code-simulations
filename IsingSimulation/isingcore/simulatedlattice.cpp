#include "simulatedlattice.h"

using namespace ising;

SimulatedLattice::SimulatedLattice(Lattice *latt, const std::string &filename,
                                   const unsigned int updates,
                                   const unsigned int preupdates,
                                   const unsigned int trials)
    : lattice(latt), updates(updates), preupdates(preupdates), trials(trials) {
    indLattice = numLattices++;
    q = 2 * ising::PI / lattice->getSize();
    initTempFiles(filename);
    temperatures.push_back(-1);
    magnetizations.push_back(-1);
    binderCumulants.push_back(-1);
    correlationLengths.push_back(-1);
}

void SimulatedLattice::initTempFiles(const std::string &filename) {
    tempDirectory = fs::path(filename);
    std::ostringstream latticeName;
    latticeName << lattice->getTemp() << tempDirectory.filename();

    tempDirectory.remove_filename();
    tempDirectory.replace_filename("temp/");
    tempData = tempDirectory / latticeName.str();

    std::unique_lock<std::mutex> lock(file_mutex);
    fs::create_directory(tempDirectory);
    lock.unlock();

    if (fs::exists(tempData)) {
        std::ifstream file(tempData);
        std::string line;
        std::getline(file, line);

        while (std::getline(file, line) && trials > 0) {
            --trials;

            std::stringstream stream(line);
            std::string field;
            dvector fields;

            while (std::getline(stream, field, ',')) {
                fields.push_back(std::stod(field));
            }

            addAvgMag(fields[0]);
            addAvgMag2(fields[1]);
            addAvgMag4(fields[2]);
            addChi0(cdouble(fields[3], fields[4]));
            addChiq(cdouble(fields[5], fields[6]));
        }

        file.close();
    } else {
        std::ofstream file(tempData);
        std::ostringstream header;
        header << "avg_mag,avg_mag2,avg_mag4,chi0_re,chi0_im,chiq_re,chiq_im\n";
        file << header.str();
        file.close();
    }
}

void SimulatedLattice::updateTempFile() {
    std::ofstream file(tempData, std::ofstream::out | std::ofstream::app);
    std::ostringstream row;
    row << avgMag.back() << "," << avgMag2.back() << "," << avgMag4.back()
        << "," << chi0.back().real() << "," << chi0.back().imag() << ","
        << chiq.back().real() << "," << chiq.back().imag() << "\n";
    file << row.str();
    file.close();
}

void SimulatedLattice::runLatticeSimulation() {
    runTrials();

    auto currT = lattice->getTemp();
    auto currM = getAvgMag();
    auto currBC = getBinderCumulant();
    auto currCL = getCorrelationLength();

    std::lock_guard<std::mutex> guard(data_mutex);
    temperatures[indLattice] = currT;
    magnetizations[indLattice] = currM;
    binderCumulants[indLattice] = currBC;
    correlationLengths[indLattice] = currCL;
}

void SimulatedLattice::runPreupdates() {
    for (unsigned int i = 0; i < preupdates; ++i) {
        lattice->updateLattice();
    }
}

void SimulatedLattice::runTrials() {
    for (unsigned int trial = 0; trial < trials; ++trial) {
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

void SimulatedLattice::runUpdates() {
    double runningMag = 0;
    double runningMag2 = 0;
    double runningMag4 = 0;
    std::map<int, imap> runningCorr;

    auto distances = lattice->getDistances();
    auto indices = lattice->getIndices();

    for (unsigned int num1 = 0; num1 < updates; ++num1) {
        for (unsigned int num2 = 0; num2 < SKIP; ++num2) {
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
            runningCorr[i][j] /= (int)updates;
            sumCorrK0 += runningCorr[i][j];
            sumCorrKq += cdouble(runningCorr[i][j]) *
                         std::exp(cdouble(0, q * distances[i][j]));
        }
    }

    addAvgMag(fabs(runningMag) / updates);
    addAvgMag2(fabs(runningMag2) / updates);
    addAvgMag4(fabs(runningMag4) / updates);
    addChi0(sumCorrK0 / cdouble(lattice->getNumIndices()));
    addChiq(sumCorrKq / cdouble(lattice->getNumIndices()));
}

void SimulatedLattice::runUpdatesStable() {
	double runningMag = 0;
	double runningMag2 = 0;
	double runningMag4 = 0;
	std::map<int, imap> runningCorr;

	unsigned int powerMax = 5;
	unsigned int power = reachStability();
	if (power > powerMax) {
		power = powerMax;
	}
	unsigned int cycleUpdates = BASEUPDATES 
			* static_cast<unsigned int>(std::pow(2, power));

	auto distances = lattice->getDistances();
	auto indices = lattice->getIndices();

	for (unsigned int num1 = 0; num1 < cycleUpdates; ++num1) {
		for (unsigned int num2 = 0; num2 < SKIP; ++num2) {
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
				std::exp(cdouble(0, q * distances[i][j]));
		}
	}

	addAvgMag(fabs(runningMag) / cycleUpdates);
	addAvgMag2(fabs(runningMag2) / cycleUpdates);
	addAvgMag4(fabs(runningMag4) / cycleUpdates);
	addChi0(sumCorrK0 / cdouble(lattice->getNumIndices()));
	addChiq(sumCorrKq / cdouble(lattice->getNumIndices()));
}

int SimulatedLattice::reachStability() {
	unsigned int cycleUpdates;
	unsigned int cycle;
	unsigned int maxCycles = 10;
	dvector runningMagBins;
	unsigned int binsToCompare = 3;
	bool continueUpdating = true;

	for (cycle = 1; cycle < binsToCompare; ++cycle) {
		cycleUpdates = BASEUPDATES * static_cast<unsigned int>(std::pow(2, cycle));
		double binMag = 0;

		for (unsigned int num1 = 0; num1 < cycleUpdates; ++num1) {
			for (unsigned int num2 = 0; num2 < SKIP; ++num2) {
				lattice->updateLattice();
			}

			binMag += lattice->getMagnetism();
		}

		runningMagBins.push_back(fabs(binMag) / cycleUpdates);
	}

	while (continueUpdating) {
		cycleUpdates = BASEUPDATES * static_cast<unsigned int>(std::pow(2, cycle));
		dvector cycleMags;

		for (unsigned int num1 = 0; num1 < cycleUpdates; ++num1) {
			for (unsigned int num2 = 0; num2 < SKIP; ++num2) {
				lattice->updateLattice();
			}

			double magnetization = lattice->getMagnetism();
			cycleMags.push_back(magnetization);
		}

		double mean = std::accumulate(cycleMags.begin(), cycleMags.end(), 0.0)
			/ cycleUpdates;
		double std = std::sqrt(std::accumulate(cycleMags.begin(), cycleMags.end(), 0.0,
			[mean](double lhs, double rhs) {
			return rhs + std::pow(lhs - mean, 2);
		}) / cycleUpdates);

		continueUpdating = false;
		for (int i = 1; i < static_cast<int>(binsToCompare); ++i) {
			double difference =
				abs((mean - runningMagBins.end()[-i]) / mean);

			if (difference >= std) {
				continueUpdating = true;
			}
		}

		if (cycle >= maxCycles) {
			continueUpdating = false;
		}

		if (continueUpdating) {
			++cycle;
		}
	}

	return cycle;
}

void SimulatedLattice::addAvgMag(double mag) {
    if (mag >= 0 && mag <= 1) {
        avgMag.push_back(mag);
    } else {
        std::cout << "\nInvalid average magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::addAvgMag2(double mag2) {
    if (mag2 >= 0 && mag2 <= 1) {
        avgMag2.push_back(mag2);
    } else {
        std::cout << "\nInvalid average squared magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

void SimulatedLattice::addAvgMag4(double mag4) {
    if (mag4 >= 0 && mag4 <= 1) {
        avgMag4.push_back(mag4);
    } else {
        std::cout << "\nInvalid average fourth-power magnetization!";
        std::cout << "Must be between 0.0 and 1.0\n\n";
        exit(EXIT_FAILURE);
    }
}

double SimulatedLattice::findAverage(dvector &v) {
	return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

cdouble SimulatedLattice::findAverage(cvector &v) {
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

double SimulatedLattice::findAverageNoOutliers(dvector &v, double percentile1,
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
    double std =
        std::sqrt(std::accumulate(mid.begin(), mid.end(), 0.0,
                                  [mean](double lhs, double rhs) {
                                      return rhs + std::pow(lhs - mean, 2);
                                  }) /
                  mid.size());

    dvector noOutliers;
    for (auto &d : v) {
        if (std::abs(d - mean) <= std) {
            noOutliers.push_back(d);
        }
    }

    return std::accumulate(noOutliers.begin(), noOutliers.end(), 0.0) /
           noOutliers.size();
}

cdouble SimulatedLattice::findAverageNoOutliers(cvector &v, double percentile1,
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

double SimulatedLattice::getBinderCumulant() {
    return 1 - getAvgMag4() / (3 * pow(getAvgMag2(), 2));
}

cdouble SimulatedLattice::getCorrelationLength() {
    return cdouble(.5 * asin(q)) * sqrt((getChi0() / getChiq()) - cdouble(1));
}

dvector SimulatedLattice::getRealCorrelationLengths() {
    dvector realLengths;
    for (auto value : correlationLengths) {
        realLengths.push_back(value.real());
    }
    return realLengths;
}