#ifndef LATTICES_H_
#define LATTICES_H_

#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <numeric>
#include <random>
#include <utility>
#include "hamiltonian.h"

typedef std::vector<double> dvector;
typedef std::map<int, int> imap;
typedef std::map<int, double> dmap;

namespace ising {
	enum { PLUS = 0, MINUS = 1 };
	enum { ALL = 'a', PSEUDO = 'p', RANDOM = 'r' };
	enum { RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't', STRIANGLE = 'v' };
	const double E = 2.71828182845904523536;
	const double PI = 3.14159265358979323846;

	class Lattice {
	public:
		Lattice(Hamiltonian h, double t, char m = 'p', bool init = true);
		virtual ~Lattice() = default;

		std::string getType() const { return type; }
		char getShape() const { return shape; }
		double getTemp() const { return temp; }
		char getMode() const { return mode; }
		int getSize() const { return size; }
		const Hamiltonian getHamiltonian() const { return hamiltonian; }
		const ivector2 getHFunction() const { return hFunction; }
		const ivector getIndices() const { return indices; }
		const int getNumIndices() const { return numIndices; }
		const i2arraymap getLocations() const { return locations; }
		const ivectormap getLocalTerms() const { return localTerms; }
		const std::map<int, ivector2> getIndInteractions() const {
			return indInteractions;
		}
		imap getSpins() const { return spins; }
		std::map<int, imap> getXDisplacements() const { return xDisplacements; }
		std::map<int, imap> getYDisplacements() const { return yDisplacements; }
		std::map<int, dmap> getDistances() const { return distances; }

		void updateLattice();
		void switchMode(char m);
		int getTotalEnergy() { return findTotalEnergy(); }
		double getMagnetism() { return findMagnetism(); }

		virtual int findXDisplacement(int, int) { return 1; };
		virtual int findYDisplacement(int, int) { return 1; };
		virtual double findDistance(int, int) { return 1; }
		void reinit() { initSpins(); }
		void flipSpins();
		virtual void print(int cols = -1) const;

	protected:
		const Hamiltonian hamiltonian;
		const ivector2 hFunction;
		const ivector indices;
		ivector randomizedIndices;
		const int numIndices;
		const i2arraymap locations;
		const ivectormap localTerms;
		const std::map<int, ivector2> indInteractions;
		imap spins;
		std::map<int, imap> xDisplacements;
		std::map<int, imap> yDisplacements;
		std::map<int, dmap> distances;

		virtual void checkShape() const {}
		void shapeError() const;
		void setType(std::string t) { type = t; }
		void setMode(char m) { mode = m; }
		void setSize(int s) { size = s; }
		void generateDistances();

		void updateAll();
		void updatePseudo();
		void updateRandom();
		double findProbability(int index);
		int findTotalEnergy();
		inline int findIndexEnergy(int index);
		double findMagnetism();

		float asFloat(uint32_t i);
		float randFloatCO();
		inline unsigned int MWC() { return (zNew() << 16) + wNew(); }

	private:
		std::string type;
		const char shape;
		const double temp;
		char mode;
		int size;
		unsigned int zSeed;
		unsigned int wSeed;

		void initSpins();

		inline unsigned int zNew() {
			return zSeed = 36969 * (zSeed & 65535) + (zSeed >> 16);
		}
		inline unsigned int wNew() {
			return wSeed = 18000 * (wSeed & 65535) + (wSeed >> 16);
		}
	};

	class LatticeFast : public virtual Lattice {
	public:
		LatticeFast(Hamiltonian h, double t, char m = 'p', bool init = true)
			: Lattice(h, t, m, init), coupling((char)h.getHamiltonian()[0][0]) {};
		char getCoupling() const { return coupling; }

	protected:
		void updateAll();
		void updateRandom();
		int findTotalEnergy();
		inline int findIndexEnergy(int index);
		double findMagnetism();

	private:
		const char coupling;
	};

	class RectangularLattice : public virtual Lattice {
	public:
		RectangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
			int r = -1, int c = -1);
		int getRows() const { return rows; }
		int getCols() const { return cols; }
		void print() const { Lattice::print(getCols()); }
		int findXDisplacement(int i, int j);
		int findYDisplacement(int i, int j);
		double findDistance(int i, int j);

	protected:
		void setRows(int r) { rows = r; }
		void setCols(int c) { cols = c; }
		void checkShape() const;

	private:
		int rows;
		int cols;

		void guessRowsCols();
	};

	class RectangularLatticeFast : public RectangularLattice, public LatticeFast {
	public:
		RectangularLatticeFast(Hamiltonian h, double t, char m = 'p',
			bool init = true, int r = -1, int c = -1)
			: Lattice(h, t, m, init),
			RectangularLattice(h, t, m, init, r, c),
			LatticeFast(h, t, m, init) {};
		int findXDisplacement(int i, int j) {
			return RectangularLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return RectangularLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return RectangularLattice::findDistance(i, j);
		}

	protected:
		void checkShape() const { RectangularLattice::checkShape(); }
	};

	class SquareLattice : public RectangularLattice {
	public:
		SquareLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
			int s = -1);
		int getSide() const { return side; }
		void print() const { Lattice::print(getSide()); }
		int findXDisplacement(int i, int j) {
			return RectangularLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return RectangularLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return RectangularLattice::findDistance(i, j);
		}

	protected:
		void setSide(int s) { side = s; }
		void checkShape() const;

	private:
		int side;

		void guessSide();
	};

	class SquareLatticeFast : public SquareLattice, public LatticeFast {
	public:
		SquareLatticeFast(Hamiltonian h, double t, char m = 'p', bool init = true,
			int s = -1)
			: Lattice(h, t, m, init),
			SquareLattice(h, t, m, init, s),
			LatticeFast(h, t, m, init) {};
		int findXDisplacement(int i, int j) {
			return SquareLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return SquareLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return SquareLattice::findDistance(i, j);
		}

	protected:
		void checkShape() const { SquareLattice::checkShape(); }
	};

	class TriangularLattice : public virtual Lattice {
	public:
		TriangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
			int r = -1, int c = -1);
		int getRows() const { return rows; }
		int getCols() const { return cols; }
		void print() const { Lattice::print(getCols()); }
		int findXDisplacement(int i, int j);
		int findYDisplacement(int i, int j);
		double findDistance(int i, int j);

	protected:
		void setRows(int r) { rows = r; }
		void setCols(int c) { cols = c; }
		void checkShape() const;

	private:
		int rows;
		int cols;

		void guessRowsCols();
	};

	class TriangularLatticeFast : public TriangularLattice, public LatticeFast {
	public:
		TriangularLatticeFast(Hamiltonian h, double t, char m = 'p',
			bool init = true, int r = -1, int c = -1)
			: Lattice(h, t, m, init),
			TriangularLattice(h, t, m, init, r, c),
			LatticeFast(h, t, m, init) {};
		int findXDisplacement(int i, int j) {
			return TriangularLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return TriangularLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return TriangularLattice::findDistance(i, j);
		}

	protected:
		void checkShape() const { TriangularLattice::checkShape(); }
	};

	class STriangularLattice : public TriangularLattice {
	public:
		STriangularLattice(Hamiltonian h, double t, char m = 'p', bool init = true,
			int s = -1);
		int getSide() const { return side; }
		void print() const { Lattice::print(getSide()); }
		int findXDisplacement(int i, int j) {
			return TriangularLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return TriangularLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return TriangularLattice::findDistance(i, j);
		}

	protected:
		void setSide(int s) { side = s; }
		void checkShape() const;

	private:
		int side;

		void guessSide();
	};

	class STriangularLatticeFast : public STriangularLattice, public LatticeFast {
	public:
		STriangularLatticeFast(Hamiltonian h, double t, char m = 'p',
			bool init = true, int s = -1)
			: Lattice(h, t, m, init),
			STriangularLattice(h, t, m, init, s),
			LatticeFast(h, t, m, init) {};
		int findXDisplacement(int i, int j) {
			return STriangularLattice::findXDisplacement(i, j);
		}
		int findYDisplacement(int i, int j) {
			return STriangularLattice::findYDisplacement(i, j);
		}
		double findDistance(int i, int j) {
			return STriangularLattice::findDistance(i, j);
		}

	protected:
		void checkShape() const { STriangularLattice::checkShape(); }
	};
}

#endif /* LATTICES_H_ */

/**
Credits to Andy Gainey
(https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers)
and George Marsaglia (http://www.cse.yorku.ca/~oz/marsaglia-rng.html)
for algoritms behind randFloatCO and MWC, respectively, used in random
number generation.
*/

/**
Variables for better random generator:
std::random_device rd;
std::mt19937 gen;
std::uniform_int_distribution<> randSpin;
std::uniform_real_distribution<> randProb;
std::uniform_int_distribution<> randInd;
*/