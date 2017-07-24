#include "stdafx.h"
#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
namespace fs = std::experimental::filesystem;

namespace UnitTest
{		
	TEST_CLASS(UnitTest)
	{
	public:
		
		TEST_METHOD(CheckDistances)
		{
			std::vector<ising::Lattice*> lattices;
			for (auto& p : fs::directory_iterator("../../isingcore/UnitTest/.")) {
				if (p.path().extension() == ".csv") {
					std::ifstream f(p.path());
					char shape = '\0';
					char c = f.peek();

					if (isalpha(c)) {
						shape = c;
						f.ignore(256, '\n');
					}

					ising::Hamiltonian h = ising::Hamiltonian(ising::importHamiltonianVector(f), shape);
					f.close();

					double t = 1;
					char m = ising::PSEUDO;
					ising::Lattice *l = ising::chooseLattice(shape, h, t, m);
					lattices.push_back(l);
				}
			}

			for (auto& l : lattices) {
				auto distances = l->getDistances();
				for (auto const &iIndex : distances) {
					int i = iIndex.first;
					for (auto const &jIndex : iIndex.second) {
						int j = jIndex.first;
						std::ostringstream result;
						result << i << "," << j << ":\t" << jIndex.second;
						Logger::WriteMessage(result.str().c_str());
						double delta = .1;
						Assert::AreEqual(jIndex.second, l->findDistance(i, j), delta);
					}
				}
				std::ostringstream buffer;
				buffer << "Shape: " << l->getType() << "\n\n\n\n\n\n\n\n";
				Logger::WriteMessage(buffer.str().c_str());
			}

			for (auto& l : lattices) {
				delete l;
			}
		}
	};	
}