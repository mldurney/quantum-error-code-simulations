CXX		 = g++
CXXFLAGS = -std=c++14 -Wall -ggdb

all : ising ising_simulation test_hamiltonian

ising : ising.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) ising.o lattices.o hamiltonian.o -o ising

ising.o : ising.cpp lattices.h
	$(CXX) $(CXXFLAGS) -c ising.cpp

ising_simulation : ising_simulation.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) ising_simulation.o lattices.o hamiltonian.o -o ising_simulation

ising_simulation.o : ising_simulation.cpp lattices.h
	$(CXX) $(CXXFLAGS) -c ising_simulation.cpp

lattices.o : lattices.cpp lattices.h hamiltonian.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

test_hamiltonian : test_hamiltonian.o hamiltonian.o
	$(CXX) $(CXXFLAGS) test_hamiltonian.o hamiltonian.o -o test_hamiltonian

test_hamiltonian.o : test_hamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c test_hamiltonian.cpp

hamiltonian.o : hamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

memcheck.o: memcheck.cpp memcheck.h
	$(CXX) $(CXXFLAGS) -c memcheck.cpp

clean :
	rm -f ising ising_simulation lattices tests test_hamiltonian *.o *.gch

.PHONY : all clean
