CXX		 = g++
CXXFLAGS = -std=c++14 -Wall -ggdb -pthread -Ofast

all : isingsimulation ising testhamiltonian

isingsimulation : isingsimulation.o isinghelpers.o simulatedlattice.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) isingsimulation.o isinghelpers.o simulatedlattice.o lattices.o hamiltonian.o -o isingsimulation

isingsimulation.o : isingsimulation.cpp isingsimulation.h isinghelpers.h simulatedlattice.h
	$(CXX) $(CXXFLAGS) -c isingsimulation.cpp

ising : ising.o isinghelpers.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) ising.o isinghelpers.o lattices.o hamiltonian.o -o ising

ising.o : ising.cpp ising.h isinghelpers.h
	$(CXX) $(CXXFLAGS) -c ising.cpp

isinghelpers.o : isinghelpers.cpp isinghelpers.h lattices.h
	$(CXX) $(CXXFLAGS) -c isinghelpers.cpp

simulatedlattice.o : simulatedlattice.cpp simulatedlattice.h lattices.h
	$(CXX) $(CXXFLAGS) -c simulatedlattice.cpp

lattices.o : lattices.cpp lattices.h hamiltonian.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

testhamiltonian : testhamiltonian.o hamiltonian.o
	$(CXX) $(CXXFLAGS) testhamiltonian.o hamiltonian.o -o testhamiltonian

testhamiltonian.o : testhamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c testhamiltonian.cpp

hamiltonian.o : hamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

memcheck.o: memcheck.cpp memcheck.h
	$(CXX) $(CXXFLAGS) -c memcheck.cpp

clean :
	rm -f ising isingsimulation testhamiltonian *.o *.gch

.PHONY : all clean
