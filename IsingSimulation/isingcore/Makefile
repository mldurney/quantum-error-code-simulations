CXX		 = g++
CXXFLAGS = -std=c++1z -Wall -pthread -Ofast

all : isingsimulation ising testhamiltonian

debug: CXXFLAGS = -std=c++1z -Wall -pthread -g -ggdb -O0
debug: isingsimulation ising testhamiltonian

isingsimulation : isingsimulation.o isinghelpers.o simulatedlattice.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) isingsimulation.o isinghelpers.o simulatedlattice.o lattices.o hamiltonian.o -o isingsimulation -lstdc++fs

isingsimulation.o : isingsimulation.cpp isingsimulation.h isinghelpers.h simulatedlattice.h threadpool.h
	$(CXX) $(CXXFLAGS) -c isingsimulation.cpp -lstdc++fs

ising : ising.o isinghelpers.o lattices.o hamiltonian.o
	$(CXX) $(CXXFLAGS) ising.o isinghelpers.o lattices.o hamiltonian.o -o ising -lstdc++fs

ising.o : ising.cpp ising.h isinghelpers.h
	$(CXX) $(CXXFLAGS) -c ising.cpp -lstdc++fs

isinghelpers.o : isinghelpers.cpp isinghelpers.h lattices.h
	$(CXX) $(CXXFLAGS) -c isinghelpers.cpp -lstdc++fs

simulatedlattice.o : simulatedlattice.cpp simulatedlattice.h lattices.h
	$(CXX) $(CXXFLAGS) -c simulatedlattice.cpp -lstdc++fs

lattices.o : lattices.cpp lattices.h hamiltonian.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

testhamiltonian : testhamiltonian.o hamiltonian.o
	$(CXX) $(CXXFLAGS) testhamiltonian.o hamiltonian.o -o testhamiltonian

testhamiltonian.o : testhamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c testhamiltonian.cpp

hamiltonian.o : hamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

clean :
	rm -f ising isingsimulation testhamiltonian *.o *.gch *.exe

.PHONY : all clean
