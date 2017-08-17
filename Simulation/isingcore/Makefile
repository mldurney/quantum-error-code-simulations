CXX		 = g++
CXXFLAGS = -std=c++1z -Wall -Wextra -pedantic -pthread -Ofast

all : testhamiltonian testreplica ising isingsimulation

debug: CXXFLAGS = -std=c++1z -Wall -Wextra -pedantic -pthread -g -ggdb -O0
debug: testhamiltonian testreplica ising isingsimulation

isingsimulation : isingsimulation.o simulation.o threadpoolhelpers.o simulatedlattice.o isinghelpers.o lattices.o replica.o hamiltonian.o
	$(CXX) $(CXXFLAGS) isingsimulation.o simulation.o threadpoolhelpers.o simulatedlattice.o isinghelpers.o lattices.o replica.o hamiltonian.o -o isingsimulation -lstdc++fs

isingsimulation.o : isingsimulation.cpp isingsimulation.h simulation.h threadpoolhelpers.h threadpool.h isinghelpers.h lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c isingsimulation.cpp -lstdc++fs

simulation.o : simulation.cpp simulation.h threadpoolhelpers.h threadpool.h isinghelpers.h lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c simulation.cpp -lstdc++fs

threadpoolhelpers.o : threadpoolhelpers.h threadpool.h
	$(CXX) $(CXXFLAGS) -c threadpoolhelpers.cpp

simulatedlattice.o : simulatedlattice.cpp simulatedlattice.h lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c simulatedlattice.cpp -lstdc++fs

ising : ising.o isinghelpers.o lattices.o replica.o hamiltonian.o
	$(CXX) $(CXXFLAGS) ising.o isinghelpers.o lattices.o replica.o hamiltonian.o -o ising -lstdc++fs

ising.o : ising.cpp ising.h isinghelpers.h lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c ising.cpp -lstdc++fs

isinghelpers.o : isinghelpers.cpp isinghelpers.h lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c isinghelpers.cpp -lstdc++fs

testreplica : testreplica.o lattices.o replica.o hamiltonian.o
	$(CXX) $(CXXFLAGS) testreplica.o lattices.o replica.o hamiltonian.o -o testreplica

testreplica.o : testreplica.cpp lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c testreplica.cpp

lattices.o : lattices.cpp lattices.h replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

replica.o : replica.cpp replica.h properties.h hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c replica.cpp

testhamiltonian : testhamiltonian.o hamiltonian.o
	$(CXX) $(CXXFLAGS) testhamiltonian.o hamiltonian.o -o testhamiltonian

testhamiltonian.o : testhamiltonian.cpp hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c testhamiltonian.cpp

hamiltonian.o : hamiltonian.cpp hamiltonian.h common.h randomgenerator.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

clean :
	rm -f testhamiltonian testreplica ising isingsimulation *.o *.gch *.exe

.PHONY : all clean
