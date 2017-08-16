CXX		 = g++
CXXFLAGS = -std=c++1z -Wall -Wextra -pedantic -Ofast

all : testhamiltonian testreplica ising

debug: CXXFLAGS = -std=c++1z -Wall -Wextra -pedantic -g -ggdb -O0
debug: testhamiltonian testreplica ising

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
	rm -f testhamiltonian testreplica ising *.o *.gch *.exe

.PHONY : all clean
