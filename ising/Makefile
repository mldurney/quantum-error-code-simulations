CXX		 = g++
CXXFLAGS = -std=c++14 -Wall -ggdb

all : lattices tests

lattices : lattices.o
	$(CXX) $(CXXFLAGS) lattices.o -o lattices

lattices.o : lattices.cpp lattices.h hamiltonian.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

tests : tests.o hamiltonian.o
	$(CXX) $(CXXFLAGS) tests.o hamiltonian.o -o tests

tests.o : tests.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c tests.cpp

hamiltonian.o : hamiltonian.cpp hamiltonian.h
	$(CXX) $(CXXFLAGS) -c hamiltonian.cpp

clean :
	rm -f lattices tests *.o

.PHONY : all clean
