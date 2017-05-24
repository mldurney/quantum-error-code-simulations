CXX		 = g++
CXXFLAGS = -std=c++14 -Wall

all : ising ising_simulation

ising : ising.o lattice_near.o
	$(CXX) $(CXXFLAGS) ising.o lattice_near.o -o ising

ising.o : ising.cpp
	$(CXX) $(CXXFLAGS) -c ising.cpp

ising_simulation : ising_simulation.o lattice_near.o
	$(CXX) $(CXXFLAGS) ising_simulation.o lattice_near.o -o ising_simulation

ising_simulation.o : ising_simulation.cpp
	$(CXX) $(CXXFLAGS) -c ising_simulation.cpp

lattice_near.o : lattice_near.cpp lattice_near.h
	$(CXX) $(CXXFLAGS) -c lattice_near.cpp

clean :
	rm -f ising ising_simulation *.o

.PHONY : all clean
