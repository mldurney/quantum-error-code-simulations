CXX		 = g++
CXXFLAGS = -std=c++14 -Wall -ggdb

all : ising ising_simulation

ising : ising.o lattices.o
	$(CXX) $(CXXFLAGS) ising.o lattices.o -o ising

ising.o : ising.cpp
	$(CXX) $(CXXFLAGS) -c ising.cpp

ising_simulation : ising_simulation.o lattices.o
	$(CXX) $(CXXFLAGS) ising_simulation.o lattices.o -o ising_simulation

ising_simulation.o : ising_simulation.cpp
	$(CXX) $(CXXFLAGS) -c ising_simulation.cpp

lattices.o : lattices.cpp lattices.h
	$(CXX) $(CXXFLAGS) -c lattices.cpp

clean :
	rm -f ising ising_simulation *.o

.PHONY : all clean
