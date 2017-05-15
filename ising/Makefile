CXX		 = g++
CXXFLAGS = -std=c++14 -Wall

all : ising

ising : ising.o lattice.o
	$(CXX) $(CXXFLAGS) ising.o lattice.o -o ising

ising.o : ising.cpp
	$(CXX) $(CXXFLAGS) -c ising.cpp

lattice.o : lattice.cpp lattice.h
	$(CXX) $(CXXFLAGS) -c lattice.cpp

clean :
	rm -f ising *.o

.PHONY : all clean
