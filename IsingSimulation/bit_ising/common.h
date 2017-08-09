#ifndef COMMON_H_
#define COMMON_H_

#include <array>
#include <map>
#include <vector>

// Typedefs
typedef std::vector<char> cvector;
typedef std::vector<cvector> cvector2;
typedef std::vector<int> ivector;
typedef std::vector<ivector> ivector2;
typedef std::vector<ivector2> ivector3;
typedef std::vector<double> dvector;
typedef std::vector<dvector> dvector2;
typedef std::map<int, ivector> ivectormap;
typedef std::array<int, 2> i2array;
typedef std::map<int, i2array> i2arraymap;
typedef std::vector<i2array> i2arrayvector;

namespace ising {
// Constants
const double E = 2.71828182845904523536;
const double PI = 3.14159265358979323846;
enum { PLUS = 0, MINUS = 1 };
enum { ALL = 'a', PSEUDO = 'p', RANDOM = 'r' };
enum { RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't', STRIANGLE = 'v' };
}

#endif /* COMMON_H_ */