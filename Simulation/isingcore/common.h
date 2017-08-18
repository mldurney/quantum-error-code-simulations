#ifndef COMMON_H_
#define COMMON_H_

#include <array>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <complex>
#include <map>
#include <set>
#include <vector>
#include "randomgenerator.h"

namespace mp = boost::multiprecision;

// Typedefs
typedef unsigned int uint;
typedef mp::number<mp::cpp_dec_float<100>> ldouble;
typedef std::complex<ldouble> cdouble;
typedef std::vector<char> cvector;
typedef std::vector<cvector> cvector2;
typedef std::vector<int> ivector;
typedef std::vector<ivector> ivector2;
typedef std::vector<ivector2> ivector3;
typedef std::vector<ldouble> dvector;
typedef std::vector<dvector> dvector2;
typedef std::vector<cdouble> cdvector;
typedef std::map<int, int> imap;
typedef std::map<int, imap> imap2;
typedef std::map<int, imap2> imap3;
typedef std::map<int, ldouble> dmap;
typedef std::map<int, cdouble> cdmap;
typedef std::map<int, ivector> ivectormap;
typedef std::map<int, dvector> dvectormap;
typedef std::map<int, cdvector> cdvectormap;
typedef std::vector<imap> imapvector;
typedef std::vector<dmap> dmapvector;
typedef std::array<int, 2> i2array;
typedef std::map<int, i2array> i2arraymap;
typedef std::vector<i2array> i2arrayvector;

namespace ising {
// Constants
const double E = 2.71828182845904523536;
const double PI = 3.14159265358979323846;
const double KB = 1.38064852;
enum { PLUS = 0, MINUS = 1 };
enum { ALL = 'a', PSEUDO = 'p', RANDOM = 'r' };
enum { RECTANGLE = 'r', SQUARE = 's', TRIANGLE = 't', STRIANGLE = 'v' };

#if defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32) && !defined(__CYGWIN__)
static const std::string SLASH = "\\";
#else
static const std::string SLASH = "/";
#endif
}

#endif /* COMMON_H_ */