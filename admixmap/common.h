// *-*-C++-*-*
#ifndef COMMON_H
#define COMMON_H 1

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>

typedef std::vector<std::string> Vector_s; //std vector of strings 
typedef std::vector<Vector_s>    Matrix_s; // std vector of std vectors

typedef std::vector<std::vector<unsigned short> > genotype;

enum RegressionType {Linear, Logistic, None, Both};
enum DataType {Continuous, Binary};
enum Sex {unknown, male, female};

//uncomment the next line for parallel version
//#define PARALLEL
#ifdef PARALLEL
#include "mpi++.h"
#endif

#endif /* !defined COMMON_H */
