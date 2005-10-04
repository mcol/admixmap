// *-*-C++-*-*
#ifndef COMMON_H
#define COMMON_H 1

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>
#include "vector.h"
#include "vector_i.h"
#include "vector_d.h"
#include "matrix.h"
#include "matrix_i.h"
#include "matrix_d.h"

typedef std::vector<std::string> Vector_s; //std vector of strings 
typedef std::vector<Vector_s>    Matrix_s; // std vector of std vectors

enum RegressionType {Linear, Logistic, None};
enum DataType {Continuous, Binary};
enum Sex {unknown, male, female};

#endif /* !defined COMMON_H */
