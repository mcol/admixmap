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

///enum for the various regression types.
/// Mlinear = Multiple linear, Mlogistic = multiple logistic
enum RegressionType {None, Linear, Logistic, Cox, Mlinear, Mlogistic, Multiple};
const std::string RegressionString[] ={"None", "Linear", "Logistic", "Cox", "Mlinear", "Mlogistic", "Multiple"};

///enum for continuous/binary datatypes
enum DataType {Continuous, Binary, CoxData};

#endif /* !defined COMMON_H */
