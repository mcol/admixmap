// *-*-C++-*-*

/// \file common.h
///
/// @warning
/// <B>WARNING!</B> this file duplicates, in the global namespace, symbols that
/// are also defined in @a <bclib/common.h>


#ifndef COMMON_H
#define COMMON_H 1

#include <vector>
#include <string>

typedef std::vector<std::string> Vector_s; //std vector of strings 
typedef std::vector<Vector_s>    Matrix_s; // std vector of std vectors

///enum for the various regression types.
/// Mlinear = Multiple linear, Mlogistic = multiple logistic
enum RegressionType {None, Linear, Logistic, Cox, Mlinear, Mlogistic, Multiple};
const std::string RegressionString[] ={"None", "Linear", "Logistic", "Cox", "Mlinear", "Mlogistic", "Multiple"};

///enum for continuous/binary datatypes
enum DataType {Continuous, Binary, CoxData};

#ifdef PARALLEL
#include "mpi2c++/mpi++.h"
#endif

#endif /* !defined COMMON_H */
