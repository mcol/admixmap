// *-*-C++-*-*
#ifndef COMMON_H
#define COMMON_H 1

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>

typedef std::vector<std::string> Vector_s; //std vector of strings 
typedef std::vector<Vector_s>    Matrix_s; // std vector of std vectors

enum RegressionType {Linear, Logistic, None, Both};
enum DataType {Continuous, Binary};
enum Sex {unknown, male, female};

#endif /* !defined COMMON_H */
