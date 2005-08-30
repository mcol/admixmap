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


// typedef double ** dmatrix;
// typedef int ** imatrix;


// typedef struct 
// {
//    int alleles[2];
// } genotype; // genotype contains array of 2 integers

// typedef std::vector<genotype> Vector_g; // vector of genotypes over loci

// typedef std::vector<Vector_g> Matrix_g; // vector of Vector g objects over individuals 

#endif /* !defined COMMON_H */
