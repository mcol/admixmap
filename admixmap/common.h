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
#include "MatrixArray_d.h"
#include "MatrixArray_i.h"

typedef std::vector<std::string> Vector_s;
typedef std::vector<Vector_s>    Matrix_s;

#endif /* !defined COMMON_H */
