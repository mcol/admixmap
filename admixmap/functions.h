// *-*-C++-*-*
#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

#include <vector>
#include "vector_d.h"

double getGammaLogDensity(double alpha, double beta, double x);

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x);

double AverageOfLogs(const std::vector<double>& vec, double max);

#endif /* !FUNCTIONS_H */
