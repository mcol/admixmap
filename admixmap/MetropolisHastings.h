// *-*-C++-*-*
#ifndef ADMIX_METROPOLISHASTINGS_H
#define ADMIX_METROPOLISHASTINGS_H 1

#include "rand.h"
#include "vector_d.h"
#include "MatrixArray_i.h"
#include "MatrixArray_d.h"

class MetropolisHastings
{
private: // members
  Vector_d parameters;
  MatrixArray_i data_i;
  MatrixArray_d data_d;
  double newnum;
  double ddf;
  double (*function)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);
  double (*dfunction)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);
  double (*ddfunction)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);

  // HELPER FUNCTIONS
  void NewtonRaphson();
  static double LogNormalDensity
  (double x, double mu, double lambda);


  // NOT IMPLEMENTED!!
  // - to prevent use
  // default constructor
  MetropolisHastings();
  // private copy constructor
  MetropolisHastings(const MetropolisHastings&);
  // private assignment operator
  MetropolisHastings& operator=(const MetropolisHastings&);

public:
  MetropolisHastings(const Vector_d &inparameters,
                     double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
                     double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
                     double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
                     const MatrixArray_i&, const MatrixArray_d&);
   
  ~MetropolisHastings();
  
  int Sample(double*);
  void UpdateParameters(const Vector_d& inparameters);
  void UpdateIntegerData(const MatrixArray_i&);
  void UpdateDoubleData(const MatrixArray_d&);
};

#endif /* !defined ADMIX_METROPOLISHASTINGS_H */
