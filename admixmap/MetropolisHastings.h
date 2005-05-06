// *-*-C++-*-*
#ifndef ADMIX_METROPOLISHASTINGS_H
#define ADMIX_METROPOLISHASTINGS_H 1

#include "rand.h"
#include "vector_d.h"
#include "matrix_i.h"
#include "matrix_d.h"

class MetropolisHastings
{
private: // members
  Vector_d parameters;
  Matrix_i data_i;
  Matrix_d data_d;
  double newnum;
  double ddf;
  double (*function)(Vector_d &, Matrix_i &, Matrix_d &, double xin);
  double (*dfunction)(Vector_d &, Matrix_i &, Matrix_d &, double xin);
  double (*ddfunction)(Vector_d &, Matrix_i &, Matrix_d &, double xin);

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
                     double (*funct)(Vector_d&, Matrix_i&, Matrix_d&, double),
                     double (*dfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
                     double (*ddfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
                     const Matrix_i&, const Matrix_d&);
   
  ~MetropolisHastings();
  
  int Sample(double*);
  void UpdateParameters(const Vector_d& inparameters);
  void UpdateIntegerData(const Matrix_i&);
  void UpdateDoubleData(const Matrix_d&);
};

#endif /* !defined ADMIX_METROPOLISHASTINGS_H */
