// *-*-C++-*-*
#ifndef ADMIX_METROPOLISHASTINGS_H
#define ADMIX_METROPOLISHASTINGS_H 1

#include "rand.h"
#include "matrix_i.h"
#include "matrix_d.h"

class MetropolisHastings
{
public:
  MetropolisHastings(const double* inparameters,
                     double (*funct)(const double*, Matrix_i&, Matrix_d&, double),
                     double (*dfunct)(const double*, Matrix_i&, Matrix_d&, double),
                     double (*ddfunct)(const double*, Matrix_i&, Matrix_d&, double),
                     const Matrix_i&, const Matrix_d&);
   
  ~MetropolisHastings();
  
  int Sample(double*);
  void UpdateParameters(const double* inparameters);
  void UpdateIntegerData(const Matrix_i&);
  void UpdateDoubleData(const Matrix_d&);

private: // members
  unsigned dim;
  const double *parameters;
  Matrix_i data_i;
  Matrix_d data_d;
  double newnum;
  double ddf;
  double (*function)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*dfunction)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*ddfunction)(const double*, Matrix_i &, Matrix_d &, double xin);

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

};

#endif /* !defined ADMIX_METROPOLISHASTINGS_H */
