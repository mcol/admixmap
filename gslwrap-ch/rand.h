#include <iostream>
class Vector_d;
///
class Vector_i;
///
double myrand();
///
double myrandRange( double Min, double Max );
///
void smyrand( long seed );
///
double gengam(double aa,double bb);
///
double genbet(double aa,double bb);
///
double gennor(double av,double sd);
//
int genbinomial( int n, double p );
//
unsigned int genpoi( double );
//
Vector_i genmultinomial( int n, Vector_d p );
//
double MultinomialLikelihood( Vector_i r, Vector_d theta );
///
///Poisson generator
long ignpoi(double mu);
///
int SampleFromDiscrete( double probs[] , int numberofelements);
///
Vector_d gendirichlet( Vector_d );
void gendirichlet(const size_t K, const double alpha[], double theta[] );
///
void ddigam( double *, double * );
///
void trigam( double *, double * );
