#include <iostream>
#include <vector>

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
std::vector<int> genmultinomial2( int n, std::vector<double> p );
//
double MultinomialLikelihood( std::vector<int> r, std::vector<double> theta );
///
///Poisson generator
long ignpoi(double mu);
///
int SampleFromDiscrete( double probs[] , int numberofelements);
///
void gendirichlet(const size_t K, const double alpha[], double theta[] );
///
void ddigam( double *, double * );
///
void trigam( double *, double * );
