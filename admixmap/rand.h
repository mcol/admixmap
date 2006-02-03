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
std::vector<int> genmultinomial2( int n, const std::vector<double> p );
//
double MultinomialPDF( const std::vector<int> r, const std::vector<double> theta );
///
///Poisson generator
long ignpoi(double mu);
///
int SampleFromDiscrete( const double probs[] , int numberofelements);
///
void gendirichlet(const size_t K, const double alpha[], double theta[] );
///
void ddigam( double *, double * );
///
void trigam( double *, double * );
