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
// multinomial probability of counts r given probability vector theta
double MultinomialLikelihood( Vector_i r, Vector_d theta );
///
///Poisson generator
long ignpoi(double mu);
///
// sampler with cumulative probs as argument
int SampleFromDiscrete( Vector_d *cdf );
///
// sampler with probs as argument
int SampleFromDiscrete2( Vector_d *probs );
///
/// sample proportion vector from Dirichlet distribution 
// Dirichlet parameter vector may contain zeroes
Vector_d gendirichlet( Vector_d );
///
// these next two functions don't really belong in the rand.h file
// 
// digamma function: first parameter is function argument, second parameter is result 
void ddigam( double *, double * );
///
// trigamma function
void trigam( double *, double * );
