#include "chib.h"

using namespace std;

#define PR(x) cout << #x << " = " << x << endl;

chib::chib()
{
   LogPrior = 0;
   LogLikelihood = 0;
   MaxLogPosterior = -9999999;
}

void chib::setLogLikelihood(double x)
{
   LogLikelihood = x;
}

void chib::setLogPrior(double x)
{
   LogPrior = x;
}

void chib::addLogPosteriorObs( double f )
{
   if( f > MaxLogPosterior )
      MaxLogPosterior = f;
   VecLogPosterior.push_back(f);
}

double chib::getLogPosterior()
{
   double x = AverageOfLogs( VecLogPosterior, MaxLogPosterior );
   if( isnan(x) ){
      for( unsigned int i = 0; i < VecLogPosterior.size(); i++ )
         cout << VecLogPosterior[i] << " ";
      cout << endl;
      exit(0);
   }
   return AverageOfLogs( VecLogPosterior, MaxLogPosterior );
}

void chib::Output(ofstream *outputstream)
{
   double LogPosterior = AverageOfLogs( VecLogPosterior, MaxLogPosterior );
   cout << "Log likelihood:          " << LogLikelihood << endl;
   cout << "Log prior:               " << LogPrior      << endl;
   cout << "Log posterior:           " << LogPosterior  << endl;
   cout << "Log marginal likelihood: " << LogLikelihood + LogPrior - LogPosterior << endl;
   *outputstream << "Log likelihood:          " << LogLikelihood << endl;
   *outputstream << "Log prior:               " << LogPrior      << endl;
   *outputstream << "Log posterior:           " << LogPosterior  << endl;
   *outputstream << "Log marginal likelihood: " << LogLikelihood + LogPrior - LogPosterior << endl;
}
