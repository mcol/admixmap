// *-*-C++-*-*
#ifndef CHIB_H
#define CHIB_H 1

#include <iostream>
#include "functions.h"

class chib
// this is meant to be a generic class for the Chib algorithm 
// calculating the marginal likelihood of the model from MCMC output 
{
  private:
   double LogLikelihood;
   double LogPrior;
   std::vector<double> VecLogPosterior;
   double MaxLogPosterior;

  public:
   chib();
   void setLogLikelihood( double );
   void setLogPrior( double );
   void addLogPosteriorObs( double );
   double getLogPosterior();
   double getMarginalLikelihood();
   void Output( std::ofstream * );
};

#endif /* !defined CHIB_H */
