// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   chib.h
 *   a generic class for the Chib algorithm 
 *   calculating the marginal likelihood of the model from MCMC output 
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef CHIB_H
#define CHIB_H 1

#include "utils/LogWriter.h"
#include <vector>

class chib
{
public:
  chib();
  void Reset();
  void setLogLikelihood( const double x );
  void setLogPrior( const double x );
  void addLogPosteriorObs( const double f );

  double getLogPrior()const;
  double getLogLikelihood()const;
  double getLogPosterior()const;
  double getLogMarginalLikelihood()const;
  void outputResults(LogWriter &Log) const;

//   double getLogLikelihood()const{
//     return LogLikelihood;
//   };

private:
  double LogLikelihood;
  double LogPrior;
  std::vector<double> VecLogPosterior;
  double MaxLogPosterior;
  
};

#endif /* !defined CHIB_H */
