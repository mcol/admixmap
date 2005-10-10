// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   chib.h
 *   this is meant to be a generic class for the Chib algorithm 
 *   calculating the marginal likelihood of the model from MCMC output 
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef CHIB_H
#define CHIB_H 1

#include <iostream>
#include "functions.h"

class chib

{
public:
  chib();
  void setLogLikelihood( double );
  void addLogPrior( double );
  void addLogPosteriorObs( double );
  double getLogPosterior();
  double getMarginalLikelihood();
  double getLogMarginalLikelihood();
  double getLogPrior();
  double getLogLikelihood(){
    return LogLikelihood;
  };

private:
  double LogLikelihood;
  std::vector<double> LogPrior;
  std::vector<double> VecLogPosterior;
  double MaxLogPosterior;
  double MaxLogPrior;
  
};

#endif /* !defined CHIB_H */
