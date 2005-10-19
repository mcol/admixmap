/** 
 *   ADMIXMAP
 *   chib.cc
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
#include "chib.h"

using namespace std;

#define PR(x) cout << #x << " = " << x << endl;

chib::chib()
{
   LogLikelihood = 0;
   MaxLogPosterior = -9999999;
   MaxLogPrior = -9999999;
}

void chib::setLogLikelihood(double x)
{
   LogLikelihood = x;
}

void chib::addLogPrior(double x)
{
  if(x > MaxLogPrior)
    MaxLogPrior = x;
  LogPrior.push_back(x);
}

void chib::addLogPosteriorObs( double f )
{
   if( f > MaxLogPosterior )
      MaxLogPosterior = f;
   VecLogPosterior.push_back(f);
}

double chib::getLogPosterior()const
{
  //double x = AverageOfLogs( VecLogPosterior, MaxLogPosterior );
  //if( isnan(x) ){
  //  for( unsigned int i = 0; i < VecLogPosterior.size(); i++ )
  //     cout << VecLogPosterior[i] << " ";
  //  cout << endl;
  //  exit(0);
  //}
   return AverageOfLogs( VecLogPosterior, MaxLogPosterior );
}

double chib::getLogPrior()const{
  return AverageOfLogs(LogPrior, MaxLogPrior);
}

double chib::getLogMarginalLikelihood()const{
  double LogPosterior = AverageOfLogs( VecLogPosterior, MaxLogPosterior );
  double logprior = AverageOfLogs( LogPrior, MaxLogPrior );
  return LogLikelihood + logprior - LogPosterior;
}

