//=============================================================================
//
// Copyright (C) 2002-2006  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file chib.cc
/// Implementation of the chib class.
//=============================================================================

#include "chib.h"
#include "bclib/misc.h"
#include "bclib/LogWriter.h"
//#include <limits>

using namespace std;

#define PR(x) cout << #x << " = " << x << endl;

chib::chib()
{
   LogLikelihood = 0.0;
   LogPrior = 0.0;
}

void chib::Reset(){
  VecLogPosterior.clear();
  LogLikelihood = 0.0;
  LogPrior = 0.0;
}

void chib::setLogLikelihood(const double x)
{
   LogLikelihood = x;
}

void chib::setLogPrior(const double x)
{
  LogPrior = x;
}

void chib::addLogPosteriorObs( const double f )
{
   VecLogPosterior.push_back(f);
   if(VecLogPosterior.size() == 1) MaxLogPosterior = f;
   if( f > MaxLogPosterior ) MaxLogPosterior = f;
}

double chib::getLogLikelihood()const{
  return LogLikelihood; 
}

double chib::getLogPrior()const{
  return LogPrior; 
}

double chib::getLogPosterior()const
{
  return bclib::AverageOfLogs( VecLogPosterior, MaxLogPosterior );
}

double chib::getLogMarginalLikelihood()const{
  double LogPosterior = bclib::AverageOfLogs( VecLogPosterior, MaxLogPosterior );
  return LogLikelihood + LogPrior - LogPosterior;
}

void chib::outputResults(bclib::LogWriter & Log)const {
  Log.setDisplayMode(bclib::On);
  Log << "\nCalculation of Chib algorithm at fixed parameter values"
      << "\nDeviance\t" << -2.0*getLogLikelihood()
      << "\nLogLikelihood\t" << getLogLikelihood()
      << "\nLogPrior\t" << getLogPrior()
      << "\nLogPosterior\t" << getLogPosterior()
      << "\nLogMarginalLikelihoodFromChibAlgorithm\t" << getLogMarginalLikelihood()
      << "\n";
}
