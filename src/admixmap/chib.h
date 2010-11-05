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
/// \file chib.h
/// Definition of the chib class.
//=============================================================================

#ifndef CHIB_H
#define CHIB_H 1

#include <vector>

namespace bclib{
  class LogWriter;
}

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
  void outputResults(bclib::LogWriter &Log) const;

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
