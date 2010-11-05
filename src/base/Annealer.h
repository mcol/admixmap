//=============================================================================
//
// Copyright (C) 2007  David O'Donnell and Paul McKeigue
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
/// \file Annealer.h
/// Definition of the Annealer class.
//=============================================================================

#ifndef THERMO_H
#define THERMO_H 1


#include <fstream>

namespace bclib {
  class LogWriter;
}

/** \addtogroup base
 * @{ */


class Annealer{

public:
  Annealer();
  void Initialise(bool thermo, unsigned numAnnealedRuns, unsigned samples, unsigned burnin, const char* filename);
  ~Annealer();

  void PrintRunLengths(bclib::LogWriter& Log, bool testoneindiv);
  void SetAnnealingSchedule();
  bool SetRunLengths(int run, unsigned* samples, unsigned* burnin, double* coolness);
  void CalculateLogEvidence(int run, double coolness, double SumEnergy, double SumEnergySq, unsigned samples);
  void CalculateLogEvidence(double *SumEnergy, double*SumEnergySq, unsigned size);
  void PrintResults(bclib::LogWriter& Log, double D_hat);
  const double* GetCoolnesses()const;
private:
  bool Thermo;
  unsigned _samples, _burnin;
  double IntervalRatio;
  int NumAnnealedRuns;
  double SumEnergy, SumEnergySq, LogEvidence;
  double MeanEnergy, VarEnergy;
  double LastMeanEnergy;
  double *IntervalWidths;
  double *Coolnesses; 

  std::ofstream annealstream;//for monitoring energy when annealing
};


/** @} */


#endif /* THERMO_H */
