// *-*-C++-*-*
/* 
 *   Annealer.h 
 *   Classs to implement simulated annealing and thermodynamic integration
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef THERMO_H
#define THERMO_H 1

#include "bclib/LogWriter.h"
class Annealer{

public:
  Annealer();
  void Initialise(bool thermo, unsigned numAnnealedRuns, unsigned samples, unsigned burnin, const char* filename);
  ~Annealer();

  void PrintRunLengths(LogWriter& Log, bool testoneindiv);
  void SetAnnealingSchedule();
  bool SetRunLengths(int run, unsigned* samples, unsigned* burnin, double* coolness);
  void CalculateLogEvidence(int run, double coolness, double SumEnergy, double SumEnergySq, unsigned samples);
  void CalculateLogEvidence(double *SumEnergy, double*SumEnergySq, unsigned size);
  void PrintResults(LogWriter& Log, double D_hat);
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


#endif /* THERMO_H */
