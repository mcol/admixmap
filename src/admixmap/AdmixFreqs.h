// *-*-C++-*-*
/** 
 *   AdmixFreqs.h
 *   header file for AdmixFreqs class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ADMIXFREQS_H
#define ADMIXFREQS_H 1

#include "AlleleFreqs.h"
#include "bclib/MuSampler.h"
#include "bclib/StepSizeTuner.h"
#include "bclib/AdaptiveRejection.h"
#include "bclib/DelimitedFileWriter.h"

class AdmixOptions;
class InputAdmixData;

/// Class to hold allele/haplotype frequencies and their priors in a dispersion model.
class AdmixFreqs : public AlleleFreqs{

public:
  AdmixFreqs();
  virtual ~AdmixFreqs();

  virtual void Initialise(Options* const options, InputData* const Data, Genome *pLoci, bclib::LogWriter &Log, bool MAP=false);
  virtual void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);
  virtual void PrintPrior(const Vector_s&, bclib::LogWriter& Log)const;
  virtual void OutputErgodicAvg( int , std::ofstream *)const;
  virtual void PrintAcceptanceRates(bclib::LogWriter& Log)const;
  virtual void OutputParams();
  virtual void OutputParams(bclib::Delimitedostream& os)const;

protected:
  virtual void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  virtual void LoadAlleleFreqs(Options* const options, InputData* const data, bclib::LogWriter &Log);
  void SetDefaultPriorParams(int i, double defaultpriorparams);
  virtual void SampleAlleleFreqs(int, const double coolness);

  static double convertValueFromFile(const string s);
};

#endif
