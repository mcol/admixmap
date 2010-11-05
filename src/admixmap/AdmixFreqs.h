//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AdmixFreqs.h
/// Definition of the AdmixFreqs class.
//=============================================================================

#ifndef ADMIXFREQS_H
#define ADMIXFREQS_H 1

#include "AlleleFreqs.h"
#include "bclib/DelimitedFileWriter.h"

class AdmixOptions;
class InputAdmixData;


/** \addtogroup admixmap
 * @{ */


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


/** @} */


#endif
