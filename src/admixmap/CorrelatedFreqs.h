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
/// \file CorrelatedFreqs.h
/// Definition of the CorrelatedFreqs class.
//=============================================================================

#ifndef CORRELATEDFREQS_H
#define CORRELATEDFREQS_H

#include "AdmixFreqs.h"
#include "bclib/MuSampler.h"
#include "bclib/DispersionSampler.h"

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN


/** \addtogroup admixmap
 * @{ */


///class to implement a correlated allele frequency model
class CorrelatedFreqs : public AdmixFreqs{
public:
  CorrelatedFreqs();
  void Initialise(Options* const a_options, InputData* const a_data, 
		  Genome *pLoci, bclib::LogWriter &Log, bool MAP);
  ///print parameters of prior on frequencies and their Dirichlet parameters
  void PrintPrior(const Vector_s& PopLabels, bclib::LogWriter& Log)const;
  void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);
  ///write ergodic average of dispersion to file
  void OutputErgodicAvg( int samples, std::ofstream *avgstream)const;
  ///write dispersion parameter to file
  void OutputParams();
  ///write dispersion parameter to stream
  void OutputParams(bclib::Delimitedostream& os)const;
  void PrintAcceptanceRates(bclib::LogWriter& Log);

private:
  double eta; ///<dispersion parameter
  double SumEta;
  double psi,tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
  bclib::MuSampler *muSampler;
  //new sampler for eta
  bclib::DispersionSampler EtaSampler;
 
  ///adaptive RW sampler for eta
  bclib::StepSizeTuner TuneEtaSampler;
  int NumberOfEtaUpdates;
  int NumberAccepted; 
  double etastep;///< step size in random walk sampler
  double etastep0;///< initial step size
  int w;//The eta sampler is tuned every w updates.

  // sampler for Dirichlet proportion parameters 
  std::vector<bclib::StepSizeTuner> MuProposal;
  
  bclib::DelimitedFileWriter outputstream;//outputs eta to paramfile

  void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  void InitializeEtaOutputFile(const char* filename, const Vector_s& PopulationLabels, bclib::LogWriter &Log);
  void SetDefaultPriorParams(int i, double defaultpriorparams);
  void SampleDirichletParams();
  void SampleEtaWithRandomWalk(bool updateSumEta);
  const double* GetStatsForEta( int locus )const;
  void UpdatePriorAlleleFreqs(const vector<vector<double> >& mu);
  void SampleAlleleFreqs(int i, double coolness);
  float getEtaSamplerAcceptanceRate()const;
  float getEtaSamplerStepsize()const;
};


/** @} */


#endif
