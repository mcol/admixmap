//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AdmixOptions.h
/// Definition of the AdmixOptions class.
//=============================================================================

#ifndef ADMIX_OPTIONS_H
#define ADMIX_OPTIONS_H 1

#include "Options.h"
using namespace::std;

#include "bclib/cvector.h"
#include <string>


/** \addtogroup admixmap
 * @{ */


/// Class to hold program options
class AdmixOptions : public Options
{
private:
  bool noConjugateUpdate;
  bool printPedSummary;
  int  excludePedsOver;
  bool singleIndividualAutoHier; ///< Force hierarchical model for single individual

public:
  AdmixOptions(int, char**);
  ~AdmixOptions();

  bool checkOptions(bclib::LogWriter& Log, int NumberOfIndividuals);
  void PrintUserOptions(const char* filename);
  bool SetOptions();

  //output files
  bool outputParams()const;
  const char *getEtaOutputFilename() const;
  const char *getIndAdmixtureFilename() const;
  const char* getIndAdmixModeFilename()const;

  //input file names
  const char *getAlleleFreqFilename() const;
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getReportedAncestryFilename() const;
  const char *getEtaPriorFilename() const;

  double getRhoalpha() const;
  double getRhobeta() const;
  double getRhobetaShape()const;
  double getRhobetaRate()const;
  double getRhoPriorMean()const;

  double getLogPsiPriorMean() const;
  double getLogPsiPriorPrec() const;
  double getLogPsigammaShape() const;
  double getLogPsigammaRate() const;

  const genepi::cvector<double> & getInitAlpha(int) const;
  const genepi::cvector<genepi::cvector<double> > & getInitAlpha() const;
  int sizeInitAlpha() const;
  double getEtaMean() const;
  double getEtaVar() const;
  const genepi::cvector<float>& getPopAdmixSamplerParams() const;

  //indicators and model options
  bool getCorrelatedAlleleFreqs() const;
  bool isRandomMatingModel() const;
  bool isGlobalRho() const;
  bool isGlobalPsi() const;
  bool PopAdmixturePropsAreEqual()const;
  bool getIndAdmixHierIndicator() const;
  bool getHapMixModelIndicator() const;
  bool getChibIndicator() const;
  bool getTestOneIndivIndicator() const;
  bool isAdmixed(unsigned) const;
  bool isSymmetric()const;
  int getPopulations() const;
  void setPopulations(int num);
  bool getFixedAlleleFreqs() const;
  bool getLocusAncestryProbsIndicator() const;
  const string& getPopLabelString()const;

  bool getNoConjugateUpdate() const { return noConjugateUpdate; }
  bool getPrintPedSummary  () const { return printPedSummary  ; }
  int  getExcludePedsOver  () const { return excludePedsOver  ; }

  /// Force hierarchical model for single individual.
  bool getSingleIndividualAutoHier() const { return singleIndividualAutoHier; }

  //score test indicators 
  bool getScoreTestIndicator() const;
  bool getTestForAdmixtureAssociation() const;
  bool getTestForAffectedsOnly() const;
  void setTestForAffectedsOnly(bool);  
  bool getTestForHaplotypeAssociation() const;
  void setTestForHaplotypeAssociation(bool);
  bool getTestForLinkageWithAncestry() const;
  void setTestForLinkageWithAncestry(bool);
  bool getTestForMisspecifiedAlleleFreqs() const;
  bool getTestForMisspecifiedAlleleFreqs2() const;
  
  //other test indicators
  bool getTestForDispersion() const;
  bool getStratificationTest() const;
  void setStratificationTest(bool);
  bool getOutputFST() const;
  bool getLocusForTestIndicator() const;
  int getLocusForTest() const;
  

  /// Aggregate: are any association tests requested?
  bool hasAnyAssociationTests() const
    {
    return (getTestForHaplotypeAssociation  () ||
	    getTestForAdmixtureAssociation  () ||
	    getTestForAllelicAssociation    () ||
	    getTestForResidualAllelicAssoc  () );
    }

private:
  int Populations;
  bool locusForTestIndicator;
  int LocusForTest;
  bool correlatedallelefreqs;
  bool RandomMatingModel;//random mating model
  bool GlobalRho;//indicator for global rho
  bool GlobalPsi; // indicator for population-level psi
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  static const bool HapMixModelIndicator = false; //model haplotypes with mixture model
  bool chibIndicator;//calculate marginal likelihood using Chib method
  bool TestOneIndivIndicator;//calculate marginal likelihood for one individual only
  bool PopAdmixPropsAreEqual;
  bool LocusAncestryProbsIndicator;

  bool TestForAdmixtureAssociation;
  bool StratificationTestIndicator;
  bool OutputFST;
  bool TestForAffectedsOnly;
  bool TestForHaplotypeAssociation;
  bool TestForDispersion;
  bool TestForLinkageWithAncestry;
  bool TestForMisspecifiedAlleleFreqs;
  bool TestForMisspecifiedAlleleFreqs2;
  bool ScoreTestIndicator; //indicator for any of the score tests in ScoreTests class

  std::vector<bool> _admixed;
  bool _symmetric;

  //priors
  //double Rhoalpha, Rhobeta;//gamma parameters for sumintensities
  //double RhobetaShape, RhobetaRate;//gamma parameters for prior on rhobeta
  genepi::cvector<double> globalrhoPrior;
  genepi::cvector<double> rhoPrior;
  genepi::cvector<double> alpha0;
  genepi::cvector<double> alpha1;
  genepi::cvector< genepi::cvector<double> > initalpha;
  double etamean, etavar;//gamma parameters for dispersion parameter

  // prior parameters on log psi (mean and precision; also gamma parameters
  // for the prior on precision if globalpsi=0)
  genepi::cvector<double> logpsiPrior;

  genepi::cvector<float> popAdmixSamplerParams;//parameters for sampler of population admixture

  //string AlleleFreqPriorOutputFilename;
  string EtaOutputFilename;
  string IndAdmixtureFilename;
  string IndAdmixModeFilename;
  string alleleFreqFilename;
  string HistoricalAlleleFreqFilename;
  string EtaPriorFilename;
  string ReportedAncestryFilename;
  string PopLabels;

  void SetDefaultValues();
  void DefineOptions();
  void setInitAlpha(bclib::LogWriter &Log);
  bool CheckInitAlpha( const genepi::cvector<double> & alphatemp )const;
  void AddFilenamesToUserOptions();

  // UNIMPLEMENTED: to avoid use
  AdmixOptions();
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
};


/** @} */


#endif /* ADMIX_OPTIONS_H */
