// *-*-C++-*-*
/**
 *   Individual.h
 *   header file for Individual class
 */

/*
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */


#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H 1


#include "common.h"
#include "Genome.h"
#include "Chromosome.h"
#include "AlleleFreqs.h"
#include "bclib/LogWriter.h"
#include "GenotypeProbOutputter.h"
#include "bclib/Regression.h"
#include "PedBase.h"
#include "AdmixtureProportions.h"


using namespace std;


/** \addtogroup base
 * @{ */


class AlleleFreqs;


///Class to represent an individual and update individual-level parameters
class Individual : public genepi::PedBase
{
public:
  Individual(unsigned number);
  //Individual(const Options* const options, const InputData* const Data);
  void Initialise(const Options* const options, const InputData* const Data);
  virtual ~Individual();

  void DeleteGenotypes();
  void HMMIsBad(bool loglikisbad);
  static void SetStaticMembers( Genome & pLoci, const Options & options );

  void setOutcome(double*);
  void setCovariates(double*);
  void setGenotypesToMissing();

  const double* getAdmixtureProps(bool isXChrom = false) const;
  const std::vector<hapPair > &getPossibleHapPairs(unsigned int locus)const;
  const int* getSampledHapPair(int locus)const;
  bool GenotypeIsMissing(unsigned int locus)const;///< locus is a comp locus
  bool simpleGenotypeIsMissing(unsigned locus)const;///< locus is a simple locus
  bool isHaploidatLocus(unsigned j)const;
  bool isHaploidIndividual()const;

  virtual double getLogLikelihood(const Options& , bool forceUpdate, bool store);
  virtual double getLogLikelihoodXChr(const Options&, bool forceUpdate,
                                      bool store);
  void storeLogLikelihood(const bool setHMMAsOK); ///< to call if a Metropolis proposal is accepted
  virtual double getLogLikelihoodAtPosteriorMeans(const Options& options);

  //--------------------------------------------------------------------------
  // Public rho-proposal methods, overridden from PedBase, ignored for individuals.
  virtual void setRho( double nv );
  virtual void startRhoProposal();
  virtual void acceptRhoProposal();
  virtual void rejectRhoProposal();
  //--------------------------------------------------------------------------

  void GetLocusAncestry(int locus, int Ancestry[2])const;
  void GetLocusAncestry(int chrm, int locus, int Ancestry[2])const;
  int GetLocusAncestry(int, int, int)const;

  void SampleHiddenStates(const Options& options);

  void SampleHapPair(unsigned chr, unsigned jj, unsigned locus, AlleleFreqs *A,
		    bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts);

  void UpdateAlleleCounts(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool annealthermo)const;

  static Genome & getGenome() { return *Loci; }

  void SampleMissingOutcomes(bclib::DataMatrix *Outcome, const std::vector<bclib::Regression*>& R);

  unsigned int getMyNumber() const { return myNumber  ; } ///< number of this individual, counting from 1
  unsigned int getIndex	  () const { return myNumber-1; } ///< number of this individual, counting from 0
  unsigned int getNumObs  () const { return NumGametes; }

private:
  const unsigned	    myNumber	     ; ///< number of this individual, counting from 1
protected:
  bool			    SexIsFemale	     ;
  bool			    isHaploid	     ; ///< indicates if individual is haploid at all
					       ///< loci or only at X loci
  static unsigned int	    numChromosomes   ;
  static int		    NumHiddenStates  ; ///< Not the number of hidden states in the Hidden Markov Model;
					       ///< apparently, this is actually K, the number of populations:
					       ///< the number of hiddon states for unrelated individuals is
					       ///< typically K^2.
  static Genome *	    Loci	     ;
  static bool		    Xdata	     ; ///< indicates if there is an X chromosome
  static unsigned int	    X_posn	     ; ///< number of X chromosome
  double		    EffectiveL[2]    ; ///< indexed on "g" (0 or 1)
  unsigned		    NumGametes	     ; ///< 1 if assortative mating or haploid data,
					       ///< 2 if random mating and diploid data
  std::vector<unsigned int> gametes	     ; ///< number of gametes on each chromosome
  std::vector<genotype>	    genotypes	     ; ///< indexed on <<something>>
  std::vector<hapPair> *    PossibleHapPairs ; ///< possible haplotype pairs compatible with genotype
					       ///< pointer to array of size Loci->GetNumberOfCompositeLoci()
  bool **		    GenotypesMissing ; ///< indicators for missing genotypes at comp loci
  bool *		    missingGenotypes ; ///< indicators for missing genotypes at simple loci
  std::vector<hapPair>	    sampledHapPairs  ;

  AdmixtureProportions	    Theta	     ; ///< admixture proportions

  int **		    LocusAncestry    ;
  genepi::RhoType	    _rho	     ; ///< sum of intensities
  double *		    Outcome	     ;
  double *		    Covariates	     ;

  struct {
    double value    ; ///< loglikelihood at current parameter values, annealed if coolness < 1.
                      ///< Valid iff 'ready' is true
    double tempvalue; ///< to store values temporarily: holds unnanealed value (-energy),
                      ///< or value at proposed update
    bool   ready    ; ///< true iff value is the loglikelihood at the current parameter values
    bool   HMMisOK  ; ///< true iff values in HMM objects correspond to current parameter values
                      ///< for this individual
  } logLikelihood;

  /// Odds ratios for the X chromosome admixtures
  genepi::cvector<double> psi;

  genepi::cvector<bclib::StepSizeTuner> TunePsiSampler;
  genepi::cvector<double> psistep;
  genepi::cvector<double> SumLogPsi;
  int NumberOfPsiUpdates;

  Individual();
  void SetUniformAdmixtureProps();

  virtual void UpdateHMMInputs(unsigned int j, const Options& options,
                               const AdmixtureProportions& theta,
                               const genepi::RhoType& rho) = 0;

  virtual double getLogLikelihood(const Options& options,
                                  const AdmixtureProportions& theta,
                                  const genepi::RhoType& rho, bool updateHMM);

  virtual double getLogLikelihoodXChr(const Options& options,
                                      const AdmixtureProportions& theta,
                                      const genepi::RhoType& rho, bool updateHMM);

  public:
    const genepi::RhoType & getRho() const;

  void setOddsRatios(const genepi::cvector<double>& psi);

};


/** @} */


#endif /* INDIVIDUAL_H */
