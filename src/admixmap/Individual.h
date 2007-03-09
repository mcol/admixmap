// *-*-C++-*-*
/** 
 *   Individual.h 
 *   header file for Individual class
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
#include "utils/LogWriter.h"
#include "GenotypeProbOutputter.h"
#include "regression/Regression.h"

#ifdef HAVE_CPPUNIT
#include <cppunit/Exception.h>
#endif

using namespace std;

class AlleleFreqs;
class Regression;

///Class to represent an individual and update individual-level parameters
class Individual
{
public:
  Individual();
  Individual(int number, const Options* const options, const InputData* const Data);
  void Initialise(int number, const Options* const options, const InputData* const Data);
  virtual ~Individual();

  void DeleteGenotypes();
  void HMMIsBad(bool loglikisbad);
  static void SetStaticMembers(Genome* const pLoci, const Options* const options);
  static void setGenome(Genome * const);
  static void setPopulations(const int);
  
  static const int getNumberOfHiddenStates();
  
  void setOutcome(double*);
  void setCovariates(double*);
  void setGenotypesToMissing();
  virtual void SetMissingGenotypes() = 0;

  const double* getAdmixtureProps()const;
  const vector<hapPair > &getPossibleHapPairs(unsigned int locus)const;
  const int* getSampledHapPair(int locus)const;
  bool GenotypeIsMissing(unsigned int locus)const;//< locus is a comp locus
  bool simpleGenotypeIsMissing(unsigned locus)const;//< locus is a simple locus
  bool isHaploidatLocus(unsigned j)const;
  bool isHaploidIndividual()const;

  virtual double getLogLikelihood(const Options* const , const bool forceUpdate, const bool store);
  void storeLogLikelihood(const bool setHMMAsOK); //< to call if a Metropolis proposal is accepted
  virtual double getLogLikelihoodAtPosteriorMeans(const Options* const options);

  void GetLocusAncestry(int locus, int Ancestry[2])const;
  void GetLocusAncestry(int chrm, int locus, int Ancestry[2])const;
  int GetLocusAncestry(int, int, int)const;
   
  void SampleLocusAncestry(const Options* const options);
  void AccumulateAncestry(int* SumAncestry);
#ifdef PARALLEL
  void SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts,
		     const double* const AlleleProbs);
#else
  void SampleHapPair(unsigned chr, unsigned jj, unsigned locus, AlleleFreqs *A, bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts);
#endif
  void UpdateAlleleCounts(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool annealthermo)const;

  void SampleMissingOutcomes(DataMatrix *Outcome, const vector<Regression*>& R);
  vector<vector<double> >& getUnorderedProbs(const unsigned int);
  void calculateUnorderedGenotypeProbs(void);
  vector<double> getStateProbs(const bool, int, int) const;
  void calculateUnorderedGenotypeProbs(unsigned);
protected:
  unsigned myNumber;//< number of this individual, counting from 1
  bool SexIsFemale;
  bool isHaploid;//< indicates if individual is haploid at all loci or only at X loci
  static unsigned int numChromosomes;
  static int Populations; //< Number of hidden states in the hidden Markov model
  static Genome *Loci;
  static bool Xdata;//< indicates if there is an X chromosome
  static unsigned int X_posn;  //number of X chromosome
  double EffectiveL[2];
  unsigned NumGametes; //< 1 if assortative mating or haploid data, 2 if random mating and diploid data
  std::vector< unsigned int > gametes;//< number of gametes on each chromosome
  std::vector<genotype> genotypes;
  std::vector<hapPair> *PossibleHapPairs;//possible haplotype pairs compatible with genotype
  bool **GenotypesMissing;//< indicators for missing genotypes at comp loci
  bool *missingGenotypes;//< indicators for missing genotypes at simple loci
  vector<hapPair> sampledHapPairs;

  double *Theta;//< admixture proportions

  int **LocusAncestry;
  std::vector<double> _rho;//< sum of intensities
  double* Outcome;
  double* Covariates;
  vector<vector<vector<double> > > UnorderedProbs;

  struct {
    double value; //< loglikelihood at current parameter values, annealed if coolness < 1.  Valid iff 'ready' is true
    double tempvalue; //< to store values temporarily: holds unnanealed value (-energy), or value at proposed update   
    bool ready;//< true iff value is the loglikelihood at the current parameter values
    bool HMMisOK;//< true iff values in HMM objects correspond to current parameter values for this individual
  } logLikelihood;
  
  void SetUniformAdmixtureProps();

  virtual void UpdateHMMInputs(unsigned int j, const Options* const options, 
			     const double* const theta, const vector<double> rho) = 0;
  virtual double getLogLikelihood(const Options* const options, 
			  const double* const theta, const vector<double > rho, bool updateHMM);
};

#endif /* INDIVIDUAL_H */

