#ifndef ICOMPOSITELOCUS_H_
#define ICOMPOSITELOCUS_H_

#include <string>
#include <vector>
#include "../Haplotype.h"

using std::vector;

class hapPair;

class ICompositeLocus
{
public:
//	ICompositeLocus();
	virtual ~ICompositeLocus();

  /* virtual */ static void SetRandomAlleleFreqs(bool) {
    throw string("static void ICompositeLocus::SetRandomAlleleFreqs(bool) is unimplemented");
  };
  /* virtual */ static void SetNumberOfPopulations( int ) {
    throw string("static void ICompositeLocus::SetNumberOfPopulations( int ) is unimplemented");
  };

  virtual void AccumulateAlleleProbs() = 0;
  virtual void AddLocus( int, std::string ) = 0;
  virtual const std::string GetLabel(int)const = 0;
  virtual int GetNumberOfAllelesOfLocus( int )const = 0;
  virtual int GetNumberOfStates()const = 0;
  virtual int GetNumberOfLoci()const = 0;
  virtual void InitialiseHapPairProbs(const double* const allelefreqs) = 0;
  virtual void setAlleleProbsMAP(const double* const Freqs) = 0;
  virtual void InitialiseHapPairProbsMAP() = 0;
#ifndef PARALLEL
  // Those methods are not implemented in the parallel code.
  // We can't leave the prototypes here for the parallel version,
  // because linking fails if we do.
  virtual void SetHapPairProbs() = 0;
  virtual void SampleHapPair(hapPair*, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const = 0;
#endif
  virtual void getConditionalHapPairProbs(std::vector<double>& Probs, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const = 0;
  virtual const std::vector<int> getAlleleCounts(int a, const int* happair)const = 0;
  virtual const int *GetHapLabels( int ) const = 0;
  virtual const vector<int> getHaplotypeCounts(const int* happair) = 0;
  virtual void getLocusAlleleProbs(double **P, int k)const = 0;
  virtual void SetDefaultMergeHaplotypes( const double* const alpha) = 0;
  virtual int GetNumberOfMergedHaplotypes()const = 0;
  virtual void SetHapPairProbsToPosteriorMeans(int iterations) = 0;
  virtual void GetGenotypeProbs(double *Probs, const std::vector<hapPair> &HaplotypePairs, bool chibindicator)const = 0;
  virtual void GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
  Haplotype HaplotypeSetter;
};

#endif /*ICOMPOSITELOCUS_H_*/
