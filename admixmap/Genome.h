// *-*-C++-*-*
#ifndef GENETIC_ARRAY_H
#define GENETIC_ARRAY_H 1

#include "AbstractCompLocus.h"
#include "LocusVisitor.h"
#include "CompositeLocus.h"
#include "vector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "LogWriter.h"
#include "AdmixOptions.h"
#include "InputData.h"


class Genome : public AbstractCompLocus
// Object is an array of objects of class CompositeLocus
// Used to loop over composite loci on a single chromosome, or an entire genome  
{
private: // members
  Vector Distances;
  double LengthOfGenome;
  double LengthOfXchrm;
  int NumberOfCompositeLoci;
  bool X_data;
  AbstractCompLocus **TheArray;
  std::vector< std::vector< int > > _chrmandlocus;

  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);

public:

  Genome();
  Genome(int);
  virtual ~Genome();

  std::vector< int > GetChrmAndLocus( int );

   bool isX_data();

   virtual void
  accept(LocusVisitor&);

  // composite-level methods
  virtual AbstractCompLocus*&
  operator()(int) const;

  void SetLabels(const std::vector<std::string> &labels, Vector_d temp);

  virtual Vector
  GetDistances();

  virtual void
  SetNumberOfCompositeLoci(int);

  virtual int
  GetNumberOfCompositeLoci();

  virtual float
  GetDistance(int);

  virtual void
  SetDistance(int,float);

  virtual Genome*
  GetChromosomes( int, std::vector<std::string> );

  virtual int
  size();

  // individual-level methods
  virtual void
  AddLocus( int );
  
  virtual Matrix_d
  GetAlleleFreqs();
  
  virtual Vector_d
  GetFst();
  
  virtual Matrix_d
  GetInfo();
  
  virtual std::string
  GetLabel(int);
  
  virtual Vector_i
  GetHapLabels(int);
  
  virtual Matrix_i
  GetLikelihoodAlleleFreqs();
  
  virtual int
  GetMergedHaplotype( int i );
  
  virtual int
  GetNumberOfAllelesOfLocus( int );
  
  virtual int
  GetNumberOfLoci();
  
  virtual int
  GetNumberOfMergedHaplotypes();
  
  virtual int
  GetNumberOfStates();
  
  virtual Vector_d
  GetPriorAlleleFreqs( int );
  
  virtual Matrix_d
  GetScore();
  
  virtual Matrix_d
  GetScoreSq();
  
  virtual int
  GetSize();

  virtual Vector_d
  GetStatsForEta( int );
  
  virtual Matrix_d
  GetSumAlleleFreqs();
  
  virtual int
  IsRandom();
  
  virtual void
  ResetLikelihoodAlleleFreqs();
  
  virtual void
  ResetScoreForMisSpecOfAlleleFreqs();

  virtual void
  ResetSumAlleleFreqs();
  
  virtual void
  SampleAlleleFreqs( int );
  
  virtual void
  SamplePriorAlleleFreqs( Vector_d eta );
  
  virtual void
  SetAlleleFreqs( Matrix_d );
  
  virtual void
  SetDefaultAlleleFreqs( int );
  
  virtual void
  SetDefaultMergeHaplotypes( Vector_d alpha );
  
  virtual void
  SetLabel( int, std::string );
  
  virtual void
  SetHistoricalAlleleFreqs( Matrix_d );
  
  virtual void
  SetNumberOfAllelesOfLocus( int, int );
  
  virtual void
  SetNumberOfLoci( int );
  
  virtual void
  SumScoreForMisSpecOfAlleleFreqs();
  
  virtual void
  UpdateFst();
  
  virtual void
  UpdatePriorAlleleFreqs( int, const Vector_d& ){};
  
  virtual void
  UpdatePriorAlleleFreqsGlobal( int, const std::vector<Vector_d>& );
  
   double GetLengthOfGenome();
   double GetLengthOfXchrm();
};

#endif /* !GENETIC_ARRAY_H */
