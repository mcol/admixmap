// *-*-C++-*-*
#ifndef ABSTRACT_COMP_LOCUS_H
#define ABSTRACT_COMP_LOCUS_H 1

#include <string>
#include "LocusVisitor.h"
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "vector_i.h"
#include "matrix_i.h"
#include "AdmixOptions.h"
#include "MatrixArray_i.h"
#include <vector>
#include <assert.h>

class Individual;
class AbstractCompLocus
{
public:
  virtual ~AbstractCompLocus() {}

  virtual Matrix_d
  getExpectedAncestry( int )
      {
         assert(2+2==5);
         Matrix_d null;
         return null;
      }
  
   virtual double
   getLogLikelihood(){return 2;}
   
   virtual Matrix_i
  SampleForLocusAncestry(Individual*)
  {
    assert(2+2==5);
    Matrix_i null;
    return null;
  }
  
   virtual void
  UpdateParameters(Individual*,Matrix_d&,AdmixOptions*,std::vector<Vector_d>&,bool){}
  
   virtual void
  UpdateParametersHaploid(Individual*,Matrix_d&,AdmixOptions*,
                          std::vector<Vector_d>&,bool){}
  
  virtual Vector_i
  SampleForHaploidLocusAncestry(Individual*)
  {
    assert(2+2==5);
    Vector_i null;
    return null;
  }
  
  virtual void
  accept(LocusVisitor&) = 0;

  virtual void
  AddLocus( int ) = 0;
  
  virtual Vector_i
  GetAlleleCountsInHaplotype(const std::vector<unsigned int>&)
  {
    assert(2+2==5);
    Vector_i null_vector_i;
    return null_vector_i;
  }
  
  virtual Matrix_d
  GetAlleleFreqs() = 0;
  
  virtual Vector_d
  GetFst() = 0;
  
  virtual Matrix_d
  GetInfo() = 0;
  
  virtual std::string
  GetLabel(int) = 0;
  
  virtual Vector_i
  GetHapLabels(int) = 0;
  
  virtual Matrix_d
  GetLikelihood( std::vector<unsigned int>, bool, bool )
  {
    assert(2+2==5);
    Matrix_d null_matrix_d;
    return null_matrix_d;
  }
  
  virtual Matrix_i
  GetLikelihoodAlleleFreqs() = 0;
  
//   virtual Vector_i
//   GetLikelihoodAlleleFreqs(int) = 0;
  
  virtual int
  GetLocus(int i) { return i; }

  virtual int
  GetMergedHaplotype( int  ) = 0;
  
  virtual int
  GetNumberOfAllelesOfLocus( int ) = 0;
  
  virtual Vector_i
  GetNumberOfAlleles() { Vector_i null_vector(1); return null_vector; }
  
  virtual int
  GetNumberOfLoci() = 0;
  
  virtual int
  GetNumberOfMergedHaplotypes() = 0;
  
  virtual int
  GetNumberOfStates() = 0;
  
  virtual Vector_d
  GetPriorAlleleFreqs( int ) = 0;
  
  virtual Matrix_d
  GetScore() = 0;
  
  virtual Matrix_d
  GetScoreSq() = 0;
  
  virtual int
  GetSize() = 0;

  virtual Vector_d
  GetStatsForEta( int ) = 0;
  
  virtual Matrix_d
  GetSumAlleleFreqs() = 0;
  
  virtual int
  IsRandom() = 0;
  
  virtual void
  ResetLikelihoodAlleleFreqs() = 0;
  
  virtual void
  ResetScoreForMisSpecOfAlleleFreqs() = 0;

  virtual void
  ResetSumAlleleFreqs() = 0;
  
  virtual void
  SampleAlleleFreqs( int ) = 0;
  
  virtual Vector_i
  SampleHaplotype(const std::vector<unsigned int>&, Vector_i )
  {
    assert(2+2==5);
    Vector_i null_vector_i;
    return null_vector_i;
  }
  
  virtual void
  SamplePriorAlleleFreqs( Vector_d eta ) = 0;
  
  virtual void
  SetAlleleFreqs( Matrix_d ) = 0;
  
  virtual void
  SetDefaultAlleleFreqs( int ) = 0;
  
  virtual void
  SetDefaultMergeHaplotypes( Vector_d alpha ) = 0;
  
  virtual void
  SetLabel( int, std::string ) = 0;
  
  virtual void
  SetHistoricalAlleleFreqs( Matrix_d ) = 0;
  
  virtual void
  SetNumberOfAllelesOfLocus( int, int ) = 0;
  
  virtual void
  SetNumberOfLoci( int ) = 0;

  virtual void
  SetPriorAlleleFreqs(Matrix_d,bool) 
  {
    assert(2+2==5);
  }

  virtual void
  SumScoreForMisSpecOfAlleleFreqs() = 0;
  
  virtual void
  UpdateFst() = 0;
  
  virtual void
  UpdateLikelihoodAlleleFreqs( const std::vector<unsigned int>&, Vector_i )
  {
    assert(2+2==5);
  }
  
  virtual void
  UpdateLikelihoodAlleleFreqs_HaploidData( const std::vector<unsigned int>&, int )
  {
    assert(2+2==5);
  }
  
   virtual void
   UpdatePriorAlleleFreqs( int, const Vector_d& ) = 0;
  
   virtual void
   UpdatePriorAlleleFreqsGlobal( int, const std::vector<Vector_d>& ) = 0;
  
  virtual void
  UpdateScoreForMisSpecOfAlleleFreqs( Matrix_d, std::vector<unsigned int>)
  {
    assert(2+2==5);
  }
  
};

#endif /* !ABSTRACT_COMP_LOCUS_H */
