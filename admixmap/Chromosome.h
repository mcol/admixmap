// *-*-C++-*-*
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "Genome.h"
#include "LocusVisitor.h"
#include "AdmixOptions.h"
#include "matrix_d.h"
#include "matrix_i.h"
#include "MatrixArray_d.h"
#include "MatrixArray_i.h"
#include "vector_d.h"
#include "HMM.h"
#include "Latent.h"
#include <vector>

class Individual;
class Chromosome:public Genome
{
private:
   int _startLoci;
   int populations;
   int D;
   int L;
   std::string _Label;
   HMM SampleStates;
   Matrix_d StationaryDist;
   MatrixArray_d TransitionProbs, Likelihood;

  // UNIMPLEMENTED
  // to avoid use
  Chromosome(); // Private default constructor
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);

public:
   Chromosome(int size,int start, int);
   void ResetStuffForX();
  virtual ~Chromosome();
  virtual void
  SetLabel( int, std::string );
  virtual std::string
  GetLabel( int );
  virtual void
  accept(LocusVisitor&);
  virtual int
  GetLocus(int);
  virtual int
  GetSize();
  virtual Vector_i
  SampleForHaploidLocusAncestry(Individual*);
  virtual void
  UpdateParameters(Individual*,Matrix_d&,AdmixOptions*,std::vector<Vector_d>&,bool);
  virtual void
  UpdateParametersHaploid(Individual*,Matrix_d&,AdmixOptions*,std::vector<Vector_d>&,bool);
  virtual Matrix_i
  SampleForLocusAncestry(Individual*);
  virtual Matrix_d
  getExpectedAncestry( int );
  virtual double
  getLogLikelihood();

};

#endif /* !defined CHROMOSOME_H */
