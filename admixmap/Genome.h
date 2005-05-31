// *-*-C++-*-*
#ifndef GENETIC_ARRAY_H
#define GENETIC_ARRAY_H 1

#include "CompositeLocus.h"
//#include "Chromosome.h"
#include "vector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "LogWriter.h"
#include "AdmixOptions.h"
#include "InputData.h"

class Chromosome;//declared here to avoid circular includes
                 //not necessary if getChromosomes function moved out
class Genome 
// Object is an array of objects of class CompositeLocus
// Used to loop over composite loci on a single chromosome, or an entire genome  
{

public:

  Genome();
  Genome(int);
  virtual ~Genome();

  std::vector< int > GetChrmAndLocus( int );

  std::vector<  std::vector< int > >GetChrmAndLocus( );

   bool isX_data();

  CompositeLocus *&
  operator()(int) const;

  void SetLabels(const std::vector<std::string> &labels, Vector_d temp);

  void loadAlleleStatesAndDistances(vector<string> * ChrmLabels,AdmixOptions *options,InputData *data_, LogWriter *Log);

  Vector GetDistances();

  unsigned int GetNumberOfCompositeLoci();

  unsigned int GetNumberOfChromosomes();

  unsigned int GetTotalNumberOfLoci();

  int getNumberOfLoci(int);

  unsigned int *GetSizesOfChromosomes();

  float GetDistance(int);

  void SetDistance(int,float);

  Chromosome **GetChromosomes(int, std::vector<std::string> &);

  void SetSizes();
 
  int GetNumberOfStates();
  
  double GetLengthOfGenome();
  double GetLengthOfXchrm();

protected:// to make available to Chromosome
  Vector Distances;
  unsigned int NumberOfCompositeLoci;
  CompositeLocus **TheArray;

private: 

  double LengthOfGenome;
  double LengthOfXchrm;

  unsigned int NumberOfChromosomes;
  unsigned int TotalLoci;//number of simple loci;
  unsigned int *SizesOfChromosomes;
  bool X_data;
  std::vector< std::vector< int > > _chrmandlocus;

  void InitialiseCompositeLoci();

  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);
 
};

#endif /* !GENETIC_ARRAY_H */
