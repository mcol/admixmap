// *-*-C++-*-*
#ifndef GENETIC_ARRAY_H
#define GENETIC_ARRAY_H 1

//#include "AbstractCompLocus.h"
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
private: // members
  Vector Distances;
  double LengthOfGenome;
  double LengthOfXchrm;
  int NumberOfCompositeLoci;
  int NumberOfChromosomes;
  bool X_data;
  CompositeLocus **TheArray;
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

  // composite-level methods
  //AbstractCompLocus*&
  CompositeLocus *&
  operator()(int) const;

  void SetLabels(const std::vector<std::string> &labels, Vector_d temp);

  Vector GetDistances();

  void SetNumberOfCompositeLoci(int);

  int GetNumberOfCompositeLoci();

  int GetNumberOfChromosomes();

  float GetDistance(int);

  void SetDistance(int,float);

  Chromosome **GetChromosomes(int, std::vector<std::string> );

  int size();

  // individual-level methods

  int GetNumberOfLoci();
  
  int GetNumberOfStates();
  
  virtual int GetSize();

  virtual void SetLabel( int, std::string );
  
  double GetLengthOfGenome();
  double GetLengthOfXchrm();

};

#endif /* !GENETIC_ARRAY_H */
