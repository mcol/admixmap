// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Chromosome.h 
 *   header file for Chromosome class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "Genome.h"
#include "HMM.h"
#include <vector>

class Individual;
class AdmixOptions;

class Chromosome:public Genome
{
public:
  Chromosome(int size,int start, int);
  void ResetStuffForX();
  ~Chromosome();
  void SetLabel(std::string );
  std::string GetLabel( int );
  int GetLocus(int);
  unsigned int GetSize();
  void isDiploid(bool b){Diploid = b;};
  bool isDiploid(){return Diploid;};

  void InitialiseLociCorr(const double rho);
  void SetLociCorr(const double rho);

  void UpdateHMMForwardProbs(Individual* ind, double *Admixture, AdmixOptions* options,  
			     std::vector< double > _rho, bool chibindicator, bool diploid);
  void UpdateHMMBackwardProbs(double* hapAdmixture);  //call only after a call to UpdateHMMForwardProbs
  //this is ok as whenever we need backward probs we also need forward probs but not vice versa


  void SampleLocusAncestry(int *OrderedStates, double *Admixture);
  void getAncestryProbs(int, double[][3]);
  double getLogLikelihood();
  void SampleJumpIndicators(int *LocusAncestry, const unsigned int gametes, 
			    int *SumLocusAncestry, int *SumLocusAncestry_X, bool isX, 
			    unsigned int SumN[], unsigned int SumN_X[], bool RhoIndicator);
private:
  int _startLocus;
  int populations;
  int D; 
  std::string _Label;
  HMM SampleStates;
  double *Lambda;
  bool Diploid;
  
  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f[2]; 
  int *CodedStates;//used to sample hidden states from HMM
  
  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome();
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);
};

#endif /* !defined CHROMOSOME_H */
