// *-*-C++-*-*
#ifndef HAPLOID_PROBABILITIES_H
#define HAPLOID_PROBABILITIES_H 1

#include <cassert>
#include "Haploid.h"
#include "vector_d.h"
#include "matrix_d.h"

/**
 * The HaploidProbabilities Class.
 *
 * A class that holds the probabilities of each haploid state
 * on paternal and maternal gametes. With this information,
 * we can calulate the probabilities of each diploid state,
 * or the probability that the maternal and paternal states
 * are in a particular phase, given the diploid state.
 *
 * Construct in one of two ways. Either give the number
 * of possible haploid states on each gamete:
 * 
 *   unsigned int H = 6;
 *   HaploidProbabilities hp(H);
 *
 * or, supply vectors for probabilities for each haploid
 * state on the two gametes:
 *
 *   HaploidProbabilities hp(vec1,vec2);
 *
 * To find the probabilities for a certain diploid state,
 * use the getDiploidProbability() method:
 *
 *   unsigned int d=3;
 *   double dP=hp.getDiploidProbability(d);
 *
 * Diploid probability can also be found from a haploid pair:
 *
 *   Haploid h(2,3);
 *   double dP=hp.getDiploidProbability(h);
 *
 * To find the probabilities for the maternal and paternal
 * states being in a particular phase, use the getPhaseProbability()
 * method:
 *
 *   double pP=hp.GetPhaseProbability(d);
 * or,
 *   double pP=hp.GetPhaseProbability(h);
 *
 * (This gives the probability that maternal gamete is h[0] and
 * paternal gamete is h[1]).
 *
 */
class HaploidProbabilities
{
private: // members
  unsigned int _hapStates;
  double *_maternal;
  double *_paternal;


  // UNIMPLEMENTED
  // - to prevent use
  HaploidProbabilities(); 
  HaploidProbabilities(const HaploidProbabilities&);
  HaploidProbabilities& operator=(const HaploidProbabilities&);
  
public:

  // CONSTRUCTORS
  HaploidProbabilities(unsigned int);
  HaploidProbabilities(Vector_d,Vector_d);

  // DESTRUCTOR
  ~HaploidProbabilities();

  double getDiploidProbability(Haploid&);
  double getDiploidProbability(unsigned int);
  Matrix_d getDiploidProbabilities();
  double getPhaseProbability(Haploid&);
  double getPhaseProbability(unsigned int);
  Vector_d getPhaseProbabilities();
  void setMaternal(unsigned int,double);
  void setMaternal(Vector_d);
  void setPaternal(unsigned int,double);
  void setPaternal(Vector_d);
  
  unsigned int size();
};


#endif /* HAPLOID_PROBABILITIES_H !defined */
