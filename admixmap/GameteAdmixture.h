// *-*-C++-*-*
#ifndef GAMETE_ADMIXTURE_H
#define GAMETE_ADMIXTURE_H 1

#include "vector_d.h"
#include <cassert>

/**
 * The GameteAdmixture Class.
 * 
 * Represents the admixture proportions on a single gamete.
 * From this, it is possible to calculate the transition matrix
 * for all possible haploid states.
 *
 * To construct a GameteAdmixture object, either
 * specify the number of possible haploid states,
 * or supply a vector of proportions:
 * 
 *   int H = 6;
 *   GameteAdmixture ga(H);
 *   ga[0] = 0.1;
 *   ga[1] = 0.2;
 *   ga[2] = 0.05;
 * 
 * and so on. Or, using a vector of doubles,
 *
 *   GameteAdmixture ga(vec1);
 *
 */

class GameteAdmixture
{
private: // members
  unsigned int _H; // possible haploid states
  double *_M;      // Gamete admixture

  // UNIMPLEMENTED
  // - to prevent use
  GameteAdmixture();
  GameteAdmixture(const GameteAdmixture&);
  GameteAdmixture& operator=(const GameteAdmixture&);

public:

  // CONSTRUCTORS
  GameteAdmixture(unsigned int);
  GameteAdmixture(Vector_d);

  //DESTRUCTOR
  ~GameteAdmixture();

  unsigned int size();
  double& operator[](unsigned int);
};

#endif /* !defined GAMETE_ADMIXTURE_H */
