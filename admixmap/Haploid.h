// *-*-C++-*-*
#ifndef HAPLOID_H
#define HAPLOID_H 1

#include <iostream>
#include <cmath>
#include <vector>
#include "vector_i.h"

/**
 * The Haploid Class.
 *
 * A class that holds a pair of haploid states.
 * This is a little different than PMcK's desciption,
 * so I'll explain a little.
 * 
 * If we have two haploid states (h0 and h1), we
 * can draw a table of possible diploid states thus:
 *
 *                   h0
 *
 *        |  0  1  2  3  4  5  6
 *     --- ---------------------
 *      0 |  0
 *        |
 *      1 |  1  2
 *        |
 *      2 |  3  4  5
 *        |
 *  h1  3 |  6  7  8  9
 *        |
 *      4 | 10 11 12 13 14
 *        |
 *      5 | 15 16 17 18 19 20
 *        |
 *      6 | 21 22 23 24 25 26 27
 *
 *
 * Construct using either both haploid states:
 *
 *   Haploid h1(3,4);
 *   Haploid h2(4,3); // same as h1
 *
 * or with just the corresponding diploid state:
 * 
 *   Haploid h3(13);
 *
 * 
 * A Haploid object can be sent to an output stream:
 *
 *   cout << h1; // ALL THESE
 *   cout << h2; // RESULT IN
 *   cout << h3; // OUTPUT OF '3,4'
 *
 *
 * The diploid state can be discovered using 'toDiploid()':
 *
 *   cout << h1.toDiploid(); // PRINTS '13'
 *
 * Each haploid state can be found by subscripting:
 *
 *   cout << h3[0]; // prints '3'
 *   cout << h3[1]; // prints '4'
 *
 * Equality and inequality are overridden in a sensible way.
 *
 *   (h1 == h2) // true
 *   (h2 != h3) // false
 *
 */

class HaploidStates;

class Haploid
{
public:
  Haploid();
  Haploid(unsigned int);
  Haploid(unsigned int,unsigned int);
  virtual unsigned int toDiploid();
  virtual Vector_i toVector();
  virtual unsigned int operator[](unsigned const int x) const;

  friend bool operator==(const Haploid left, const Haploid right);
  friend bool operator!=(const Haploid left, const Haploid right);
  friend std::ostream& operator<<(std::ostream& os, const Haploid& hs);

private:
  unsigned int _hap[2];
  static HaploidStates _haploidStates;
};

// Contains some precalculated values for use in Haploid constructor.
// This improves the speed of admixmap

class HaploidStates{
public:
  HaploidStates(int size);  // Bug in g++ means have to have argument as g++ won't initialize a static object with the default constructor.
  unsigned int getState(unsigned int d, unsigned int i);
private:
  std::vector<unsigned int> _haploidStates0;
  std::vector<unsigned int> _haploidStates1;
};

#endif /* !defined HAPLOID_H */
