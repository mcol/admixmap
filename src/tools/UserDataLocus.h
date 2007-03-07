#ifndef USERDATALOCUS_H_
#define USERDATALOCUS_H_

#include <iostream>
#include <string>
#include <set>
#include "Legend.h"
#include "UserDataIndividual.h"
// #include "UserData.h"

class UserData;
using std::string;
using std::set;
using std::ostream;
using std::cerr;

enum LocusCoding { FORWARD_STRAND, REVERSE_STRAND };

/**
 * Class representing a locus in user data. It doesn't necessarily hold
 * the data itself, but it provides interface to deal with given locus.
 */

class UserDataLocus
{
private:
  string id;
  unsigned int idx;
  LocusCoding coding;
  bool flipped;
  bool monomorphic;
  unsigned char number_alleles;
  bool faulty;
  bool analyzed;
  set<char> alleles;
  
  /** Pointer to the data this locus is a member of. */
  UserData *vpUserData;
  const locus_t *locusData;
public:
  UserDataLocus(unsigned int, locus_t *);
  void setUserData(UserData *);
  void setLegend(Legend *);
  void analyze();
  void confess();
  virtual ~UserDataLocus();
  string& getId();
  unsigned int getPosition();
  int getChromosome();
  char getAllele1();
  char getAllele2();
  
  bool isPresent(vector<UserDataIndividual>::iterator&);
  char getBaseAllele(
      vector<UserDataIndividual>::iterator&,
      unsigned char);
  short int getNumberAllele(
      vector<UserDataIndividual>::iterator&,
      unsigned char, Legend&);
  int compare(const UserDataLocus&) const;
  bool operator < (const UserDataLocus&) const;
  bool operator > (const UserDataLocus&) const;
  bool operator == (const UserDataLocus&) const;
  bool operator != (const UserDataLocus&) const;
  string toString();
};

#include "UserData.h"

#endif /*USERDATALOCUS_H_*/
