#ifndef USERDATAINDIVIDUAL_H_
#define USERDATAINDIVIDUAL_H_

#include <iostream>
#include <vector>

using std::string;
using std::vector;
using std::cout;
using std::endl;

/**
 * Structure to hold a single genotype.
 */
struct base_genotype_pair_t {
        unsigned char present;
        char g[2];
};

class UserDataIndividual
{
private:
  string id;
  int outcome;
  /** List of pairs of genotypes. */
  vector<base_genotype_pair_t> base_gpairs;
public:
  UserDataIndividual();
  UserDataIndividual(string, int);
  base_genotype_pair_t getBaseGenotype(unsigned int);
  virtual ~UserDataIndividual();
  void setBaseGenotypes(vector<base_genotype_pair_t>);
  void confess();
  string& getId();
  vector<base_genotype_pair_t>& getBaseGenotypes();
  short int getOutcome();
};

#endif /*USERDATAINDIVIDUAL_H_*/
