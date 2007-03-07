#include "UserDataIndividual.h"

UserDataIndividual::UserDataIndividual()
{
}

UserDataIndividual::UserDataIndividual(string _id, int _outcome)
{
  id = _id;
  outcome = _outcome;
}

UserDataIndividual::~UserDataIndividual()
{
}

void UserDataIndividual::setBaseGenotypes(vector<base_genotype_pair_t> gpairs) {
  base_gpairs = gpairs;
}

void UserDataIndividual::confess(void) {
  cout << id << ", " << outcome << ", ";
  vector<base_genotype_pair_t>::iterator i;
  for (i = base_gpairs.begin(); i != base_gpairs.end(); i++) {
    if (i->present) {
      cout << i->g[0] << i->g[1] << " ";
    } else {
      cout << "?? ";
    }
  }
  cout << endl;
}

base_genotype_pair_t UserDataIndividual::getBaseGenotype(unsigned int i) {
  if (i >= base_gpairs.size()) {
    // FIXME: perhaps raise an exception
  }
  return base_gpairs[i];
}

string& UserDataIndividual::getId()
{
  return id;
}

vector<base_genotype_pair_t>& UserDataIndividual::getBaseGenotypes()
{
  return base_gpairs;
}

short int UserDataIndividual::getOutcome()
{
  return outcome;
}
