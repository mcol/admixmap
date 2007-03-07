#include "UserData.h"

UserData::UserData()
{
}

UserData::~UserData()
{
  loci.clear();
}

void UserData::setLoci(vector<UserDataLocus> tmp_loci)
{
  loci = tmp_loci;
  vector<UserDataLocus>::iterator i;
  // Every locus should be able to access its parent UserData.
  for (i = loci.begin(); i != loci.end(); i++) {
    i->setUserData(this);
  }
  // Sort the loci
  sort(loci.begin(), loci.end());
}

void UserData::setIndivs(vector<UserDataIndividual> _indivs) {
  individuals = _indivs;
}

void UserData::confess(void) {
  
  vector<UserDataIndividual>::iterator i;
  for (i = individuals.begin(); i != individuals.end(); i++) {
    i->confess();
  }
  
  vector<UserDataLocus>::iterator j;
  for (j = loci.begin(); j != loci.end(); j++) {
    j->confess();
  }
}

vector<UserDataIndividual>& UserData::getIndividuals() {
  return individuals;
}

vector<UserDataLocus>& UserData::getLoci()
{
  return loci;
}
