#ifndef USERDATA_H_
#define USERDATA_H_

#include <iostream>
#include <string>
#include <vector>

#include "UserDataIndividual.h"
#include "UserDataLocus.h"

using std::string;
using std::vector;

class UserData
{
private:
  string id;
  vector<UserDataLocus> loci;
  vector<UserDataIndividual> individuals;
public:
  UserData();
  virtual ~UserData();
  void setLoci(vector<UserDataLocus>);
  void setLegend(Legend *);
  void printSnps();
  vector<UserDataIndividual>& getIndividuals();
  void setIndivs(vector<UserDataIndividual>);
  void confess(void);
  vector<UserDataLocus>& getLoci();
};

#endif /*USERDATA_H_*/
