#ifndef OPERATORS_H_
#define OPERATORS_H_

#include <mockpp/mockpp.h>

String& operator<<(String& s, const vector<double>&)
{
  s << "vector<double>";
  return s;
}

#endif /*OPERATORS_H_*/
