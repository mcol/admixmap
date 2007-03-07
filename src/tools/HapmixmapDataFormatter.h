#ifndef HAPMIXMAPDATAFORMATTER_H_
#define HAPMIXMAPDATAFORMATTER_H_

#include <cstdlib>
#include <string>
#include <ostream>
#include <iostream>
#include <vector>
#include "DataFormatter.h"
#include "UserDataLocus.h"

using std::string;
using std::ostream;
using std::endl;
using std::vector;

class HapmixmapDataFormatter : public DataFormatter
{
public:
	HapmixmapDataFormatter();
  HapmixmapDataFormatter(UserData *, Legend *);
	virtual ~HapmixmapDataFormatter();
  int getNoFiles();
  void renderFile(ostream&, int);
  void render_genotypes(ostream&, Legend&);
  void render_loci(ostream&);
  void render_outcome(ostream&);
};

#endif /*HAPMIXMAPDATAFORMATTER_H_*/
