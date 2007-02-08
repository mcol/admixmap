#ifndef LEGEND_H_
#define LEGEND_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

/**
 * Locus structure
 */
struct locus_t {
  char *rs;
  long position;
  char a1;
  char a2;
  char *chromosome;
};

/**
 * Structure containing a pointer to a locus
 */
struct locus_p_t {
  locus_t *p;
};


class Legend
{	
public:
  Legend();
  Legend(string);
  char *get_chromosome_by_snp(string);
  virtual ~Legend();
protected:
  map<const string, locus_p_t> rs_map;
  map<const long, locus_p_t> rs_by_position;
  void dbg_locus_t(locus_t);
  string file_name;
};

#endif /*LEGEND_H_*/
