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
struct locus_t
{
  unsigned long rs;
  unsigned long position;
  char          a1;
  char          a2;
  char          chromosome;
};

/**
 * Structure containing a pointer to a locus
 */
struct locus_p_t {
  locus_t *p;
};

/**
 * Class representing the legend. Capable of loading the legend of the
 * whole genotype.
 */

class Legend
{	
public:
  Legend();
  Legend(string);
  string get_chromosome_by_snp(string);
  virtual ~Legend();
protected:
  map<const unsigned long, locus_p_t> rs_map;
  map<const unsigned long, locus_p_t> rs_by_position;
  void dbg_locus_t(locus_t);
  string file_name;
};

#endif /*LEGEND_H_*/
