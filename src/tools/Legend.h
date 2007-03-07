#ifndef LEGEND_H_
#define LEGEND_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::stringstream;
using std::ifstream;
using std::istream;

/**
 * Locus structure. It needs to be tight because the whole genome's
 * legend needs to fit in the memory.
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
struct locus_pointer_t {
  locus_t *p;
};

/**
 * Struct to hold encoding of forward and reverse strands
 */
struct strand_coding_t {
  map<char, char> c;
};

/**
 * Struct to hold chromosome number and position
 */
struct locus_pos_t {
  signed char chromosome;
  unsigned long position;
};

/**
 * Comparator for loci positions, based on the chromosome
 * and base (position).
 */
struct locus_pos_t_cmp {
  bool operator() (const locus_pos_t s1, const locus_pos_t s2) const
  {
    if (s1.chromosome < s2.chromosome) {
      return -1;
    } else if (s1.chromosome > s2.chromosome) {
      return 1;
    } else {
      if (s1.position < s2.position) {
        return -1;
      } else if (s1.position > s2.position) {
        return 1;
      } else {
        return 0;
      }
    }
  }
};

/**
 * Legend, a class which represents the Hapmap project's legend,
 * supporting SNP lookups by rs identifier or chromosome and position.
 * 
 * Capable of loading the legend of the whole genom.
 *
 * Author: Maciej Blizinski
 */

class Legend
{	
public:
  Legend(string);
  Legend(ifstream&);
  Legend(stringstream&);
  void initialize(istream&);
  
  string getChromosomeBySnp(const string&);
  virtual ~Legend();
  short int getSnpNumberValueByBase(const string&, char);
  locus_t getLocusBySnp(const string &);
  locus_t *getLocusPointerBySnp(const string& snp_id);
  bool isInitialized();
protected:
  bool initialized;
  string file_name;
  map<const unsigned long, locus_pointer_t> snp_by_id;
  map<const locus_pos_t, locus_pointer_t, locus_pos_t_cmp> snp_by_pos;
  map<string, strand_coding_t> strandCoding;
  
  unsigned long getSnpNumberByString(const string&);
  
  void dbg_locus_t(locus_t);
  short int getSnpNumberValueByBaseAndEncoding(
    const locus_t *, char, const string&);
};

#endif /*LEGEND_H_*/
