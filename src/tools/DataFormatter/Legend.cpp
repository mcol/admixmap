#include "Legend.h"
#include <iostream>
#include <fstream>
#include <sstream>

#define RS_LIMIT_READ 0

using namespace std;

Legend::Legend()
{
}

/**
 * Constructor, reads a file with data, tab-separated columns
 * rsXXXX <position>  <allele1> <allele2> <chromosome>
 */
Legend::Legend(string _file_name)
{
  file_name = _file_name;
  string line, buf;
  stringstream ss;
  string rs, chromosome;
  unsigned long position, rs_no;
  char a1, a2;
  int counter = 0;
  locus_t *l;
  locus_p_t lp;

  ifstream data_file(file_name.c_str());

  if (data_file.is_open()) {
    while (!data_file.eof()) {
      counter++;

      // Read data and insert a value into rs map.
      getline(data_file, line);
      stringstream ss(line);
      // Read the identifier
      ss >> rs;
      ss >> position;
      ss >> a1;
      ss >> a2;
      ss >> chromosome;
      l = (locus_t*)malloc(sizeof(locus_t));
      if (l == NULL) {
      	// memory not allocated
      	exit(EXIT_FAILURE);
      }
      rs_no =atoi(rs.substr(2, rs.size() - 1).c_str());
      l->rs = rs_no;
      l->position = position;
      l->a1 = a1;
      l->a2 = a2;
      l->chromosome = (char)atoi(chromosome.substr(3, chromosome.size() - 2).c_str());
      lp.p = l;
      rs_map[rs_no] = lp;
      rs_by_position[position] = lp;
      // dbg_locus_t(*(rs_map[rs_no].p));
      
      if (RS_LIMIT_READ && (counter % (int)1e5 == 0)) {
      	cerr << "Rs no. " << counter << endl;
      }
      
      if (RS_LIMIT_READ && (counter > (int)1e5)) {
      	break;
      }

    }
  } else {
    cerr << "Can't open the data file." << endl;
    exit(EXIT_FAILURE);
  }
}

Legend::~Legend()
{
  map<const unsigned long, locus_p_t>::iterator i;
  for(i = rs_map.begin(); i != rs_map.end(); i++) {
    free(i->second.p);
  }
  rs_map.clear();
  rs_by_position.clear();
}

string Legend::get_chromosome_by_snp(string _snp)
{
  unsigned long snp_no;
  stringstream ss;
  string ret;
  snp_no = atoi(_snp.substr(2, _snp.size() - 1).c_str());
  // cerr << "'" << _snp.substr(2, _snp.size() - 1).c_str() << "'" << endl;
  // cerr << "Looking up snp number " << snp_no << endl;
  map<const unsigned long, locus_p_t>::iterator i;
  i = rs_map.find(snp_no);
  if (i == rs_map.end()) {
    cerr << "SNP '" << _snp << "' doesn't exist!" << endl;
    return("");
  } else {
    ss << "chr" << (unsigned int)i->second.p->chromosome;
    ss >> ret;
    return ret;
  }
}

/**
 * Print locus contents
 */
void Legend::dbg_locus_t(locus_t l) {
  cout << l.rs << "\t"
    << l.position << "\t"
    << l.a1 << "\t"
    << l.a2 << "\t"
    << "chr" << (int)l.chromosome << endl;
}

