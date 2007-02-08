#include "Legend.h"

Legend::Legend()
{
}

Legend::Legend(string _file_name)
{
  file_name = _file_name;
  string line, buf;
  stringstream ss;
  string rs, chromosome;
  long position;
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
      	// not allocated the memory
      	exit(1);
      }
      l->rs = (char *)malloc(sizeof(char) * rs.size() + 1);
      if (l->rs == NULL) {
        exit(1);
      }
      strcpy(l->rs, rs.c_str());
      l->position = position;
      l->a1 = a1;
      l->a2 = a2;
      l->chromosome = (char *)malloc(sizeof(char) * chromosome.size() + 1);
      if (l->chromosome == NULL) {
        exit(1);
      }
      strcpy(l->chromosome, chromosome.c_str());
      // dbg_locus_t(*l);
      lp.p = l;
      rs_map[rs] = lp;
      rs_by_position[position] = lp;
      // dbg_locus_t(rs_map[rs]);
      
      if (counter % (int)1e5 == 0) {
      	cerr << "Rs no. " << counter << endl;
      }
      
      if (1 && (counter > (int)1e5)) {
      	break;
      }
    }
  }
}

Legend::~Legend()
{
  map<const string, locus_p_t>::iterator i;
  for(i = rs_map.begin(); i != rs_map.end(); i++) {
    free(i->second.p);
  }
  rs_map.clear();
  rs_by_position.clear();
}

char *Legend::get_chromosome_by_snp(string _snp)
{
  if (rs_map.find(_snp) == rs_map.end()) {
    cerr << "SNP '" << _snp << "' doesn't exist!" << endl;
    return("");
  } else {
    return rs_map.find(_snp)->second.p->chromosome;
  }
}

/**
 * Print locus contents
 */
void dbg_locus_t(locus_t l) {
  cout << l.rs << "\t"
    << l.position << "\t"
    << l.a1 << "\t"
    << l.a2 << "\t"
    << l.chromosome << endl;
}

