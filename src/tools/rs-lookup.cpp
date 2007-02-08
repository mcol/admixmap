/**
 * 
 * Fast look-up, determines chromosome given a locus identifier.
 *
 * Maciej Blizinski <maciej.blizinski@ucd.ie>
 *
 */

#include <cstdlib>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

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

/**
 * Print the chromosome given the rs identifier.
 */
void print_chr(
    const map<const string, locus_p_t>& rs,
          string                      key)
{
    if (rs.find(key) == rs.end()) {
      cerr << "Key '" << key << "' doesn't exist!" << endl;
    }  else {
      cout << rs.find(key)->second.p->chromosome << endl;
      // cerr << key << ": " << rs.find(key)->second.chromosome
      //   << " (" << &rs.find(key)->second.chromosome << ")" << endl;
    }
}

void read_rs_data(
    string                        file_name,
    map<const string, locus_p_t>& rs_map,
    map<const long, locus_p_t>&   rs_by_position)
{
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
      // Locus locus = Locus(rs, position, a1, a2, chromosome);
      // rs_map[rs] = Locus(rs, position, a1, a2, chromosome);
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

int main(void) {

  string buf;
  map<const string, locus_p_t> rs_map;
  map<const long, locus_p_t> rs_by_position;
  vector<string> in_rs;

  read_rs_data("rs_data.csv", rs_map, rs_by_position);

  // Read the input and write the output
  while (cin >> buf) {
	print_chr(rs_map, buf);
        in_rs.push_back(buf);
  }

  map<const string, locus_p_t>::iterator i;
  for(i = rs_map.begin(); i != rs_map.end(); i++) {
    free(i->second.p);
  }
  rs_map.clear();
  rs_by_position.clear();

  return(EXIT_SUCCESS);
}
