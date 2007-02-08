/**
 * 
 * Fast look-up, determines chromosome given a locus identifier.
 *
 * Maciej Blizinski <maciej.blizinski@ucd.ie>
 *
 */

#include <cstdlib>
#include <map>
#include "Legend.h"

using namespace std;

int main(void) {

  string buf;
  vector<string> in_rs;

  Legend legend = Legend("rs_data.csv");

  // Read the input and write the output
  while (cin >> buf) {
        cout << legend.get_chromosome_by_snp(buf) << endl;
        in_rs.push_back(buf);
  }

  return(EXIT_SUCCESS);
}
