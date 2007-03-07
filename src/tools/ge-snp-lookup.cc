/**
 * Fast look-up, determines chromosome given a locus identifier.
 *
 * Maciej Blizinski <maciej.blizinski@ucd.ie>
 */

#include "ge-snp-lookup.h"

int main(int argc, char **argv) {

  string buf, rs_str;
  vector<string> in_rs;

  if (argc < 2) {
    cout << endl;
    cout << "Usage: " << argv[0] << " <file>" << endl;
    cout << endl;
    exit(EXIT_FAILURE);
  }

  string file_name = string(argv[1]);

  Legend legend = Legend(string(file_name));

  // Read the input and write the output
  while (cin >> buf) {
    if (buf.size() == 0) {
      continue;
    }
    cout << legend.getChromosomeBySnp(buf) << endl;
    in_rs.push_back(buf);
  }

  return(EXIT_SUCCESS);
}
