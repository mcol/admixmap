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

void print_rs(
    const map<const string, string>&  rs,
          string                      key)
{
    if (rs.find(key) == rs.end()) {
      cerr << "Key '" << key << "' doesn't exist!" << endl;
    }  else {
      // cout << "'" << key << "', '" << rs.find(key)->second << "'" << endl;
      cout << rs.find(key)->second << endl;
    }
}

int main(void) {

  string line, buf;
  stringstream ss;
  vector<string> tokens;
  map<const string, string> rs;
  int counter = 0;

  ifstream data_file("rs_data.csv");

  if (data_file.is_open()) {
    while (!data_file.eof()) {
      counter++;
      // Empty the tokens vector
      tokens.clear();

      // Read data and insert a value into rs map.
      getline(data_file, line);
      stringstream ss(line);
      while (ss >> buf) {
        tokens.push_back(buf);
      }
      if (tokens.size() == 2) {
        rs[tokens[0]] = tokens[1];
      } else {
        cerr << "Wrong line: '" << line << "'." << endl;
      }

      // if (counter % (int)1e5 == 0) {
      //  cerr << "Rs no. " << counter << endl;
      // }
    }
  }

  // Read the input and write the output
  while (cin >> buf) {
    print_rs(rs, buf);
  }

  rs.clear();

  return(EXIT_SUCCESS);
}
