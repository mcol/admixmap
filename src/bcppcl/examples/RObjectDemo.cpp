#include "../bcppcl/RObjectWriter.h"
#include <vector>
#include <string>

using std::vector;
using std::string;
int main(){

  RObjectWriter R("testR.txt");

  R << "label1" << 1 << 2 << 3 << 4 << Rcomment("line1") << newline
    << Rcomment("blank line") << newline
    << "label2" << 5 << 6 << 7 << 8 << Rcomment("line2") << newline;

  vector<int> dims;
  dims.push_back(5);
  dims.push_back(2);
  vector<string> colnames;
  colnames.push_back("Labels");
  colnames.push_back("col1");
  colnames.push_back("col2");
  colnames.push_back("col3");
  colnames.push_back("col4");
  vector<vector<string> >  dimnames;
  dimnames.push_back(colnames);

  R.close(dims, dimnames);


  return 0;
}
