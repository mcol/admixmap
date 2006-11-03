/*
  Purpose: to sort the rows of the hapmap data files by map position. ~2500 are out of sequence.
*/

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include "ranker.h"
#include "gsl/gsl_permutation.h"

int main(int argc, char** argv)
{
  string method = "default";  // Can also be "min" or "max" or "average"
  std::ifstream infile(argv[1]);
  std::ofstream outfile(strcat(argv[1],"_sorted"));
    
  vector<string >  data;
  string line;
  vector<double> pos;
  string s;
  double d;

  getline(infile, s);//read header
  outfile << s;//and output
  infile >> s;

  while(!infile.eof()){
    for(unsigned i = 1; i < 4; ++i){
      line.append(s+" ");
      infile >> s;
    }
    pos.push_back(atol(s.c_str()));
    line.append(s+" ");
    getline(infile, s);//rest of line
    line.append(s);
    data.push_back(line);
    line.clear();
    infile >> s;
  }

  uint num = pos.size();
  vector<double> ranks;
  rank(pos, ranks, method);

  gsl_permutation * p = gsl_permutation_alloc (num);
  for (uint i = 0; i < num; ++i){
    p->data[i] = (int)ranks[i]-1;
  }

  //if(gsl_permutation_valid (p)){
    gsl_permutation * q = gsl_permutation_alloc (num);
    gsl_permutation_inverse (q, p);
    //gsl_permutation_fprintf (stdout, q, " %u");
    for (uint i = 0; i < num; ++i){
      ranks[i] = q->data[i];
    }

    //std::cout << data[(uint)ranks[0]] << std::endl;
    for(uint i = 0; i < num; ++i){
      outfile << data[(uint)ranks[i]] << std::endl;
    }

    //}
    //else std::cout << "Invalid permutation" << std::endl;


  return 0;

} // main
