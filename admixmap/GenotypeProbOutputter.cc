#include "GenotypeProbOutputter.h"
#include <string>

using namespace::std;
void GenotypeProbOutputter::Initialise(unsigned Nindivs, unsigned Nloci){
  NumMaskedIndivs = Nindivs;
  NumMaskedLoci = Nloci;
  NumIterations = 0;
  SumProbs.assign(Nindivs*Nloci*3, 0.0);
}

void GenotypeProbOutputter::Update(unsigned i, unsigned j, const CompositeLocus* Locus, const std::vector<hapPair > &HapPairs, const int ancestry[2]){
#ifdef PARALLEL
  //TODO: shortcut - HapPairs will always be the same
  Locus->getConditionalHapPairProbs(Probs, HapPairs, ancestry);

  //  cout << Probs[0] << " " << Probs[1] << " " << Probs[2] << " " << Probs[3] << endl;
  //cout << "j= " << j << " i= " << i << " " << (j*NumMaskedIndivs +i )*3 << " ";
  //cout << "J = " << NumMaskedLoci << " I = " << NumMaskedIndivs << " " << SumProbs.size() << " "<< SumProbs[(j*NumMaskedIndivs +i )*3 ] << endl;
  SumProbs[(j*NumMaskedIndivs +i )*3 ] += Probs[0];//genotype "1,1"
  SumProbs[(j*NumMaskedIndivs +i )*3 +1] += Probs[1] + Probs[2];//genotype "1,2"
  SumProbs[(j*NumMaskedIndivs +i )*3 +2] += Probs[3];//genotype "2,2"
  ++NumIterations;
#endif
}

///Output Probs as R object
void GenotypeProbOutputter::Output(const char* filename){
  outfile.open(filename);
  outfile << "structure(.Data=c(" << endl;
  //average over iterations
  NumIterations /= NumMaskedIndivs * NumMaskedLoci;
  for(vector<double>::iterator i = SumProbs.begin(); i != SumProbs.end(); ++i)*i /= (double)NumIterations;
  //output to file, followed by comma
  copy(SumProbs.begin(), SumProbs.end()-1, ostream_iterator<double>(outfile, ", "));
  //output last element
  outfile << *(SumProbs.end()-1) << endl;

  //output dimensions
  std::vector<int> dim(3,0);
  dim[0] = 3;
  dim[1] = NumMaskedIndivs;
  dim[2] = NumMaskedLoci;
  
  std::vector<std::string> labels(dim[0],"");
  labels[0] = "Genotype1";
  labels[1] = "Genotype2";
  labels[2] = "Genotype3";
  outfile << ")," << endl;
  outfile << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    outfile << dim[i];
    if(i != dim.size() - 1){
      outfile << ",";
    }
  }
  outfile << ")," << endl;
  outfile << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    outfile << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      outfile << ",";
    }
  }
  outfile << "), numeric(" << dim[1] << "), numeric(" << dim[2] << ")))" << endl;
  outfile.close();
}
