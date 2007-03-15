#include "GenotypeProbOutputter.h"
#include <string>
#include <cmath>

using namespace::std;
void GenotypeProbOutputter::Initialise(unsigned Nindivs, unsigned Nloci){
  NumMaskedIndivs = Nindivs;
  NumMaskedLoci = Nloci;
  NumIterations = 0;
  Probs.assign(4, 0);
  //NB: assuming only SNPs.
  //TODO:?? extend to arbitrarily-sized loci (requires number of haplotype states for each locus) 
  SumProbs.assign(Nloci*Nindivs*3, 0.0);
}

void GenotypeProbOutputter::Update(unsigned i, unsigned j, const CompositeLocus* Locus, const std::vector<hapPair > &HapPairs, const int ancestry[2]){
#ifndef PARALLEL
  //TODO: shortcut - HapPairs will always be the same
  Locus->getConditionalHapPairProbs(Probs, HapPairs, ancestry);

//   if(fabs(Probs[0]+Probs[1]+Probs[2]+Probs[3] - 1.0) >1e-6 ){
//     cout << Probs[0] << " " << Probs[1] << " " << Probs[2] << " " << Probs[3] << endl;
//     cout << "j= " << j << " i= " << i << " " << (j*NumMaskedIndivs +i )*3 << " ";
//   }


  //cout << "J = " << NumMaskedLoci << " I = " << NumMaskedIndivs << " " << SumProbs.size() << " "<< SumProbs[(j*NumMaskedIndivs +i )*3 ] << endl;
  SumProbs[(j*NumMaskedIndivs +i )*3 ] += Probs[0];//genotype "1,1"
  SumProbs[(j*NumMaskedIndivs +i )*3 +1] += Probs[1] + Probs[2];//genotype "1,2"
  SumProbs[(j*NumMaskedIndivs +i )*3 +2] += Probs[3];//genotype "2,2"
  ++NumIterations;

#endif
}

//! Output Probs as R object to be read with dget() function
void GenotypeProbOutputter::Output(const char* filename){
  outfile.open(filename);
  outfile << "structure(.Data=c(" << endl;
  //  outfile.setf(ios::fixed); 
  outfile.precision(6);
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

  //output dimensions
  outfile << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    outfile << dim[i];
    if(i != dim.size() - 1){
      outfile << ",";
    }
  }
  outfile << ")," << endl;

  //output dimnames
  outfile << ".Dimnames=list(";

  outfile << "c(";
  for(unsigned int i=0;i<labels.size();i++){
    outfile << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      outfile << ", ";
    }
  }
  outfile << "), ";

  // Individual labels
  outfile << "1:" << dim[1];
  // outfile << "INDIVIDUAL_LABELS";
  outfile << ", ";
  // Locus labels
  outfile << "1:" << dim[2] << "";
  // outfile << "LOCI_LABELS";

  // Closing brackets, one to end list, one to end object.
  outfile << "))" << endl;
  outfile.close();
}
