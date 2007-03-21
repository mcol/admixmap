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
  sumProbs.assign(Nindivs, vector<vector<double> >(Nloci, vector<double>(3, 0)));
}

/// Accumulates genotype probabilities for individual i, locus j
void GenotypeProbOutputter::Update(
    const unsigned i, //< Individual number
    const unsigned j, //< Locus number
    const vector<vector<double> >& unorderedGenotypeProbs)
{
#ifndef PARALLEL
#define GPO_SPIT cout << "GenotypeProbOutputter::Update(): (" \
    << i << ", " << j \
    << ") (" \
    << unorderedGenotypeProbs[0][0] << ", " \
    << unorderedGenotypeProbs[1][0] << ", " \
    << unorderedGenotypeProbs[2][0] << "), (" \
    << sumProbs[i][j][0] << ", " \
    << sumProbs[i][j][1] << ", " \
    << sumProbs[i][j][2] << ")" << endl

//  if (isnan(unorderedGenotypeProbs[0][0])
//   || isnan(unorderedGenotypeProbs[1][0])
//   || isnan(unorderedGenotypeProbs[2][0])) {
//    GPO_SPIT;
//    throw("unorderedGenotypeProbs are nan");
//  }

  //TODO: shortcut - HapPairs will always be the same

//   if(fabs(Probs[0]+Probs[1]+Probs[2]+Probs[3] - 1.0) >1e-6 ){
//     cout << Probs[0] << " " << Probs[1] << " " << Probs[2] << " " << Probs[3] << endl;
//     cout << "j= " << j << " i= " << i << " " << (j*NumMaskedIndivs +i )*3 << " ";
//   }

//  //cout << "J = " << NumMaskedLoci << " I = " << NumMaskedIndivs << " " << SumProbs.size() << " "<< SumProbs[(j*NumMaskedIndivs +i )*3 ] << endl;

  for (int g = 0; g < 3; ++g) {
    sumProbs[i][j][g] += unorderedGenotypeProbs[g][0];
  }
  
//  if (i == 1) {
//    GPO_SPIT;
//  }
//  if (isnan(sumProbs[i][j][0])
//   || isnan(sumProbs[i][j][1])
//   || isnan(sumProbs[i][j][2])) {
//    GPO_SPIT;
//    throw("sumProbs are nan");
//  }

  ++NumIterations;

#endif
}

//! Output Probs as R object to be read with dget() function
void GenotypeProbOutputter::Output(const char* filename){
  outfile.open(filename);
  outfile << "# $Id$" << endl;
  outfile << "structure(.Data=c(" << endl;
  outfile << setfill(' ');
  //  outfile.setf(ios::fixed); 
  outfile.precision(6);
  //average over iterations
  NumIterations /= NumMaskedIndivs;
  NumIterations /= NumMaskedLoci;
  double numIterationsReciprocal = 1.0 / NumIterations;
  
  vector<vector<vector<double> > >::iterator indiv_i;
  vector<vector<double> >::iterator locus_i;
  vector<double>::iterator gt_i;
  
  unsigned counter = 0;
  unsigned icounter1 = 1; //< Individual counter (1-based)
  unsigned lcounter1 = 1; //< Locus counter (1-based)
  unsigned last = NumMaskedLoci * NumMaskedIndivs * 3 - 1;
  for(indiv_i = sumProbs.begin(); indiv_i != sumProbs.end(); ++indiv_i) {
    outfile << "# Individual " << icounter1 << endl;
    lcounter1 = 1;
    for (locus_i = indiv_i->begin(); locus_i != indiv_i->end(); ++locus_i) {
      for (gt_i = locus_i->begin(); gt_i != locus_i->end(); ++gt_i) {
        double tmp_gt = *gt_i;
        *gt_i *= numIterationsReciprocal;
        // Round the probabilites to 6 decimal places
        // Makes them sum up to exactly 1.0 after outputting
        // (assuming they do sum up to 1.0 before rounding)
        double tmp_gt2 = *gt_i;
        *gt_i = floor((*gt_i) * 1e6 + 0.5) * 1e-6;
        // Output the number
        outfile << setw(10);
        outfile << *gt_i;
        // Write a comma if it wasn't the last element
        if (!(counter == last)) { outfile << ", "; }
        ++counter;
        if (isnan(*gt_i)) {
          cerr << "*gt_i is nan. "
              << tmp_gt << " * "
              << numIterationsReciprocal << " -> " 
              << tmp_gt2 << " -> "
              << "nan" << endl;
          throw("*gt_i is nan");
        }
      }
      outfile << "# masked locus no. " << lcounter1;
      // Keep one line per genotype triplet
      outfile << endl;
      ++lcounter1;
    }
    ++icounter1;
  }

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
