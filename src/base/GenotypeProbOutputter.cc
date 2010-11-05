//=============================================================================
//
// Copyright (C) 2006  David O'Donnell
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file GenotypeProbOutputter.cc
/// Implementation of the GenotypeProbOutputter class.
//=============================================================================

#include "GenotypeProbOutputter.h"
#include <string>
#include <cmath>
#include "bclib/Exceptions.h"

//define to turn on checking for nans
//#define GPO_DEBUG
 
//define to output indiv and locus labels
//#define GPO_LABELS

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
void GenotypeProbOutputter::Update( unsigned indiv, unsigned locus, 
				    const vector<vector<double> >& unorderedGenotypeProbs){
#define GPO_SPIT cout << "GenotypeProbOutputter::Update(): (" \
    << indiv << ", " << locus \
    << ") (" \
    << unorderedGenotypeProbs[0][0] << ", " \
    << unorderedGenotypeProbs[1][0] << ", " \
    << unorderedGenotypeProbs[2][0] << "), (" \
    << sumProbs[indiv][locus][0] << ", " \
    << sumProbs[indiv][locus][1] << ", " \
    << sumProbs[indiv][locus][2] << ")" << endl

//  if (isnan(unorderedGenotypeProbs[0][0])
//   || isnan(unorderedGenotypeProbs[1][0])
//   || isnan(unorderedGenotypeProbs[2][0])) {
//    GPO_SPIT;
//    throw("unorderedGenotypeProbs are nan");
//  }

  //TODO: shortcut - HapPairs will always be the same

  for (int g = 0; g < 3; ++g) {
    sumProbs[indiv][locus][g] += unorderedGenotypeProbs[g][0];
  }
  
  ++NumIterations;

}

/// Output Probs as R object to be read with dget() function
void GenotypeProbOutputter::Output(const char* filename, const vector<string>& LocusLabels){
  outfile.open(filename);
  outfile.setDecimalPrecision(6);

  NumIterations /= NumMaskedIndivs;
  NumIterations /= NumMaskedLoci;
  double numIterationsReciprocal = 1.0 / NumIterations;
  
  vector<vector<vector<double> > >::iterator indiv_i;
  vector<vector<double> >::iterator locus_i;
  vector<double>::iterator gt_i;
  
  unsigned icounter1 = 1; ///< Individual counter (1-based)
  unsigned lcounter1 = 1; ///< Locus counter (1-based)
  double SumProb = 0.0, UnRoundedMeanProb = 0.0;

  for(indiv_i = sumProbs.begin(); indiv_i != sumProbs.end(); ++indiv_i) {
#ifdef GPO_LABELS
    stringstream label;
    label << "Individual " << icounter1;
    outfile << Rcomment(ss.str().c_str()) << bclib::newline;
#endif
    lcounter1 = 1;
    for (locus_i = indiv_i->begin(); locus_i != indiv_i->end(); ++locus_i) {
      for (gt_i = locus_i->begin(); gt_i != locus_i->end(); ++gt_i) {
        SumProb = *gt_i;
	//average over iterations
        *gt_i *= numIterationsReciprocal;
        // Round the probabilities to 6 decimal places
        // Makes them sum up to exactly 1.0 after outputting
        // (assuming they do sum up to 1.0 before rounding)
        UnRoundedMeanProb = *gt_i;
        *gt_i = floor((*gt_i) * 1e6 + 0.5) * 1e-6;

        // Output the number, preceded by a comma if not the first element
        outfile << *gt_i;

#ifdef GPO_DEBUG
         if (isnan(*gt_i)) {
          cerr << "*gt_i is nan. "
              << SumProb << " * "
              << numIterationsReciprocal << " -> " 
              << UnRoundedMeanProb << " -> "
              << "nan" << endl;
          throw nanException("*gt_i is nan ", __FILE__, __LINE__);
        }
#endif
      }//end loop over genotypes
#ifdef GPO_LABELS
      if(icounter1< NumMaskedIndivs || lcounter1 < NumMaskedLoci){
	stringstream label;
	label << " masked locus no. "  << lcounter1;
	outfile << Rcomment(label.str().c_str());
      }
#endif
      ++lcounter1;
      // Keep one line per genotype triplet
      outfile << bclib::newline;
    }//end loop over individuals
    ++icounter1;
  }//end loop over loci

  //output dimensions
  std::vector<int> dim(3,0);
  dim[0] = 3;  
  dim[1] = NumMaskedLoci;
  dim[2] = NumMaskedIndivs;

  std::vector<std::string> dim1labels(dim[0],"");
  dim1labels[0] = "Genotype1";
  dim1labels[1] = "Genotype2";
  dim1labels[2] = "Genotype3";

  vector<vector<string> > dimnames;
  dimnames.push_back(dim1labels);
  dimnames.push_back(LocusLabels);
  outfile.close(dim, dimnames);
}
