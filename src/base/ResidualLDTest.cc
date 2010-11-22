//=============================================================================
//
// Copyright (C) 2006, 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file ResidualLDTest.cc
/// Implementation of the ResidualLDTest class.
//=============================================================================

#include "ResidualLDTest.h"
#include "FreqArrays.h"
#include "IndividualCollection.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_errno.h"
#include "bclib/linalg.h"//for HH_solve to compute chi-sq
#include "bclib/LogWriter.h"
#include "bclib/DelimitedFileWriter.h"
#include <cmath>

ResidualLDTest::ResidualLDTest(){
  individuals = 0;
  numUpdates = 0;
  numPrintedIterations = 0;
  NumIntervals = 0;
}

ResidualLDTest::~ResidualLDTest(){
  //Tcount.clear();
  if(test){
    
    vector<vector<string> >labels(1);
    labels[0].push_back("Loci");
    labels[0].push_back("log10Pvalue");

    vector<int> dimensions(3,0);
    dimensions[0] = labels[0].size();
    dimensions[1] = NumIntervals;
    dimensions[2] = (int)(numPrintedIterations);

    //R_output3DarrayDimensions(&outputfile, dimensions, labels);
    R.close(dimensions, labels);
  }
}

void ResidualLDTest::Initialise(const char* filename, const IndividualCollection* const indiv, const Genome* const Loci){
  individuals = indiv;
  chrm = Loci->getChromosomes();
  Lociptr = Loci;
  test = true;
  NumIntervals = Loci->GetNumberOfCompositeLoci() - Loci->GetNumberOfChromosomes();

  const int NumOfChromosomes = Lociptr->GetNumberOfChromosomes();

  //open cumulative output file
  R.open(filename);

  SumScore.resize(NumOfChromosomes);
  SumScore2.resize(NumOfChromosomes);
  SumInfo.resize(NumOfChromosomes);

  Score.resize(NumOfChromosomes);
  Info.resize(NumOfChromosomes);
  
  int locus = 0;
  for(int j = 0; j < NumOfChromosomes; ++j){
    unsigned NumberOfLoci = Lociptr->GetSizeOfChromosome(j);
    
    SumScore[j].resize(NumberOfLoci-1);
    SumScore2[j].resize(NumberOfLoci-1);
    SumInfo[j].resize(NumberOfLoci-1);
    
    Score[j].resize(NumberOfLoci-1);
    Info[j].resize(NumberOfLoci-1);
    
    for(unsigned k = 0; k < NumberOfLoci-1; ++k){
      unsigned dim = (Lociptr->GetNumberOfStates(locus)-1) * (Lociptr->GetNumberOfStates(locus+1)-1);
      
      SumScore[j][k].assign(dim, 0.0);
      SumScore2[j][k].assign(dim*dim, 0.0);
      SumInfo[j][k].assign(dim*dim, 0.0);
      //Tcount.push_back(0);
      
      Score[j][k].assign(dim, 0.0);
	Info[j][k].assign(dim*dim, 0.0);
	++locus;
    }
    ++locus;//for last locus on chrm
  }
}

void ResidualLDTest::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.
  if(test){
    int locus = 0;
    const unsigned numberOfChromosomes = Lociptr->GetNumberOfChromosomes();
    for (unsigned j = 0; j < numberOfChromosomes; ++j) {
      for(unsigned k = 0; k < Lociptr->GetSizeOfChromosome(j)-1; ++k){
	fill(Score[j][k].begin(), Score[j][k].end(), 0.0);
	fill(Info[j][k].begin(), Info[j][k].end(), 0.0);
	++locus;
      }
      ++locus;//for the last locus on chromosome
    }
  }
}

void ResidualLDTest::Update(const FreqArray& AlleleFreqs, bool ){
  int abslocus = 0;
  const unsigned numberOfChromosomes = Lociptr->GetNumberOfChromosomes();
  const unsigned *sizeOfChromosome = Lociptr->GetSizesOfChromosomes();

  for (unsigned c = 0; c < numberOfChromosomes; ++c) {
    for (unsigned j = 0; j < sizeOfChromosome[c] - 1; ++j) {
      //if(ishapmixmodel)
      //UpdateScoresForResidualAllelicAssociation2(c, j, AlleleFreqs[abslocus], AlleleFreqs[abslocus+1]);
      //else
      UpdateScoresForResidualAllelicAssociation(c, j, AlleleFreqs[abslocus], AlleleFreqs[abslocus+1]);
      ++abslocus;
    }
    ++abslocus;//for last locus on chrm
  }

  //vector<unsigned>::iterator T_iter = Tcount.begin();
  //accumulate score, square of score and info
  for (unsigned c = 0; c < numberOfChromosomes; ++c)
    for (unsigned k = 0; k < sizeOfChromosome[c] - 1; ++k) {
      int locus = chrm[c]->GetLocus(k);
      unsigned dim = (Lociptr->GetNumberOfStates(locus)-1) * (Lociptr->GetNumberOfStates(locus+1)-1);

// 	if(ishapmixmodel){//assuming all SNPs in hapmixmodel
// 	  const double Tobs = Score[c][k][0];
// 	  const double Trep = Rand::gennor(0.0, Info[c][k][0]);
// 	  if(fabs(Trep)> fabs(Tobs))++(*T_iter);
// 	  ++T_iter;
// 	}

      for(unsigned j = 0; j < dim; ++j){
	SumScore[c][k][j] += Score[c][k][j];
	for(unsigned jj = 0; jj < dim; ++jj){
	  SumScore2[c][k][j*dim +jj] += Score[c][k][j]*Score[c][k][jj];
	  SumInfo[c][k][j*dim +jj] += Info[c][k][j*dim + jj];
	}
      }
    }

  ++numUpdates;
}


static int delta(int i, int j){
  return (i == j);
}

/**
 // M alleles at locus A, N alleles at locus B
 // notation is for arrays subscripted to run from 1 to M, 1 to N
 // r_mn is indicator variable for haplotype mn
 // phi_mn is expected frequency of haplotype mn given ancestry states at A, B
 // log-linear model: log (lambda_mn = hazard rate) = alpha_m + beta_n + gamma_mn
 // alpha_M = 1 - sum_{alpha_1 to alpha_{M-1}), similarly for beta, gammaM, gammaN
 //
 // log-likelihood is sum of terms of form r_mn * log(lambda_mn) - lambda_mn 
 // test null hypothesis that gamma_{11} to gamma_{M-1, N-1} are all 0
 // score U_mn = r_mn - r_mN - r_Mn + r_MN - (phi_mn - phi_mN - phi_Mn + phi_MN) 
 // info V_mnm'n' = delta_mm' * delta_nn' * phi_mn 
 //               + delta_mm' * phi_mN 
 //               + delta_nn' * phi_Mn 
 //               + phi_MN
 // 
 // r_mn = delta(hA[g], m) * delta(hB[g], n)
 // phi_mn = AlleleFreqsA[ancA[g]*M + m] * AlleleFreqsB[ancB[g]*N + n]
 
 // info depends on locus ancestry but not on indicator variables: could calculate matrix for all possible 
 // 2-locus ancestry states, then sum over counts of ancestry states
 */
void ResidualLDTest::UpdateScoresForResidualAllelicAssociation(int c, int locus,  
							       const double* const AlleleFreqsA, 
							       const double* const AlleleFreqsB) {
  int abslocus = chrm[c]->GetLocus(locus);//number of this locus
  int M = Lociptr->GetNumberOfStates(abslocus);
  int N = Lociptr->GetNumberOfStates(abslocus+1);

  int dim = (M-1)*(N-1);
  int ancA[2];//ancestry at A
  int ancB[2];//ancestry at B
  for(int i = 0; i < individuals->getSize(); i++) {
    const PedBase & ind = individuals->getElement(i);
    if( !ind.GenotypeIsMissing(abslocus) && !ind.GenotypeIsMissing(abslocus+1) ) {
      //skip missing genotypes as hap pairs not sampled
      ind.GetLocusAncestry(c, locus, ancA);
      ind.GetLocusAncestry(c, locus+1, ancB);
      const int* hA = ind.getSampledHapPair(abslocus);//realized hap pair at locus A: values from 0 to M-1
      const int* hB = ind.getSampledHapPair(abslocus+1);//realized hap pair at locus B: values from 0 to N-1
      int numGametes = 2;
      if(hA[1] < 0) numGametes = 1;//haplotypes are coded as happair with second hap as -1
      for(int g = 0; g < numGametes; ++g) {
	if(dim == 1) { // diallelic version
	  int h = (hA[g] == hB[g]);//indicator for coupling: 1 if 2-locus haplotype is 1-1 or 2-2
	  double phiA = AlleleFreqsA[ancA[g]*2];//frequency of first allele in this gamete's ancestry state at locus A
	  double phiB = AlleleFreqsB[ancB[g]*2];//   "           "    "         "     "         "       "       "    B
	  double ProbCoupling = phiA*phiB + (1 - phiA)*(1 - phiB);
	  Score[c][locus][0] += 2.0 * (h - ProbCoupling); 
	  Info[c][locus][0] += 4.0 * ProbCoupling * (1.0 - ProbCoupling);
	} else { // multi-allelic version
	  for(int m = 0; m < M-1; ++m) {   //update score vector of length dim = (M-1)(N-1)
	    for(int n = 0; n < N-1; ++n) { //m and n index elements of score vector, rows of info matrix
	      Score[c][locus][m*(N-1)+n] +=  delta(hA[g], m) * delta(hB[g], n) // r
		-                        delta(hA[g], m) * delta(hB[g], N-1) 
		-                        delta(hA[g], M-1) * delta(hB[g], n) 
		+                        delta(hA[g], M-1) * delta(hB[g], N-1) // observed
		- AlleleFreqsA[ancA[g]*M + m] * AlleleFreqsB[ancB[g]*N + n]
		+ AlleleFreqsA[ancA[g]*M + m] * AlleleFreqsB[ancB[g]*N + N-1]
		+ AlleleFreqsA[ancA[g]*M + M-1] * AlleleFreqsB[ancB[g]*N + n]
		- AlleleFreqsA[ancA[g]*M + M-1] * AlleleFreqsB[ancB[g]*N + N-1]; // minus expected
	      for(int mp = 0; mp < M-1 ; ++mp) {  //update info
		for(int np = 0; np < N-1; ++np) { //mp and np index cols of info matrix
		  Info[c][locus][(m*(N-1)+n) * dim + (mp*(N-1)+np)] += 
		    delta(m,mp) * delta(n, np) * AlleleFreqsA[ancA[g]*M + m]   * AlleleFreqsB[ancB[g]*N + n]
		    +             delta(m, mp) * AlleleFreqsA[ancA[g]*M + m]   * AlleleFreqsB[ancB[g]*N + N-1]
		    +             delta(n, np) * AlleleFreqsA[ancA[g]*M + M-1] * AlleleFreqsB[ancB[g]*N + n]
		    +                            AlleleFreqsA[ancA[g]*M + M-1] * AlleleFreqsB[ancB[g]*N + N-1];
		}
	      }
	    }
	  }
	} // end else block for multi-allelic version
      } //end loop over gametes
    } // end block conditional on non-missing genotypes
  } //end loop over individuals
}

///this version for Bayesian p-values (SNPs only)
void ResidualLDTest::UpdateScoresForResidualAllelicAssociation2(int c, int locus,  
							       const double* const AlleleFreqsA, 
							       const double* const AlleleFreqsB) {
  int abslocus = chrm[c]->GetLocus(locus);//number of this locus
  int M = Lociptr->GetNumberOfStates(abslocus);
  int N = Lociptr->GetNumberOfStates(abslocus+1);

  int dim = (M-1)*(N-1);
  int ancA[2];//ancestry at A
  int ancB[2];//ancestry at B
  for(int i = 0; i < individuals->getSize(); i++) {
    const PedBase & ind = individuals->getElement(i);
    if( !ind.GenotypeIsMissing(abslocus) && !ind.GenotypeIsMissing(abslocus+1) ) {
      //skip missing genotypes as hap pairs not sampled
      ind.GetLocusAncestry(c, locus, ancA);
      ind.GetLocusAncestry(c, locus+1, ancB);
      const int* hA = ind.getSampledHapPair(abslocus);//realized hap pair at locus A: values from 0 to M-1
      const int* hB = ind.getSampledHapPair(abslocus+1);//realized hap pair at locus B: values from 0 to N-1
      int numGametes = 2;
      if(hA[1] < 0) numGametes = 1;//haplotypes are coded as happair with second hap as -1
      for(int g = 0; g < numGametes; ++g) {
	if(dim == 1) { // diallelic version
	  double phiA = AlleleFreqsA[ancA[g]*2];//frequency of first allele in this gamete's ancestry state at locus A
	  double phiB = AlleleFreqsB[ancB[g]*2];//   "           "    "         "     "         "       "       "    B

	  if(hA[g] == 0){
	    if(hB[g]==0)Score[c][locus][0] += (1.0 - phiA) * (1.0 - phiB) ;//1.0 / (phiA*phiB);
	    else Score[c][locus][0] -= (1.0 - phiA) * phiB;//1.0 / (phiA*phiB);
	  }
	  else if(hA[g]== 1){
	    if(hB[g]==0)Score[c][locus][0] -= phiA * (1.0 - phiB);//1.0 / ((1.0-phiA)*(1.0-phiB));
	    else Score[c][locus][0] += phiA * phiB;//1.0 / ((1.0-phiA)*(1.0-phiB));
	  }
	
	  Info[c][locus][0] += phiA*phiB*(1.0-phiA)*(1.0-phiB);
	}
      }
    }
  }
}

void ResidualLDTest::Output(const std::vector<std::string>& LocusLabels){
  if(test){
    OutputTestsForResidualAllelicAssociation(R, false, LocusLabels);
    ++numPrintedIterations;
  }
}

void ResidualLDTest::WriteFinalTable(const char* filename, const std::vector<std::string>& LocusLabels, bclib::LogWriter& Log){
  if(test){
    bclib::DelimitedFileWriter finaltable(filename);
    Log << bclib::Quiet << "Tests for residual allelic association" << " written to " 
 	<< filename << "\n";

    finaltable << "Loci" << "Score" << "CompleteInfo" << "ObservedInfo" << "PercentInfo" << "df" << "ChiSquared" << "PValue" << bclib::newline;
    
    OutputTestsForResidualAllelicAssociation(finaltable, true, LocusLabels);
    finaltable.close();
  }
}

void ResidualLDTest::OutputTestsForResidualAllelicAssociation(bclib::DelimitedFileWriter& outputstream, bool final,
							      const std::vector<std::string>& LocusLabels){
  //cannot function in base class as it output a line for each dimension of score
  double *score = 0, *ObservedInfo = 0;

  int abslocus = 0;
  //vector<unsigned>::iterator T_iter = Tcount.begin();
  for(unsigned int c = 0; c < Lociptr->GetNumberOfChromosomes(); c++ ){
    for(unsigned j = 0; j < Lociptr->GetSizeOfChromosome(c)-1; ++j){
      int M = Lociptr->GetNumberOfStates(abslocus)-1;
      int N = Lociptr->GetNumberOfStates(abslocus+1)-1;

      const string label1 = LocusLabels[abslocus];
      const string label2 = LocusLabels[abslocus+1];

      int dim = M*N;
      score = new double[dim];
      ObservedInfo = new double[dim*dim];
      double obsinfo = 0.0;
      double compinfo = 0.0;

      for(int i = 0; i < dim; ++i){
	score[i] = SumScore[c][j][i]/(double) numUpdates; //score
      }
      for(int i = 0; i < dim; ++i){
	for(int ii = 0; ii < dim; ++ii){
	  ObservedInfo[i*dim +ii] = SumInfo[c][j][i*dim+ii]/ (double) numUpdates//complete info
	    + score[i]*score[ii] - SumScore2[c][j][i*dim+ii]/(double)numUpdates;//-missing info
	}
	obsinfo += ObservedInfo[i*dim +i];//trace of Observed Info
	compinfo += SumInfo[c][j][i*dim+i]/ (double) numUpdates;//trace of Complete Info
      }

      //output labels
      const string label = "\"" + label1 + "/" + label2 + "\"";
      outputstream << label;
      //output 3 decimal places
      outputstream.setDecimalPrecision(3);
      if(final)outputstream << score[0] 
			    << compinfo 
			    << obsinfo
			    << double2R((obsinfo*100)/compinfo) 
			    << dim;

      //outputstream << (float)*T_iter /(float)iterations << newline;
      //++T_iter; 

      //compute chi-squared statistic
      try{
	double* VinvU = new double[dim];
	bclib::HH_solve(dim, ObservedInfo, score, VinvU);
	double chisq = 0.0;
	for(int i = 0; i < dim; ++i){
	  chisq += score[i] * VinvU[i];
	}
	delete[] VinvU;
	
	if(chisq < 0.0){
	  outputstream << "NA" ;
	  if(final) outputstream << "NA";
	  outputstream << bclib::newline;
	}
	else {
	  //compute p-value
	  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
	  try{
	    double pvalue = gsl_cdf_chisq_Q (chisq, dim);
	    if(final)outputstream << double2R(chisq) << double2R(pvalue)<< bclib::newline;
	    else outputstream << double2R(-log10(pvalue)) << bclib::newline;
	  }
	  catch(...){
	    if(final)outputstream << double2R(chisq) << "NA" << bclib::newline;
	    else outputstream << "NA" << bclib::newline;
	  }
	  gsl_set_error_handler(old_handler);
	}
      }
      catch(...){//in case ObservedInfo is rank deficient
	outputstream  << "NA";
	if(final)outputstream << dim << "NA" << "NA" ;
	outputstream << bclib::newline;
      }

      delete[] score;
      delete[] ObservedInfo;
      ++abslocus;
    }//end loop over loci on chromosome
    ++abslocus;//for last locus on chromosome
  }
}
