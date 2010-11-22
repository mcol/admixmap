//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file HWTest.cc
/// Implementation of the HWTest class.
//=============================================================================

#include "HWTest.h"
#include "gsl/gsl_cdf.h"
#include "bclib/linalg.h"
#include "bclib/DelimitedFileWriter.h"
#include "bclib/LogWriter.h"

HWTest::HWTest(){
  score = 0;
  sumscore = 0;
  sumscore2 = 0;
  suminfo = 0;
  NumInd = 0;
  NumLoci = 0;
  numUpdates  = 0;
  test=false;
}

HWTest::~HWTest(){
  //free_matrix(sumscore, NumInd);
  //free_matrix(sumscore2, NumInd);
  //free_matrix(suminfo, NumInd);
  delete[] score;
  delete[] sumscore;
  delete[] sumscore2;
  delete[] suminfo;
}

//void HWTest::Initialise(Options *options, int nind, int nloci, LogWriter *Log){
void HWTest::Initialise(int nloci){

  //NumInd = nind;
  NumLoci = nloci;
  test=true;

  //sumscore = alloc2D_d(NumInd, NumLoci);
  //sumscore2 = alloc2D_d(NumInd, NumLoci);
  //suminfo = alloc2D_d(NumInd, NumLoci);
  score = new double[NumLoci];
  sumscore = new double[NumLoci];
  sumscore2 = new double[NumLoci];
  suminfo = new double[NumLoci];
  for(int i = 0; i < NumLoci; ++i){
    sumscore[i] = 0.0;
    sumscore2[i] = 0.0;
    suminfo[i] = 0.0;
  }
}    

///reset score
void HWTest::Reset(){
  for(int j = 0; j < NumLoci; ++j){
    score[j] = 0.0;
  }
}

void HWTest::Update(const IndividualCollection* const IC, const Genome* const Loci){
  double H; // prob heterozygous
  bool h;
  int locus = 0; // absolute simple locus number within loop over individuals
  int slocus = 0; // incremented with loop over composite loci 
  int complocus = 0;  // absolute comp locus number
  int Ancestry0, Ancestry1;
  double **Prob0 = 0, **Prob1 = 0;
  const unsigned numberOfChromosomes = Loci->GetNumberOfChromosomes();
  const unsigned *sizeOfChromosome   = Loci->GetSizesOfChromosomes();

  Reset();

  for (unsigned int chr = 0; chr < numberOfChromosomes; ++chr) {

    // loop over composite loci on chromosome
    for (unsigned int j = 0; j < sizeOfChromosome[chr]; ++j) {
      const CompositeLocus & cLoc = (*Loci)[ complocus ];
      const unsigned int NumSimpleLociWithinCL = cLoc.GetNumberOfLoci();
      int alleles0[ NumSimpleLociWithinCL ];
      int alleles1[ NumSimpleLociWithinCL ];
      //allocate arrays to hold marginal alleleprobs; could be done in function but easier to control here.      
      Prob0 = new double*[NumSimpleLociWithinCL];
      Prob1 = new double*[NumSimpleLociWithinCL];
      for(unsigned int jj = 0; jj < NumSimpleLociWithinCL; ++jj){
	Prob0[jj] = new double[ cLoc.GetNumberOfAllelesOfLocus(jj) ];
	Prob1[jj] = new double[ cLoc.GetNumberOfAllelesOfLocus(jj) ];
      }
      
      for(int i = 0; i < IC->getSize(); ++i) { // loop over individuals to get locus ancestry
	const PedBase & ind = IC->getElement(i);
	//bool diploid = !ind.isFemale() && Loci->isXLocus(complocus);
	bool diploid = !(ind.isHaploidatLocus(complocus));
	if(diploid){//skip Xloci in males
	  Ancestry0 = ind.GetLocusAncestry( (int)chr, 0, j);
	  Ancestry1 = ind.GetLocusAncestry( (int)chr, 1, j);
	  const int* happair = ind.getSampledHapPair(complocus);

	  cLoc.HaplotypeSetter.decodeIntAsHapAlleles(happair[0], alleles0);

	  cLoc.HaplotypeSetter.decodeIntAsHapAlleles(happair[0], alleles0);
	  cLoc.HaplotypeSetter.decodeIntAsHapAlleles(happair[1], alleles1);
	  //retrieve marginal alleleprobs for composite locus complocus, given current ancestry states on each gamete of ind
	  cLoc.getLocusAlleleProbs(Prob0, Ancestry0);// #loci x #alleles array  
	  cLoc.getLocusAlleleProbs(Prob1, Ancestry1);

	  locus = slocus;
	  
	  for(unsigned int jj = 0; jj < NumSimpleLociWithinCL; ++jj){       //loop over simple loci within comp locus
	    if( !ind.simpleGenotypeIsMissing(locus)){ //non-missing genotype, assumes second gamete missing if first is
	      h = alleles0[jj] != alleles1[jj];
	      H = 1.0;
	      //compute prob of heterozygosity by subtracting from 1 the prob of homozygosity, ie sum of diagonal products	    
	      for(int a = 0; a < cLoc.GetNumberOfAllelesOfLocus(jj); ++a){//loop over alleles
		H -= Prob0[jj][a] * Prob1[jj][a];
	      }
	      //accumulate score over individuals
	      if( h ){//heterozygous - prob H under null
		score[locus ] -= 1.0 - H; 
	      }
	      else {//homozygous - prob (1-H) under null
		score[locus ] += H;   
	      }
	      suminfo[locus ] += H * (1.0 - H); 
	    } 
	    ++locus;
	  } // ends loop over simple loci within compound locus
	}
      } // ends loop over individuals 

      //reset pointers for next compound locus	
      bclib::free_matrix(Prob0, NumSimpleLociWithinCL);
      bclib::free_matrix(Prob1, NumSimpleLociWithinCL);
      slocus += NumSimpleLociWithinCL;
      ++complocus;
    }//end loop over compound loci
  } // end loop over chromosomes
  for(int j = 0; j < NumLoci; ++j){
    sumscore[j] += score[j];
    sumscore2[j] += score[j] * score[j];
  }
  ++numUpdates;
}

void HWTest::Output(const char *filename, const Vector_s& LocusLabels,
                    bclib::LogWriter& Log) {
  try{
    bclib::DelimitedFileWriter outputfile(filename);
    Log << bclib::Quiet << "Tests for Hardy-Weinberg equilibrium writen to " << filename << "\n";
    //header line
    outputfile <<"Locus\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tz-score\tp-value" << bclib::newline;

    for(int j = 0; j < NumLoci; j++ ){
      //call function in base class
      OutputScalarScoreTest(numUpdates, outputfile, LocusLabels[j], sumscore[j], sumscore2[j], suminfo[j], true);
    }
  }
  //if something goes wrong just exit
  catch(string s){
    Log << s;
    return;
  }
  catch(...){
    return;
  }
}

