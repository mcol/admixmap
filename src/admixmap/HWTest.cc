/** 
 *   HWTest.cc 
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * Class implements a score test for deviation from Hardy-Weinberg equilibrium
 * in order to test for genotyping errors.
 * This version evaluates for each single locus and sums over individuals.
 * code for evaluation for each individual (rather than summing) is commented out.
 *
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HWTest.h"
//#include "IndividualCollection.h"
#include "gsl/gsl_cdf.h"
#include "utils/misc.h"

HWTest::HWTest(){
  score = 0;
  sumscore = 0;
  sumscore2 = 0;
  suminfo = 0;
  NumInd = 0;
  NumLoci = 0;
  samples  = 0;
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

//void HWTest::Initialise(AdmixOptions *options, int nind, int nloci, LogWriter *Log){
void HWTest::Initialise(const AdmixOptions* const options, int nloci, LogWriter &Log){

  //NumInd = nind;
  NumLoci = nloci;
  Log.setDisplayMode(Quiet);

  if( options->getHWTestIndicator() ){//not really necessary

    if ( strlen( options->getHWTestFilename() ) ){
      OpenFile(Log, &outputfile, options->getHWTestFilename(), "Tests for Hardy-Weinberg equilibrium", false);
	
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
    else{
      Log << "No hwtestfile given\n";
    }
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
  Individual *ind = 0;
  double **Prob0 = 0, **Prob1 = 0;

  Reset();

  for(unsigned int chr = 0; chr < Loci->GetNumberOfChromosomes(); ++chr){ //loop over chromosomes
    for(unsigned int j = 0; j < Loci->GetSizeOfChromosome(chr); ++j){ //loop over comp loci on chromosome
      unsigned NumSimpleLociWithinCL = Loci->getNumberOfLoci(complocus);
      int* alleles0 = new int[NumSimpleLociWithinCL];
      int* alleles1 = new int[NumSimpleLociWithinCL];
      //allocate arrays to hold marginal alleleprobs; could be done in function but easier to control here.      
      Prob0 = new double*[NumSimpleLociWithinCL];
      Prob1 = new double*[NumSimpleLociWithinCL];
      for(unsigned int jj = 0; jj < NumSimpleLociWithinCL; ++jj){
	Prob0[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];
	Prob1[jj] = new double[ (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj)];
      }
      
      for(int i = 0; i < IC->getSize(); ++i) { // loop over individuals to get locus ancestry
	ind = IC->getIndividual(i);
	//bool diploid = !ind->isFemale() && Loci->isXLocus(complocus);
	bool diploid = !(ind->isHaploidatLocus(complocus));
	if(diploid){//skip Xloci in males
	  Ancestry0 = ind->GetLocusAncestry( (int)chr, 0, j);
	  Ancestry1 = ind->GetLocusAncestry( (int)chr, 1, j);
	  const int* happair = ind->getSampledHapPair(complocus);
	  (*Loci)(complocus)->HaplotypeSetter.decodeIntAsHapAlleles(happair[0], alleles0);
	  (*Loci)(complocus)->HaplotypeSetter.decodeIntAsHapAlleles(happair[1], alleles1);
	  //retrieve marginal alleleprobs for composite locus complocus, given current ancestry states on each gamete of ind
	  (*Loci)(complocus)->getLocusAlleleProbs(Prob0, Ancestry0);// #loci x #alleles array  
	  (*Loci)(complocus)->getLocusAlleleProbs(Prob1, Ancestry1);
	  locus = slocus;
	  
	  for(unsigned int jj = 0; jj < NumSimpleLociWithinCL; ++jj){       //loop over simple loci within comp locus
	    if( !ind->simpleGenotypeIsMissing(locus)){ //non-missing genotype, assumes second gamete missing if first is
	      h = alleles0[jj] != alleles1[jj];
	      H = 1.0;
	      //compute prob of heterozygosity by subtracting from 1 the prob of homozygosity, ie sum of diagonal products	    
	      for(int a = 0; a < (*Loci)(complocus)->GetNumberOfAllelesOfLocus(jj); ++a){//loop over alleles
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
	  ind = 0;
	}
      } // ends loop over individuals 

      //reset pointers for next compound locus	
      free_matrix(Prob0, NumSimpleLociWithinCL);
      free_matrix(Prob1, NumSimpleLociWithinCL);
      delete[] alleles0;
      delete[] alleles1;
      slocus += NumSimpleLociWithinCL;
      ++complocus;
    }//end loop over compound loci
  } // end loop over chromosomes
  for(int j = 0; j < NumLoci; ++j){
    sumscore[j] += score[j];
    sumscore2[j] += score[j] * score[j];
  }
  ++samples;
}

void HWTest::Output(const Vector_s LocusLabels){
  //header line
  outputfile <<"Locus\tScore\tCompleteInfo\tMissingInfo\tObservedInfo\tPercentInfo\tz-score\tp-value"<<endl;

  for(int j = 0; j < NumLoci; j++ ){
    //call function in base class
    OutputScalarScoreTest(samples, &outputfile, LocusLabels[j], sumscore[j], sumscore2[j], suminfo[j], true);
  }
}

