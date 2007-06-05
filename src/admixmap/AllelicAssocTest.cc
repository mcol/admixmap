/** 
 *   ADMIXMAP
 *   AllelicAssocTest.cc 
 *   Class implements 

 *   (1) Score test for allelic association
 *   (2) Score test for within-halpotype association
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AllelicAssocTest.h"

AllelicAssocTest::AllelicAssocTest(){
  locusObsIndicator = 0;

  options = 0;
  individuals = 0;

  NumOutputs = 0;
  onFirstLineAllelicAssoc = true;
  onFirstLineHapAssoc = true;

  test = false;
  numUpdates = 0;
  numPrintedIterations = 0;
}

AllelicAssocTest::~AllelicAssocTest(){
  delete[] locusObsIndicator;
}

void AllelicAssocTest::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci, 
			    LogWriter &Log){
  test = true;
  options = op;
  individuals = indiv;
  chrm = Loci->getChromosomes();
  Lociptr = Loci;
  Log.setDisplayMode(Quiet);

  const int L = Lociptr->GetNumberOfCompositeLoci();

  OpenFile(Log, &outputfile, options->getAllelicAssociationScoreFilename(), 
	   "Tests for allelic association", true);
      
      locusObsIndicator = new int[L];
      
      //search for loci with no observed genotypes
      for(int j = 0; j < L; ++j){
	locusObsIndicator[j] = false;
	for(int i = indiv->getFirstScoreTestIndividualNumber(); i < indiv->getNumberOfIndividualsForScoreTests(); i++){
	  if(!indiv->getIndividual(i)->GenotypeIsMissing(j)){
	    locusObsIndicator[j] = true;
	  }
	}
      }
    
      unsigned NumCovars = indiv->GetNumCovariates() - indiv->GetNumberOfInputCovariates();
      AllelicAssocSubTest::SetNumCovars(NumCovars);

      for( int j = 0; j < L; j++ ){
	const int NumberOfStates = Lociptr->GetNumberOfStates(j);
	const int NumberOfLoci = Lociptr->getNumberOfLoci(j);
	
	if(NumberOfLoci > 1 )//haplotype
	  SubTests.push_back(new WithinHaplotypeTest(NumberOfLoci, true));
	else if(NumberOfStates == 2 )//simple diallelic locus
	  SubTests.push_back(new SNPTest(true));
	else//simple multiallelic locus
	  SubTests.push_back(new MultiAllelicLocusTest(NumberOfStates, true));
	

      }//end loop over loci

  /*----------------------
    | haplotype association |
    ----------------------*/  
  if( strlen( options->getHaplotypeAssociationScoreFilename() ) ){
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){//cannot test for SNPs in Haplotype if only simple loci
      OpenFile(Log, &HaplotypeAssocScoreStream, options->getHaplotypeAssociationScoreFilename(), 
	       "Tests for haplotype associations", true);
      for( int j = 0; j < L; j++ ){

	if( Lociptr->getNumberOfLoci(j) > 1 )
	  HaplotypeAssocTests.push_back(new HaplotypeTest(1, true));
	else
	  HaplotypeAssocTests.push_back(0);
      }
    }
    else {
      op->setTestForHaplotypeAssociation(false);
      Log << "ERROR: Cannot test for haplotype associations if all loci are simple\n" << "This option will be ignored\n";
    }
  }
}

void AllelicAssocTest::SetAllelicAssociationTest(const std::vector<double> &alpha0){
  /*
    Invokes merging of haplotypes in CompositeLocus and resizes arrays for allelic association test accordingly  
    alpha0 = alpha[0], pop admixture dirichlet params, from Latent
  */

  //first scale alphas so they sum to 1
  const unsigned K = options->getPopulations();
  double* alphaScaled = new double[K];
  double sum  = accumulate(alpha0.begin(), alpha0.end(), 0.0, std::plus<double>());//sum of alpha0 over pops
  for( int k = 0; k < options->getPopulations(); k++ )
    alphaScaled[k] = alpha0[k] / sum;

  //merge rare haplotypes
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){//skip simple loci
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( alphaScaled);
      const unsigned NumMergedHaplotypes = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();

      HaplotypeAssocTests[j]->Resize(NumMergedHaplotypes);
    }
  }

  delete[] alphaScaled;
}

void AllelicAssocTest::Reset(){
  if( test ){
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      SubTests[j]->Reset();
      if( options->getTestForHaplotypeAssociation() && Lociptr->getNumberOfLoci(j) > 1 ){
	HaplotypeAssocTests[j]->Reset();
      }
    }
  }
}

void AllelicAssocTest::Update( const Individual* const ind, double YMinusEY, double phi, double DInvLink, bool missingOutcome)
{
  int locus = 0;

  for(unsigned int j = 0; j < Lociptr->GetNumberOfChromosomes(); j++ ){
    for(unsigned int jj = 0; jj < Lociptr->GetSizeOfChromosome(j); jj++ ){
      //skip loci with missing genotypes as hap pairs have not been sampled for these
      bool condition = true;

      if(options->getHapMixModelIndicator())condition = !missingOutcome;
      else condition = !ind->GenotypeIsMissing(locus);
      //in hapmixmodel, skip individuals with observed outcome
      if(condition){
	//retrieve sampled hap pair from Individual
	const int* happair = ind->getSampledHapPair(locus);
	const unsigned numLoci = Lociptr->getNumberOfLoci(locus);
	
	SubTests[locus]->Update(happair, (*Lociptr)(locus), ind->getAdmixtureProps(), YMinusEY, phi, DInvLink);
	  
	if( options->getTestForHaplotypeAssociation() && numLoci > 1){
	  HaplotypeAssocTests[locus]->Update(happair, (*Lociptr)(locus), ind->getAdmixtureProps(), YMinusEY, phi, DInvLink);
	}

    }
    locus++;
    }
  }
}


void AllelicAssocTest::Accumulate(){
  try{
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( test ){
	SubTests[j]->Accumulate();	

	if( options->getTestForHaplotypeAssociation() && Lociptr->getNumberOfLoci(j) > 1){
	  HaplotypeAssocTests[j]->Accumulate();
	}
      }
      
    }//end comp locus loop
  }
  catch(string s){
    string error_string = "Error accumulating scores for allelicassociation or haplotype association scoretest\n";
    error_string.append(s);
    throw(error_string);
  }
  ++numUpdates;
}


void AllelicAssocTest::Output(const Vector_s& LocusLabels, bool final){

  string sep = final ? "\t" : ",";//separator
  ofstream* outfile;
  if(!final)++numPrintedIterations;

  if(final){
    string filename(options->getResultsDir());
    filename.append("/AllelicAssocTestsFinal.txt");
    outfile = new ofstream(filename.c_str(), ios::out);
    *outfile << "Locus\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare\n";
  }
  else {
    outfile = &outputfile;
  }
  for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(options->getHapMixModelIndicator() || locusObsIndicator[j]){
      const string locuslabel = LocusLabels[j];
      
      SubTests[j]->Output(outfile, locuslabel, Lociptr->GetLocus(j), final, onFirstLineAllelicAssoc, numUpdates);

      onFirstLineAllelicAssoc = false;
    }//end j loop over comp loci
  if(final)delete outfile;

  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    if(final){
      string filename(options->getResultsDir());
      filename.append("/HaplotypeAssocTestsFinal.txt");
      outfile = new ofstream(filename.c_str(), ios::out);
      *outfile << "Locus\tHaplotype\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare\n";
    }
    else outfile = &HaplotypeAssocScoreStream;

    const string dummystring;
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) 
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	HaplotypeAssocTests[j]->Output(outfile, dummystring, Lociptr->GetLocus(j), final, onFirstLineHapAssoc, numUpdates);
      
      onFirstLineHapAssoc = false;
    }
    if(final)delete outfile;
  }//end if hap assoc test
}

void AllelicAssocTest::ROutput(){
   int count = 0;
   for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(locusObsIndicator[j]){
    const int NumberOfLoci = (* Lociptr)(j)->GetNumberOfLoci();

    if( NumberOfLoci == 1 ) count += SubTests[j]->getDim();
    else count += NumberOfLoci;
  }         
  vector<int> dimensions(3,0);
  dimensions[0] = 2;
  dimensions[1] = count;
  dimensions[2] = (int)(numPrintedIterations);
  vector<string> labels(dimensions[0],"");
  labels[0] = "Locus";
  labels[1] = "log10Pvalue";
  R_output3DarrayDimensions(&outputfile, dimensions, labels);

  /** 
   * writes out the dimensions and labels of the        
   * R-matrix for score tests of SNPs in haplotypes.        
   */ 
  if( options->getTestForHaplotypeAssociation()  ){      
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	count += (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      }
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 3;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Haplotype";
    labels[2] = "log10PValue";
    R_output3DarrayDimensions(&HaplotypeAssocScoreStream,dimensions,labels);
  }
 }


