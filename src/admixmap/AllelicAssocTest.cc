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
#include "bcppcl/TableWriter.h"

AllelicAssocTest::AllelicAssocTest(){
  locusObsIndicator = 0;
  options = 0;
  individuals = 0;
  NumOutputs = 0;
  test = false;
  numUpdates = 0;
  numPrintedIterations = 0;
  NumCompositeLoci = 0;
}

AllelicAssocTest::~AllelicAssocTest(){
  if(test)
    ROutput();
  delete[] locusObsIndicator;
}

void AllelicAssocTest::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci,
				  LogWriter& Log){
  test = true;
  options = op;
  individuals = indiv;
  Lociptr = Loci;
  chrm = Loci->getChromosomes();
  NumCompositeLoci = Loci->GetNumberOfCompositeLoci();
  Log.setDisplayMode(Quiet);

  const int L = Lociptr->GetNumberOfCompositeLoci();

  AllelicAssocRObject.open(options->getAllelicAssociationScoreFilename());      

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
    NumLoci.push_back(NumberOfLoci);
    
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
    //cannot test for SNPs in Haplotype if only simple loci
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){
      HaplotypeAssocRObject.open(options->getHaplotypeAssociationScoreFilename());

      NumMergedHaplotypes.assign(NumCompositeLoci, 1);
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

void AllelicAssocTest::MergeRareHaplotypes(const std::vector<double> &alpha0){
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
  NumMergedHaplotypes.clear();
  NumMergedHaplotypes.assign(NumCompositeLoci, 1);
  for(unsigned int j = 0; j < NumCompositeLoci; j++ ){
    if( NumLoci[j] > 1 ){//skip simple loci
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( alphaScaled);
      NumMergedHaplotypes[j] = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();

      HaplotypeAssocTests[j]->Resize(NumMergedHaplotypes[j]);
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

void AllelicAssocTest::Update( const Individual* const ind, double YMinusEY, double phi, 
			       double DInvLink, bool missingOutcome)
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
	  HaplotypeAssocTests[locus]->Update(happair, (*Lociptr)(locus), 
					     ind->getAdmixtureProps(), YMinusEY, phi, DInvLink);
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


void AllelicAssocTest::Output(const Vector_s& LocusLabels){

  ++numPrintedIterations;

  for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(options->getHapMixModelIndicator() || locusObsIndicator[j]){
      SubTests[j]->Output(AllelicAssocRObject, LocusLabels[j], Lociptr->GetLocus(j), false, numUpdates);

    }//end j loop over comp loci


  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){

    const string dummystring;
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) 
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	HaplotypeAssocTests[j]->Output(HaplotypeAssocRObject, dummystring, Lociptr->GetLocus(j), false, numUpdates);
    }
  }//end if hap assoc test
}

void AllelicAssocTest::WriteFinalTables(const Vector_s& LocusLabels, LogWriter& Log){
  {
  string filename(options->getResultsDir());
  filename.append("/AllelicAssocTestsFinal.txt");
  TableWriter finaltable(filename.c_str());
  Log << Quiet << "Tests for allelic association written to " << filename << "\n";
  finaltable << "Locus\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare"
	     << newline;

  for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(options->getHapMixModelIndicator() || locusObsIndicator[j]){
       SubTests[j]->Output(finaltable, LocusLabels[j], Lociptr->GetLocus(j), true, numUpdates);
    }//end j loop over comp loci

  finaltable.close();
  }

  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    
    string filename(options->getResultsDir());
    filename.append("/HaplotypeAssocTestsFinal.txt");
    
    TableWriter finaltable(filename.c_str());
    Log << Quiet << "Tests for haplotype association written to " << filename << "\n";
    finaltable << "Locus\tHaplotype\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare"
	       << newline;

    const string dummystring;
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) 
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	HaplotypeAssocTests[j]->Output(finaltable, dummystring, Lociptr->GetLocus(j), true, numUpdates);
      
    }
    finaltable.close();
  }//end if hap assoc test
}

void AllelicAssocTest::ROutput(){
   int count = 0;
   for(unsigned int j = 0; j < NumCompositeLoci; j++ )
    if(locusObsIndicator[j]){
    if( NumLoci[j] == 1 ) count += SubTests[j]->getDim();
    else count += NumLoci[j];
  }         

  vector<vector<string> >labels(1);
  labels[0].push_back( "Locus");
  labels[0].push_back( "log10Pvalue");

  vector<int> dimensions(3,0);
  dimensions[0] = labels[0].size();
  dimensions[1] = count;
  dimensions[2] = (int)(numPrintedIterations);

  AllelicAssocRObject.close(dimensions, labels);

  /** 
   * writes out the dimensions and labels of the        
   * R-matrix for score tests of SNPs in haplotypes.        
   */ 
  if( options->getTestForHaplotypeAssociation()  ){      
    count = 0;
    for(unsigned int j = 0; j < NumCompositeLoci; j++ ){
      if( NumLoci[j] > 1 ){
	count += NumMergedHaplotypes[j];
      }
    }         

    vector<vector<string> > labels(1);
    labels[0].push_back("Locus");
    labels[0].push_back("Haplotype");
    labels[0].push_back("log10PValue");

    vector<int> dimensions(3,0);
    dimensions[0] = labels[0].size();
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);

    HaplotypeAssocRObject.close(dimensions,labels);
  }
 }


