//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AllelicAssocTest.cc
/// Implementation of the AllelicAssocTest class.
//=============================================================================

#include "AllelicAssocTest.h"
#include "bclib/DelimitedFileWriter.h"
#include "bclib/LogWriter.h"
#include "AdmixFilenames.h"
#include <numeric>   // for accumulate()


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
				  bclib::LogWriter& Log){
  test = true;
  options = op;
  individuals = indiv;
  Lociptr = Loci;
  chrm = Loci->getChromosomes();
  NumCompositeLoci = Loci->GetNumberOfCompositeLoci();
  Log.setDisplayMode(bclib::Quiet);

  const int L = Lociptr->GetNumberOfCompositeLoci();

  AllelicAssocRObject.open((options->getResultsDir() + "/" + ALLELICASSOCTEST_PVALUES).c_str());      

  locusObsIndicator = new int[L];
  
  //search for loci with no observed genotypes
  for(int j = 0; j < L; ++j){
    locusObsIndicator[j] = false;
	for(int i = indiv->getFirstScoreTestIndividualNumber(); i < indiv->getNumberOfIndividualsForScoreTests(); i++){
	  if ( ! indiv->getElement(i).GenotypeIsMissing(j) ) {
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
      SubTests.push_back(new WithinHaplotypeTest(NumberOfLoci));
    else if(NumberOfStates == 2 )//simple diallelic locus
      SubTests.push_back(new SNPTest());
    else//simple multiallelic locus
      SubTests.push_back(new MultiAllelicLocusTest(NumberOfStates));
    
    
  }//end loop over loci
  
  /*----------------------
    | haplotype association |
    ----------------------*/  
  if( options->getTestForHaplotypeAssociation() ){
    //cannot test for SNPs in Haplotype if only simple loci
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){
      const string hapassocfilename = options->getResultsDir() + "/" + HAPLOTYPEASSOCTEST_PVALUES;
      HaplotypeAssocRObject.open(hapassocfilename.c_str());

      NumMergedHaplotypes.assign(NumCompositeLoci, 1);
      for( int j = 0; j < L; j++ ){
	if( Lociptr->getNumberOfLoci(j) > 1 )
	  HaplotypeAssocTests.push_back(new HaplotypeTest(1));
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
    Invokes merging of haplotypes in CompositeLocus and resizes arrays for haplotype association test accordingly  
    alpha0 = alpha[0], pop admixture dirichlet params, from Latent
  */
  if( options->getTestForHaplotypeAssociation()){
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

void AllelicAssocTest::Update( const PedBase * const ind, double YMinusEY, double phi, 
			       double DInvLink, bool /*missingOutcome*/)
{
  int locus = 0;

  for(unsigned int j = 0; j < Lociptr->GetNumberOfChromosomes(); j++ ){

    bool isXChrom = Lociptr->isXChromosome(j);

    for(unsigned int jj = 0; jj < Lociptr->GetSizeOfChromosome(j); jj++ ){
      //skip loci with missing genotypes as hap pairs have not been sampled for these
      if(!ind->GenotypeIsMissing(locus)){
	//retrieve sampled hap pair from Individual
	const int* happair = ind->getSampledHapPair(locus);
	const unsigned numLoci = Lociptr->getNumberOfLoci(locus);

        SubTests[locus]->Update(happair, (*Lociptr)(locus),
                                ind->getAdmixtureProps(isXChrom),
                                YMinusEY, phi, DInvLink);

	if( options->getTestForHaplotypeAssociation() && numLoci > 1){
	  HaplotypeAssocTests[locus]->Update(happair, (*Lociptr)(locus), 
                                             ind->getAdmixtureProps(isXChrom),
                                             YMinusEY, phi, DInvLink);
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


void AllelicAssocTest::Output(){

  ++numPrintedIterations;

  for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(locusObsIndicator[j]){
      SubTests[j]->Output(AllelicAssocRObject, Lociptr->GetLocus(j), false, numUpdates);

    }//end j loop over comp loci


  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){

    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) 
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	HaplotypeAssocTests[j]->Output(HaplotypeAssocRObject, Lociptr->GetLocus(j), false, numUpdates);
    }
  }//end if hap assoc test
}

void AllelicAssocTest::WriteFinalTables(bclib::LogWriter& Log){
  {
  string filename = options->getResultsDir() + "/" + ALLELICASSOCTEST_FINAL;
  bclib::DelimitedFileWriter finaltable(filename.c_str());
  Log << bclib::Quiet << "Tests for allelic association written to " << filename << "\n";
  finaltable << "Locus" << "Score" << "CompleteInfo" << "ObservedInfo" 
	     << "PercentInfo" << "StdNormal" << "PValue" << "ChiSquare"
	     << bclib::newline;

  for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
    if(locusObsIndicator[j]){
      SubTests[j]->Output(finaltable, Lociptr->GetLocus(j), true, numUpdates);
    }//end j loop over comp loci

  finaltable.close();
  }

  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    
    string filename = options->getResultsDir() + "/" + HAPLOTYPEASSOCTEST_FINAL;
    bclib::DelimitedFileWriter finaltable(filename.c_str());
    Log << bclib::Quiet << "Tests for haplotype association written to " << filename << "\n";
    finaltable << "Locus" << "Haplotype" << "Score" << "CompleteInfo" << "ObservedInfo" 
	       << "PercentInfo" << "StdNormal" << "PValue" << "ChiSquare"
	       << bclib::newline;

    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) 
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	HaplotypeAssocTests[j]->Output(finaltable, Lociptr->GetLocus(j), true, numUpdates);
      
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


