/** 
 *   ADMIXMAP
 *   DispersionTest.cc 
 *   Class to implement test for dispersion of allele frequencies between the unadmixed populations sampled 
 *   and the corresponding ancestry-specific allele frequencies in the admixed population under study.  
 *   This is evaluated for each subpopulation at each locus, and as a global test over all loci.  
 *   Valid only if option priorallelefreqfile is specified. The results are "Bayesian p-values". 
 *   Copyright (c) 2005,2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "DispersionTest.h"
#include "AdmixOptions.h"
#include <numeric>
#include "bcppcl/linalg.h"
#include "bcppcl/dist.h"

using namespace ::std;

DispersionTest::DispersionTest(){
  divergentallelefreqstest = 0;
  NumberOfCompositeLoci = -1;
  NumberOfPopulations = 0;
  options = 0;
}

void DispersionTest::Initialise(const AdmixOptions* const op, LogWriter &Log, int NumLoci){
  options = op;
  NumberOfCompositeLoci = NumLoci;
  NumberOfPopulations = options->getPopulations();

  divergentallelefreqstest = alloc2D_i(NumberOfCompositeLoci + 1, NumberOfPopulations );
  
  if( options->getTestForDispersion() ){
    Log.setDisplayMode(Quiet);
    Log << "Writing dispersion test results to " << options->getDispersionTestFilename() << "\n";
    dispersionoutputstream.open( options->getDispersionTestFilename(), ios::out );
    if( !dispersionoutputstream ){
    Log.setDisplayMode(On);
    Log << "ERROR: Couldn't open dispersiontestfile\n";
   exit( 1 );
    }
  }
}

DispersionTest::~DispersionTest(){
  free_matrix(divergentallelefreqstest, NumberOfCompositeLoci + 1);
}

void DispersionTest::TestForDivergentAlleleFrequencies(const AlleleFreqs* const A, const IndividualCollection* const IC)
{
  vector<int> rep;
  vector<double> popfreqs;

  vector<int> AlleleCount; 
  vector<double> temp(NumberOfPopulations, 0.0);
  vector<vector<double> > LogLikelihood(NumberOfCompositeLoci + 1, temp);
  vector<vector<double> > RepLogLikelihood(LogLikelihood);
  temp.clear();
  vector<double> sum(NumberOfPopulations); vector<double> repsum(NumberOfPopulations);

  for( int k = 0; k < NumberOfPopulations; k++ ){
    sum[k] = 0.0; 
    repsum[k] = 0.0;
  }
  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    for( int k = 0; k < NumberOfPopulations; k++ ){
      popfreqs =  A->GetAlleleFreqs(j, k );

      // Generate replicate data conditional on locus ancestry
      //AlleleCount = A->GetAlleleCounts(j, k);
      AlleleCount = IC->getAlleleCounts(j, k, popfreqs.size());
      int sumcounts = accumulate(AlleleCount.begin(), AlleleCount.end(), 0, plus<int>());
      rep = Rand::genmultinomial( sumcounts, popfreqs );

      // Calculate likelihood of observed and repliate data.
      LogLikelihood[j][k] =
	log( MultinomialPDF( AlleleCount, popfreqs ) );
      RepLogLikelihood[j][k] =
	log( MultinomialPDF( rep, popfreqs ) );
      if(!( LogLikelihood[j][k] < RepLogLikelihood[j][k]) )
	divergentallelefreqstest[j][k] ++;
      sum[k] += LogLikelihood[j][k];
      repsum[k] += RepLogLikelihood[j][k];
    }
  }

  //test statistic for population
  for( int k = 0; k < options->getPopulations(); k++ ){
    if(!( sum[k] < repsum[k]) )
      divergentallelefreqstest[NumberOfCompositeLoci][k] ++;
  }
}

void DispersionTest::Output(int samples, const Genome& Loci, const Vector_s& PopLabels){
  //write header
  dispersionoutputstream << "Locus";
  for(int i = 0; i< NumberOfPopulations; ++i)
    dispersionoutputstream << "\t\"" << PopLabels[i]<<"\""; 
  dispersionoutputstream << endl;
  
  for(int j = 0; j < NumberOfCompositeLoci; j++ ){
    dispersionoutputstream << Loci(j)->GetLabel(0);
    for(int k=0; k < NumberOfPopulations; ++k)
      dispersionoutputstream << "\t" << (float)divergentallelefreqstest[j][k] / (float)samples;
    dispersionoutputstream  << endl;
  }
  dispersionoutputstream << "Population";
  for(int k=0; k < NumberOfPopulations; ++k)
    dispersionoutputstream << "\t" << (float)divergentallelefreqstest[NumberOfCompositeLoci][k] /(float) samples;
  dispersionoutputstream  << endl;
}
