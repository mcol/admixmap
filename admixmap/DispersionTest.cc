/** 
 *   ADMIXMAP
 *   DispersionTest.cc 
 *   Class to implement test for dispersion of allele frequencies between the unadmixed populations sampled 
 *   and the corresponding ancestry-specific allele frequencies in the admixed population under study.  
 *   This is evaluated for each subpopulation at each locus, and as a global test over all loci.  
 *   Valid only if option priorallelefreqfile is specified. The results are "Bayesian p-values". 
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "DispersionTest.h"

using namespace ::std;

DispersionTest::DispersionTest(){
  divergentallelefreqstest = 0;
  NumberOfCompositeLoci = -1;
  NumberOfPopulations = 0;
  options = 0;
}

void DispersionTest::Initialise(AdmixOptions *op,LogWriter *Log, int NumLoci){
  options = op;
  NumberOfCompositeLoci = NumLoci;
  NumberOfPopulations = options->getPopulations();

  divergentallelefreqstest = alloc2D_i(NumberOfCompositeLoci + 1, NumberOfPopulations );
  
  if( options->getTestForDispersion() ){
    Log->logmsg(true, "Writing dispersion test results to ");
    Log->logmsg(true,options->getDispersionTestFilename());
    Log->logmsg(true,"\n");
    dispersionoutputstream.open( options->getDispersionTestFilename(), ios::out );
    if( !dispersionoutputstream ){
      Log->logmsg(true,"ERROR: Couldn't open dispersiontestfile\n");
   exit( 1 );
    }
  }
}

DispersionTest::~DispersionTest(){
  free_matrix(divergentallelefreqstest, NumberOfCompositeLoci + 1);
}

void DispersionTest::TestForDivergentAlleleFrequencies(AlleleFreqs *A)
{
  int numberofstates;
  Vector_i rep;
  Vector_d popfreqs;

  Vector_i AlleleCount; 
  double LogLikelihood[NumberOfCompositeLoci + 1][ NumberOfPopulations ], 
    RepLogLikelihood[NumberOfCompositeLoci + 1][ NumberOfPopulations ];
  double sum[NumberOfPopulations], repsum[NumberOfPopulations];

  for( int k = 0; k < NumberOfPopulations; k++ ){
    sum[k] = 0.0; 
    repsum[k] = 0.0;
  }
  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    numberofstates = A->getLocus(j)->GetNumberOfStates();
    //AlleleCount = A->GetAlleleCounts(j);
    //freqs = A->GetAlleleFreqs(j);
    for( int k = 0; k < NumberOfPopulations; k++ ){
      popfreqs =  A->GetAlleleFreqs(j, k );
      popfreqs.AddElement( numberofstates - 1 );
      popfreqs( numberofstates - 1 ) = 1 - popfreqs.Sum();
      // Generate replicate data conditional on locus ancestry
      AlleleCount = A->GetAlleleCounts(j, k);
      rep = genmultinomial( AlleleCount.Sum(), popfreqs );

      // Calculate likelihood of observed and repliate data.
      LogLikelihood[j][k] =
	log( MultinomialLikelihood( AlleleCount, popfreqs ) );
      RepLogLikelihood[j][k] =
	log( MultinomialLikelihood( rep, popfreqs ) );
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

void DispersionTest::Output(int samples,Genome& Loci, std::string *PopLabels){
  //write header
  dispersionoutputstream << "Locus";
  for(int i = 0; i< NumberOfPopulations; ++i)
    dispersionoutputstream << " " <<PopLabels[i];//should use pop labels here
  dispersionoutputstream << endl;

  for(int j = 0; j < NumberOfCompositeLoci; j++ ){
    if(options->IsPedFile())
      dispersionoutputstream << "\"" << Loci(j)->GetLabel(0) << "\"" << " ";
    else
      dispersionoutputstream << Loci(j)->GetLabel(0) << " ";
    for(int k=0; k < NumberOfPopulations; ++k)
      dispersionoutputstream << (float)divergentallelefreqstest[j][k] / (float)samples << " ";
    dispersionoutputstream  << endl;
  }
  dispersionoutputstream << "Population ";
  for(int k=0; k < NumberOfPopulations; ++k)
    dispersionoutputstream << (float)divergentallelefreqstest[NumberOfCompositeLoci][k] /(float) samples << " ";
  dispersionoutputstream  << endl;
}
