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
}
void DispersionTest::Initialise(AdmixOptions *op,LogWriter *Log, int NumberOfCompositeLoci){
  options = op;

  divergentallelefreqstest = Matrix_i(NumberOfCompositeLoci + 1, options->getPopulations() );
  
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

void DispersionTest::TestForDivergentAlleleFrequencies(AlleleFreqs *A)
{
  int numberofstates;
  Vector_i locusancestry, rep, PopCounts;
  Vector_d popfreqs;

  Matrix_i AlleleCount, test( A->GetNumberOfCompositeLoci() + 1, options->getPopulations() );
  Matrix_d LogLikelihood( A->GetNumberOfCompositeLoci() + 1, options->getPopulations() ), 
    RepLogLikelihood( A->GetNumberOfCompositeLoci() + 1, options->getPopulations() ),
    freqs;

  for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
    numberofstates = A->getLocus(j)->GetNumberOfStates();
    AlleleCount = A->GetAlleleCounts(j);
    freqs = A->GetAlleleFreqs(j);
    for( int k = 0; k < options->getPopulations(); k++ ){
      popfreqs = freqs.GetColumn( k );
      popfreqs.AddElement( numberofstates - 1 );
      popfreqs( numberofstates - 1 ) = 1 - popfreqs.Sum();
      // Generate replicate data conditional on locus ancestry
      rep = genmultinomial( (AlleleCount.GetColumn(k)).Sum(), popfreqs );

      // Calculate likelihood of observed and repliate data.
      LogLikelihood( j, k ) =
	log( MultinomialLikelihood( AlleleCount.GetColumn( k ), popfreqs ) );
      RepLogLikelihood( j, k ) =
	log( MultinomialLikelihood( rep, popfreqs ) );
      if( LogLikelihood( j, k ) < RepLogLikelihood( j, k ) )
	test( j, k ) = 0;
      else
	test( j, k ) = 1;
    }
  }

  for( int k = 0; k < options->getPopulations(); k++ ){
    if( LogLikelihood.GetColumn(k).Sum() < RepLogLikelihood.GetColumn(k).Sum() )
      test( A->GetNumberOfCompositeLoci(), k ) = 0;
    else
      test( A->GetNumberOfCompositeLoci(), k ) = 1;
  }

  divergentallelefreqstest += test;
}

void DispersionTest::Output(int samples,Genome& Loci, std::string *PopLabels){
  //write header
  dispersionoutputstream << "Locus";
  for(int i = 0; i< divergentallelefreqstest.GetNumberOfCols(); ++i)
    dispersionoutputstream << " " <<PopLabels[i];//should use pop labels here
  dispersionoutputstream << endl;

  for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
    if(options->IsPedFile())
      dispersionoutputstream << "\"" << Loci(j)->GetLabel(0) << "\"" << " ";
    else
      dispersionoutputstream << Loci(j)->GetLabel(0) << " ";
    dispersionoutputstream << divergentallelefreqstest.GetRow(j).Float() / samples << endl;
  }
  dispersionoutputstream << "Population "
			 << divergentallelefreqstest.GetRow( Loci.GetNumberOfCompositeLoci() ).Float() / samples
			 << endl;
}
