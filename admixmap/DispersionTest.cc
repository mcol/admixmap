#include "DispersionTest.h"

using namespace ::std;

DispersionTest::DispersionTest(){
}
void DispersionTest::Initialise(AdmixOptions *op,LogWriter *Log, int NumberOfCompositeLoci){
  options = op;

  divergentallelefreqstest = Matrix_i(NumberOfCompositeLoci + 1, options->getPopulations() );
  
  if( options->getTestForDispersion() ){
    Log->logmsg(true,options->getDispersionTestFilename());
    Log->logmsg(true,"\n");
    dispersionoutputstream.open( options->getDispersionTestFilename(), ios::out );
    if( !dispersionoutputstream ){
      Log->logmsg(true,"ERROR: Couldn't open dispersiontestfile\n");
   exit( 1 );
    }
  }
}
void
DispersionTest::UpdateBayesianPValueTest(AlleleFreqs *A)
{
  if( options->getTestForDispersion() ){
    divergentallelefreqstest += TestForDivergentAlleleFrequencies(A);
  }
  if( options->getTestForMisspecifiedAlleleFreqs2() ){
    for( int j = 0; j < A->GetNumberOfCompositeLoci(); j++ ){
      CompositeLocus *locus = (CompositeLocus*)A->getLocus(j);
      locus->UpdateScoreForMisSpecOfAlleleFreqs2(A->GetAlleleFreqs(j), A->GetAlleleCounts(j));
    }
  }

}

Matrix_i DispersionTest::TestForDivergentAlleleFrequencies(AlleleFreqs *A)
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

  return( test );
}

void DispersionTest::Output(int samples,Genome& Loci){
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
