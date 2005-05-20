/*
  Implements score tests for mis-specified allele frequencies
 * Only used with fixed allele frequencies and only for monitoring. 

  There are currently two tests:
  1. scalar test carried out only at composite loci with a single locus
  2. a generalized test , not fully implemented yet

*/
#include "MisSpecAlleleFreqTest.h"
#include "functions.h"
#include "VectorLoop.h"

MisSpecAlleleFreqTest::MisSpecAlleleFreqTest(){
  ScoreGene = 0;
  SumScoreGene = 0;
  InfoGene = 0;
  SumInfoGene = 0;
  SumScoreGeneSq = 0;
  SumNewScore = 0;
  SumNewInfo = 0;
  SumNewScoreSq = 0;
  NumTestLoci = 0;
  NumCompLoci = 0;
  Test1 = false;
  Test2 = false;
  Populations = 1;
}

MisSpecAlleleFreqTest::~MisSpecAlleleFreqTest(){
  delete[] ScoreGene;
  delete[] SumScoreGene;
  delete[] InfoGene;
  delete[] SumInfoGene;
  delete[] SumScoreGeneSq;
  for(int i = 0; i < NumCompLoci; ++i){
    delete[] SumNewScore[i];
    delete[] SumNewInfo[i];
    delete[] SumNewScoreSq[i];
  }
  delete[] SumNewScore;
  delete[] SumNewInfo;
  delete[] SumNewScoreSq;
}


void MisSpecAlleleFreqTest::Initialise(AdmixOptions *options, Genome *Loci, LogWriter *Log )
{
 if( strlen( options->getAlleleFreqFilename() ) || options->getFixedAlleleFreqs()){

   Test1 = options->getTestForMisspecifiedAlleleFreqs();
   Test2 = options->getTestForMisspecifiedAlleleFreqs2();
   Populations = options->getPopulations();
   NumCompLoci = Loci->GetNumberOfCompositeLoci();
   
   if(Test1){
     NumTestLoci = 0;
     for(int i = 0; i < NumCompLoci; ++i)
       if((*Loci)(i)->GetNumberOfLoci() == 1) NumTestLoci++;

     ScoreGene = new Matrix_d[NumCompLoci];
     SumScoreGene = new Matrix_d[NumCompLoci];
     InfoGene = new Matrix_d[NumCompLoci];
     SumInfoGene = new Matrix_d[NumCompLoci];
     SumScoreGeneSq = new Matrix_d[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       ScoreGene[i].SetNumberOfElements( Populations, 1 );
       SumScoreGene[i].SetNumberOfElements( Populations, 1 );
       InfoGene[i].SetNumberOfElements( Populations, Populations );
       SumInfoGene[i].SetNumberOfElements( Populations, Populations );
       SumScoreGeneSq[i].SetNumberOfElements( Populations, Populations );
    }
     allelefreqscorestream.open( options->getAlleleFreqScoreFilename() );
     if( !allelefreqscorestream ){
      Log->logmsg(true,"ERROR: Couldn't open allelefreqscorefile\n");
      Log->logmsg(true,options->getAlleleFreqScoreFilename());
      Log->logmsg(true,"\n");
      exit( 1 );
     }
    else{
      Log->logmsg(true,"Writing score tests for mis-specified allele frequencies to ");
      Log->logmsg(true,options->getAlleleFreqScoreFilename());
      Log->logmsg(true,"\n");
      allelefreqscorestream << "structure(.Data=c(" << endl;
    }
   }
   
   if(Test2){
     SumNewScore = new Matrix_d*[NumCompLoci];
     SumNewInfo = new Matrix_d*[NumCompLoci];
     SumNewScoreSq = new Matrix_d*[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       SumNewScore[i] = new Matrix_d[Populations];
       SumNewInfo[i] = new Matrix_d[Populations];
       SumNewScoreSq[i] = new Matrix_d[Populations];

       int NumberOfStates = (*Loci)(i)->GetNumberOfStates();
       for(int k = 0; k < Populations; ++k){
	 SumNewScore[i][k].SetNumberOfElements(NumberOfStates - 1, 1 );
	 SumNewInfo[i][k].SetNumberOfElements(NumberOfStates - 1, NumberOfStates - 1 );
	 SumNewScoreSq[i][k].SetNumberOfElements(NumberOfStates - 1, NumberOfStates - 1 );
       }
     }
     allelefreqscorestream2.open( options->getAlleleFreqScoreFilename2() );
     if( !allelefreqscorestream2 ){
       Log->logmsg(true,"ERROR: Couldn't open allelefreqscorefile\n");
       Log->logmsg(true,options->getAlleleFreqScoreFilename2());
       Log->logmsg(true,"\n");
       exit( 1 );
    }
     else{
       Log->logmsg(true,"Writing score tests for mis-specified allele frequencies to ");
       Log->logmsg(true,options->getAlleleFreqScoreFilename2());
       Log->logmsg(true,"\n");
       allelefreqscorestream2 << "structure(.Data=c(" << endl;
     }
   }
 }
 else{
   Log->logmsg(true, "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n");
   Log->logmsg(true, "This option will be ignored\n");
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecAlleleFreqTest::Update(IndividualCollection *individuals, AlleleFreqs *A, Genome *Loci){
  if(Test1)
    for(int j = 0; j < NumCompLoci; ++j){
      ScoreGene[j].SetElements(0);
      InfoGene[j].SetElements(0);
    }

  if( Test1 ) {
    dmatrix phi = alloc2D_d(Populations, Populations);
    for( int i = 0; i < individuals->getSize(); i++ ){
      Individual* ind = individuals->getIndividual(i);
      
      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  phi[k][ kk ] = ind->getAdmixtureProps()( k, 0 ) * ind->getAdmixtureProps()( kk, 0 );
      
      for(int j = 0; j < NumCompLoci; j++ ){
	if( !(ind->IsMissing(j)) && 
	    (*Loci)(j)->GetNumberOfLoci() == 1  && !(A->IsRandom()) ){
	  UpdateScoreForMisSpecOfAlleleFreqs( j, phi, ind->getGenotype(j), A->GetAlleleFreqs(j) );
	}
      }
    }
    for(int j = 0; j < NumCompLoci; j++ ){
      SumScoreGene[j] += ScoreGene[j];
      SumInfoGene[j] += InfoGene[j];
      SumScoreGeneSq[j] += ScoreGene[j] * ScoreGene[j].Transpose();
    }
    free_matrix(phi, Populations);
  }

  if( Test2 )    
    for( int j = 0; j < NumCompLoci; j++ ){
      for( int i = 0; i < individuals->getSize(); i++ ){
	UpdateScoreForMisSpecOfAlleleFreqs2(j, (*Loci)(j)->GetNumberOfStates(), A->GetAlleleFreqs(j), A->GetAlleleCounts(j));
      }
    }

}

/**
 * N.B. This only works for a single SNiP.
 * Updates what's required for the score tests. Only used with fixed
 * allele frequencies. This function is only used for monitoring.
 */
 void MisSpecAlleleFreqTest::UpdateScoreForMisSpecOfAlleleFreqs(int j,  dmatrix phi, unsigned short **x, Matrix_d AlleleFreqs)
{
   double Score[ Populations ];
   double Pi[3] = {0.0, 0.0, 0.0};

   for( int k = 0; k < Populations; k++ ){
      for( int kk = 0; kk < Populations; kk++ ){
         Pi[0] += AlleleFreqs( 0, k ) * AlleleFreqs( 0, kk ) * phi[k][kk];
         Pi[1] += ( ( 1 - AlleleFreqs( 0, k ) ) * AlleleFreqs( 0, kk ) + ( 1 - AlleleFreqs( 0, kk ) ) * AlleleFreqs( 0, k ) ) *phi[k][kk];
         Pi[2] += ( 1 - AlleleFreqs( 0, k ) ) * ( 1 - AlleleFreqs( 0, kk ) ) * phi[k][kk];
      }
   }

   if( x[0][0] == 1 && x[0][1] == 1 ){
      for( int k = 0; k < Populations; k++ ){
         Score[k] = 2 * AlleleFreqs( 0, k ) * phi[k][k];
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
               Score[k] += AlleleFreqs( 0, kk ) * (phi[k][kk] + phi[kk][k]);
         Score[k] /= Pi[0];
         ScoreGene[j]( k, 0 ) += Score[k];
         InfoGene[j]( k, k ) += Score[k] * Score[k] - 2 * phi[k][k] / Pi[0];
      }
      for( int k = 0; k < Populations; k++ )
         for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                InfoGene[j]( k, kk ) += Score[k] * Score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[0];}
   
   else if( x[0][0] == 1 && x[0][1] != 1 ){
      for( int k = 0; k < Populations; k++ ){
         Score[k] = 2 * ( 1 - 2 * AlleleFreqs( 0, k ) ) * phi[k][k];
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
                Score[k] += ( 1 - 2 * AlleleFreqs( 0, kk ) ) * (phi[k][kk] + phi[kk][k]);
         Score[k] /= Pi[1];
         ScoreGene[j]( k, 0 ) += Score[k];
         InfoGene[j]( k, k ) += Score[k] * Score[k] + 4 * phi[k][k] / Pi[1];}
       for( int k = 0; k < Populations; k++ )
          for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                InfoGene[j]( k, kk ) += Score[k] * Score[kk] + 2*(phi[k][kk] + phi[kk][k]) / Pi[1];}
   
   else if( x[0][0] != 0 && x[0][0] != 1 && x[0][1] != 1 ){
      for( int k = 0; k < Populations; k++ ){
          Score[k] = -2 * ( 1 - AlleleFreqs( 0, k ) ) * phi[k][k];
          for( int kk = 0; kk < Populations; kk++ )
             if( k != kk )
                Score[k] -= ( 1 - AlleleFreqs( 0, kk ) ) * (phi[k][kk] + phi[kk][k]);
          Score[k] /= Pi[2];
          ScoreGene[j]( k, 0 ) += Score[k];
          InfoGene[j]( k, k ) += Score[k] * Score[k] - 2 * phi[k][k] / Pi[2];}
      for( int k = 0; k < Populations; k++ )
         for( int kk = 0; kk < Populations; kk++ )
            if( k != kk )
               InfoGene[j]( k, kk ) += Score[k] * Score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[2];}
}

// presumably this calculates score test for mis-spec allele freqs at multi-allelic loci
void MisSpecAlleleFreqTest::UpdateScoreForMisSpecOfAlleleFreqs2(const int locus, const int NumberOfStates, const Matrix_d &AlleleFreqs, 
								const Matrix_i &AlleleCounts)
{
   double rn, r, pj, pi, q;
   Matrix_d NewScore( NumberOfStates - 1, 1 ), NewInfo( NumberOfStates - 1, NumberOfStates - 1 );
   for( int k = 0; k < Populations; k++ ){
      rn = (double)( AlleleCounts( NumberOfStates - 1, k ) );
      q = 1 - AlleleFreqs.GetColumn(k).Sum();
      for( int j = 0; j < NumberOfStates - 1; j++ ){
         r = AlleleCounts( j, k );
         pj = AlleleFreqs( j, k );
         NewScore( j, 0 ) = ( r / pj - rn / q ) * pj * ( 1 - pj );
         NewInfo( j, j ) = pj * ( 1 - pj )
            * ( r - ( rn / q ) * ( 2*pj - 1.0 - pj / q + pj * pj / q ) );
         for( int i = j+1; i < NumberOfStates - 1; i++ ){
            pi = AlleleFreqs( i, k );
            NewInfo( j, i ) = rn * pj * ( 1 - pj ) * pi * ( 1 - pi ) / ( q * q );
            NewInfo( i, j ) = NewInfo( j, i );
         }
      }
      SumNewScore[locus][k] += NewScore;
      SumNewInfo[locus][k] += NewInfo;
      SumNewScoreSq[locus][k] += NewScore * NewScore.Transpose();
   }
}


void MisSpecAlleleFreqTest::Output(int samples, Genome *Loci,  std::string * PopLabels, bool IsPedFile){
    if( Test1){
      OutputTestsForMisSpecifiedAlleleFreqs(samples, Loci, PopLabels, IsPedFile);
    }
    if( Test2 ){
      OutputTestsForMisSpecifiedAlleleFreqs2(samples, Loci, PopLabels, IsPedFile);
    }
}

void MisSpecAlleleFreqTest::OutputTestsForMisSpecifiedAlleleFreqs( int samples, Genome *Loci, std::string * PopLabels,bool IsPedFile)
{
  //int samples = iteration - options->getBurnIn();
  Matrix_d ScoreMatrix, CompleteMatrix, ObservedMatrix;
  for(int j = 0; j < NumCompLoci; j++ ){
    if( (*Loci)(j)->GetNumberOfLoci() == 1 ){
      ScoreMatrix = SumScoreGene[j] / samples;
      CompleteMatrix = SumInfoGene[j] / samples;
      ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumScoreGeneSq[j] / samples;
      for( int k = 0; k < Populations; k++ ){
	// Test for mis-specification within each continental-population.
	if(IsPedFile)
	  allelefreqscorestream << "\"" << (*Loci)(j)->GetLabel(0) << "\"" << ",";
	else
	  allelefreqscorestream << (*Loci)(j)->GetLabel(0) << ",";
	allelefreqscorestream << PopLabels[k] << ",";
	allelefreqscorestream << double2R(ScoreMatrix( k, 0 ) ) << ",";
	allelefreqscorestream << double2R(CompleteMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(ObservedMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(100*ObservedMatrix( k, k ) / CompleteMatrix( k, k ) ) << ",";
	allelefreqscorestream << double2R(ScoreMatrix( k, 0 ) / sqrt( ObservedMatrix( k, k ) ) ) << ",";
	if( k < Populations - 1 )
	  allelefreqscorestream  << "\"NA\"," << endl;
      }
      // Summary chi-sq test for mis-specification in all continental-populations.
      ObservedMatrix = ObservedMatrix.SubMatrix( 0, Populations - 1, 0, Populations - 1 );
      ScoreMatrix = ScoreMatrix.SubMatrix( 0, Populations - 1, 0, 0 );
      ObservedMatrix.InvertUsingLUDecomposition();
      allelefreqscorestream << (ScoreMatrix.Transpose() * ObservedMatrix * ScoreMatrix)(0,0) << "," << endl;
    }
  }
  /**
   * writes out the dimensions and labels of the 
   * R- array for old test for mis-specified allele frequencies
   */
  vector<int> dimensions(3,0);
  dimensions[0] = 8;
  dimensions[1] = NumTestLoci * Populations;
  dimensions[2] = 1;//(int)(numPrintedIterations);
  
  vector<string> labels(8,"");
  labels[0] = "Locus";
  labels[1] = "Population";
  labels[2] = "Score";
  labels[3] = "CompleteInfo";
  labels[4] = "ObservedInfo";
  labels[5] = "PercentInfo";
  labels[6] = "StdNormal";
  labels[7] = "ChiSquared";
  
  R_output3DarrayDimensions(&allelefreqscorestream,dimensions,labels);
}

void MisSpecAlleleFreqTest::OutputTestsForMisSpecifiedAlleleFreqs2( int samples, Genome *Loci, std::string * PopLabels,bool IsPedFile)
{
  Matrix_d score, completeinfo, observedinfo;
  for(int j = 0; j < NumCompLoci; j++ ){
    //CompositeLocus *locus = (CompositeLocus*)(*Lociptr)(j);
    for( int k = 0; k < Populations; k++ ){
      score = SumNewScore[j][k] / samples;
      completeinfo = SumNewInfo[j][k] / samples;
      observedinfo = completeinfo + score * score.Transpose() - SumNewScoreSq[j][k] / samples;
      if(IsPedFile)
	allelefreqscorestream2 << "\"" << (*Loci)(j)->GetLabel(0) << "\"" << ",";
      else
	allelefreqscorestream2 << (*Loci)(j)->GetLabel(0) << ",";
      allelefreqscorestream2 << PopLabels[k] << ",";
      allelefreqscorestream2 << double2R(completeinfo.Determinant()) << ",";
      allelefreqscorestream2 << double2R(observedinfo.Determinant()) << ",";
      allelefreqscorestream2 << double2R(100*observedinfo.Determinant() / completeinfo.Determinant()) << ",";
      observedinfo.InvertUsingLUDecomposition();
      allelefreqscorestream2 << double2R((score.Transpose() * observedinfo * score)(0,0)) << "," << endl;
    }
  }
  /**
   * writes out the dimensions and labels of the 
   * R-array for new test for mis-specified allele frequencies
   */
    vector<int> dimensions(3,0);
    dimensions[0] = 6;
    dimensions[1] = NumCompLoci * Populations;
    dimensions[2] = 1;//(int)(numPrintedIterations);

    vector<string> labels(6,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "ChiSquared";

    R_output3DarrayDimensions(&allelefreqscorestream2,dimensions,labels);

}

void MisSpecAlleleFreqTest::R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels)
{
  *stream << ")," << endl;
  *stream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    *stream << dim[i];
    if(i != dim.size() - 1){
      *stream << ",";
    }
  }
  *stream << ")," << endl;
  *stream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    *stream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      *stream << ",";
    }
  }
  *stream << "), character(0), character(0)))" << endl;


}
string MisSpecAlleleFreqTest::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
