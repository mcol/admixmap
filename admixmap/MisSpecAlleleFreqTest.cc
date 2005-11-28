/*
  Implements score tests for mis-specified allele frequencies
 * Only used with fixed allele frequencies and only for monitoring. 

  There are currently two tests:
  1. scalar test carried out only at composite loci with a single locus
  2. a generalized test , not fully implemented yet

*/
#include "MisSpecAlleleFreqTest.h"
#include "functions.h"

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
  if(Test1){
    for(int i = NumCompLoci; i > 0; --i){
      delete[] ScoreGene[i-1];
      delete[] InfoGene[i-1];
      delete[] SumScoreGene[i-1];
      delete[] SumInfoGene[i-1];
      delete[] SumScoreGeneSq[i-1];  
    }
    delete[] ScoreGene;
    delete[] SumScoreGene;
    delete[] InfoGene;
    delete[] SumInfoGene;
    delete[] SumScoreGeneSq;
  }
  if(Test2){
    for(int i = NumCompLoci; i > 0; --i){
      for(int k = Populations; k >0; --k){
	delete[] SumNewScore[i-1][k-1];
 	delete[] SumNewInfo[i-1][k-1];
	delete[] SumNewScoreSq[i-1][k-1];
      }

      delete[] SumNewScore[i-1];
      delete[] SumNewInfo[i-1];
      delete[] SumNewScoreSq[i-1];
    }
    delete[] SumNewScore;
    delete[] SumNewInfo;
    delete[] SumNewScoreSq;
  }
}


void MisSpecAlleleFreqTest::Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log )
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

     ScoreGene = new double*[NumCompLoci];
     SumScoreGene = new double*[NumCompLoci];
     InfoGene = new double*[NumCompLoci];
     SumInfoGene = new double*[NumCompLoci];
     SumScoreGeneSq = new double*[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       ScoreGene[i] = new double[ Populations ];
       InfoGene[i] = new double[ Populations * Populations ]; 
       SumScoreGene[i] = new double[ Populations ];
       fill(SumScoreGene[i], SumScoreGene[i]+Populations, 0.0);
       SumInfoGene[i] = new double[ Populations * Populations ];
       fill(SumInfoGene[i], SumInfoGene[i]+ Populations*Populations, 0.0);
       SumScoreGeneSq[i] = new double[ Populations * Populations ];
       fill(SumScoreGeneSq[i], SumScoreGeneSq[i]+Populations*Populations, 0.0);
    }
     Log.setDisplayMode(On);
     allelefreqscorestream.open( options->getAlleleFreqScoreFilename() );
     if( !allelefreqscorestream ){
       Log << "ERROR: Couldn't open allelefreqscorefile " << options->getAlleleFreqScoreFilename() << "\n";
      exit( 1 );
     }
    else{
      Log.setDisplayMode(Quiet);
      Log << "Writing score tests for mis-specified allele frequencies(1) to " << options->getAlleleFreqScoreFilename() << "\n";
      allelefreqscorestream << "structure(.Data=c(" << endl;
    }
   }
   
   if(Test2){
     SumNewScore = new double**[NumCompLoci];
     SumNewInfo = new double**[NumCompLoci];
     SumNewScoreSq = new double**[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       SumNewScore[i] = new double*[Populations];
       SumNewInfo[i] = new double*[Populations];
       SumNewScoreSq[i] = new double*[Populations];

       int NumberOfStates = (*Loci)(i)->GetNumberOfStates();
       for(int k = 0; k < Populations; ++k){
	 SumNewScore[i][k] = new double[NumberOfStates - 1];
	 fill(SumNewScore[i][k], SumNewScore[i][k] + (NumberOfStates-1), 0.0);
	 SumNewInfo[i][k] = new double[ (NumberOfStates - 1)*( NumberOfStates - 1 )];
	 fill(SumNewInfo[i][k], SumNewInfo[i][k] + (NumberOfStates-1)*(NumberOfStates-1), 0.0);
	 SumNewScoreSq[i][k] = new double[(NumberOfStates - 1) * (NumberOfStates - 1 )];
	 fill(SumNewScoreSq[i][k], SumNewScoreSq[i][k] + (NumberOfStates-1)*(NumberOfStates-1), 0.0);
       }
     }
     allelefreqscorestream2.open( options->getAlleleFreqScoreFilename2() );
     if( !allelefreqscorestream2 ){
       Log << "ERROR: Couldn't open allelefreqscorefile " << options->getAlleleFreqScoreFilename2() << "\n";
       exit( 1 );
    }
     else{
       Log.setDisplayMode(Quiet);
       Log << "Writing score tests for mis-specified allele frequencies(2) to " << options->getAlleleFreqScoreFilename2() << "\n";
       allelefreqscorestream2 << "structure(.Data=c(" << endl;
     }
   }
 }
 else{
   Log << "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n"
       << "This option will be ignored\n";
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecAlleleFreqTest::Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci){
  if(Test1)
    for(int j = 0; j < NumCompLoci; ++j){
      fill(ScoreGene[j], ScoreGene[j]+Populations, 0.0);
      fill(InfoGene[j], InfoGene[j]+Populations*Populations, 0.0);
    }

  if( Test1 ) {
    double** phi = alloc2D_d(Populations, Populations);
    for( int i = 0; i < individuals->getSize(); i++ ){
      const Individual* ind = individuals->getIndividual(i);
      
      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  phi[k][ kk ] = ind->getAdmixtureProps()[k] * ind->getAdmixtureProps()[kk];
      
      for(int j = 0; j < NumCompLoci; j++ ){
	if( !(ind->IsMissing(j)) && 
	    (*Loci)(j)->GetNumberOfLoci() == 1  && !(A->IsRandom()) ){
	  UpdateScoreForMisSpecOfAlleleFreqs( j, phi, ind->getGenotype(j), A->GetAlleleFreqs(j) );
	}
      }
    }
    for(int j = 0; j < NumCompLoci; j++ ){
      //SumScoreGene[j] += ScoreGene[j];
      transform(ScoreGene[j], ScoreGene[j]+Populations, SumScoreGene[j], SumScoreGene[j], std::plus<double>());
      //SumInfoGene[j] += InfoGene[j];
      transform(InfoGene[j], InfoGene[j]+Populations*Populations, SumInfoGene[j], SumInfoGene[j], std::plus<double>());
      //SumScoreGeneSq[j] += ScoreGene[j] * ScoreGene[j].Transpose();
      double* ScoreGeneSq = new double[Populations*Populations];
      matrix_product(ScoreGene[j], ScoreGeneSq, Populations, 1);
      transform(ScoreGeneSq, ScoreGeneSq+Populations*Populations, SumScoreGeneSq[j], SumScoreGeneSq[j], std::plus<double>());
      delete[] ScoreGeneSq;
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
 void MisSpecAlleleFreqTest::UpdateScoreForMisSpecOfAlleleFreqs(int j, const double* const* phi, const vector<vector<unsigned short> > x, 
								const double* const AlleleFreqs)
{
  vector<double> Score( Populations );
  double Pi[3] = {0.0, 0.0, 0.0};
  
  for( int k = 0; k < Populations; k++ ){
    for( int kk = 0; kk < Populations; kk++ ){
      Pi[0] += AlleleFreqs[ k ] * AlleleFreqs[ kk ] * phi[k][kk];
      Pi[1] += ( ( 1 - AlleleFreqs[ k ] ) * AlleleFreqs[kk ] + ( 1 - AlleleFreqs[ kk ] ) * AlleleFreqs[ k ] ) *phi[k][kk];
      Pi[2] += ( 1 - AlleleFreqs[ k ] ) * ( 1 - AlleleFreqs[ kk ] ) * phi[k][kk];
    }
  }
  
  if( x[0][0] == 1 && x[0][1] == 1 ){
    for( int k = 0; k < Populations; k++ ){
      Score[k] = 2 * AlleleFreqs[ k ] * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  Score[k] += AlleleFreqs[ kk ] * (phi[k][kk] + phi[kk][k]);
      Score[k] /= Pi[0];
      ScoreGene[j][ k ] += Score[k];
      InfoGene[j][ k*Populations + k ] += Score[k] * Score[k] - 2 * phi[k][k] / Pi[0];
    }
    for( int k = 0; k < Populations; k++ )
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  InfoGene[j][ k*Populations + kk ] += Score[k] * Score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[0];}
  
  else if( x[0][0] == 1 && x[0][1] != 1 ){
    for( int k = 0; k < Populations; k++ ){
      Score[k] = 2 * ( 1 - 2 * AlleleFreqs[ k ] ) * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  Score[k] += ( 1 - 2 * AlleleFreqs[  kk ] ) * (phi[k][kk] + phi[kk][k]);
         Score[k] /= Pi[1];
         ScoreGene[j][ k ] += Score[k];
         InfoGene[j][ k*Populations + k ] += Score[k] * Score[k] + 4 * phi[k][k] / Pi[1];}
       for( int k = 0; k < Populations; k++ )
	 for( int kk = 0; kk < Populations; kk++ )
	   if( k != kk )
	     InfoGene[j][ k*Populations + kk ] += Score[k] * Score[kk] + 2*(phi[k][kk] + phi[kk][k]) / Pi[1];}
  
  else if( x[0][0] != 0 && x[0][0] != 1 && x[0][1] != 1 ){
    for( int k = 0; k < Populations; k++ ){
      Score[k] = -2 * ( 1 - AlleleFreqs[ k ] ) * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  Score[k] -= ( 1 - AlleleFreqs[ kk ] ) * (phi[k][kk] + phi[kk][k]);
      Score[k] /= Pi[2];
      ScoreGene[j][ k ] += Score[k];
      InfoGene[j][ k*Populations + k ] += Score[k] * Score[k] - 2 * phi[k][k] / Pi[2];}
      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  if( k != kk )
	    InfoGene[j][ k*Populations + kk ] += Score[k] * Score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[2];}
}

// presumably this calculates score test for mis-spec allele freqs at multi-allelic loci
void MisSpecAlleleFreqTest::UpdateScoreForMisSpecOfAlleleFreqs2(const int locus, const int NumberOfStates, 
								const double* const AlleleFreqs, const int* const AlleleCounts)
{
   double rn, r, pj, pi, q;
   double* NewScore = new double[ NumberOfStates - 1];
   double* NewInfo = new double[ (NumberOfStates - 1) * (NumberOfStates - 1 )];
   for( int k = 0; k < Populations; k++ ){
     rn = (double)( AlleleCounts[ (NumberOfStates - 1)*Populations + k ] );
     q = 1.0;
     for( int j = 0; j < NumberOfStates - 1; j++ )
       q -= AlleleFreqs[j*Populations+k];
      for( int j = 0; j < NumberOfStates - 1; j++ ){
         r = AlleleCounts[ j*Populations + k ];
         pj = AlleleFreqs[ j*Populations+ k ];
         NewScore[ j ] = ( r / pj - rn / q ) * pj * ( 1 - pj );
         NewInfo[ j*(NumberOfStates-1) + j ] = pj * ( 1 - pj )
            * ( r - ( rn / q ) * ( 2*pj - 1.0 - pj / q + pj * pj / q ) );
         for( int i = j+1; i < NumberOfStates - 1; i++ ){
            pi = AlleleFreqs[ i*Populations+ k ];
            NewInfo[ j*(NumberOfStates-1) + i ] = rn * pj * ( 1 - pj ) * pi * ( 1 - pi ) / ( q * q );
            NewInfo[ i*(NumberOfStates-1) + j ] = NewInfo[ j*(NumberOfStates-1) + i ];
         }
      }
      //SumNewScore[locus][k] += NewScore;
      transform(NewScore, NewScore+(NumberOfStates-1), SumNewScore[locus][k], SumNewScore[locus][k], std::plus<double>());
      //SumNewInfo[locus][k] += NewInfo;
      transform(NewInfo, NewInfo+(NumberOfStates-1)*(NumberOfStates-1), SumNewInfo[locus][k], SumNewInfo[locus][k], std::plus<double>());
      //SumNewScoreSq[locus][k] += NewScore * NewScore.Transpose();
      double* NewScoreSq = new double[(NumberOfStates - 1) * (NumberOfStates - 1 )];
      matrix_product(NewScore, NewScoreSq, NumberOfStates-1, 1);
      transform(NewScoreSq, NewScoreSq+(NumberOfStates-1)*(NumberOfStates-1), SumNewScoreSq[locus][k], SumNewScoreSq[locus][k], std::plus<double>());
      delete[] NewScoreSq;
   }
   delete[] NewScore;
   delete[] NewInfo;
}


void MisSpecAlleleFreqTest::Output(int samples, const Genome* const Loci, const std::string* const PopLabels){
    if( Test1){
      OutputTestsForMisSpecifiedAlleleFreqs(samples, Loci, PopLabels);
    }
    if( Test2 ){
      OutputTestsForMisSpecifiedAlleleFreqs2(samples, Loci, PopLabels);
    }
}

void MisSpecAlleleFreqTest::OutputTestsForMisSpecifiedAlleleFreqs( int samples, const Genome* const Loci, 
								   const std::string* const PopLabels)
{
  double* ObservedMatrix = new double[Populations*Populations];
  double* ScoreSq = new double[Populations*Populations];
  for(int j = 0; j < NumCompLoci; j++ ){
    if( (*Loci)(j)->GetNumberOfLoci() == 1 ){
 
      for(int k = 0; k < Populations; ++k)SumScoreGene[j][k] /= (double) samples;
      matrix_product(SumScoreGene[j], ScoreSq, Populations, 1);//ScoreSq = Score * t(Score)
      //CompleteMatrix = SumInfoGene[j] / samples;
      for(int k = 0; k < Populations*Populations; ++k){
	SumInfoGene[j][k] /= (double) samples;//SumInfoGene[j] = CompleteMatrix
	ObservedMatrix[k] = SumInfoGene[j][k] + ScoreSq[k] - SumScoreGeneSq[j][k] / (double)samples;
      }

      //ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumScoreGeneSq[j] / samples;



      for( int k = 0; k < Populations; k++ ){
	// Test for mis-specification within each continental-population.
	allelefreqscorestream << (*Loci)(j)->GetLabel(0) << ",";
	allelefreqscorestream << "\""<<PopLabels[k] << "\",";
	allelefreqscorestream << double2R(SumScoreGene[j][ k ] ) << ",";
	allelefreqscorestream << double2R(SumInfoGene[j][ k*Populations + k ] ) << ",";
	allelefreqscorestream << double2R(ObservedMatrix[ k*Populations + k ] ) << ",";
	allelefreqscorestream << double2R(100*ObservedMatrix[ k*Populations + k ] / SumInfoGene[j][ k*Populations + k ] ) << ",";
	allelefreqscorestream << double2R(SumScoreGene[j][ k ] / sqrt( ObservedMatrix[ k*Populations + k ] ) ) << ",";
	if( k < Populations - 1 )
	  allelefreqscorestream  << "\"NA\"," << endl;
      }
      // Summary chi-sq test for mis-specification in all continental-populations.
      double chisq = GaussianConditionalQuadraticForm(Populations, SumScoreGene[j], ObservedMatrix, Populations);
      allelefreqscorestream << double2R(chisq)<< "," << endl;
    }
  }
  delete[] ObservedMatrix;
  delete[] ScoreSq; 
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

void MisSpecAlleleFreqTest::OutputTestsForMisSpecifiedAlleleFreqs2( int samples, const Genome* const Loci, 
								    const std::string* const PopLabels)
{
  for(int j = 0; j < NumCompLoci; j++ ){
    int NumberOfStates = (*Loci)(j)->GetNumberOfStates();
    double* score = new double[NumberOfStates-1];
    double* completeinfo = new double[(NumberOfStates-1)*(NumberOfStates-1)];
    double* observedinfo = new double[(NumberOfStates-1)*(NumberOfStates-1)];
    double* scoresq = new double[(NumberOfStates-1)*(NumberOfStates-1)];
    
    //CompositeLocus *locus = (CompositeLocus*)(*Lociptr)(j);
    for( int k = 0; k < Populations; k++ ){
      for(int s = 0; s < NumberOfStates-1; ++s){
	score[s] = SumNewScore[j][k][s] / (double)samples;
      }

      matrix_product(score, scoresq, NumberOfStates-1, 1);//scoresq = score * t(score)
      for(int s = 0; s < (NumberOfStates-1)*(NumberOfStates-1); ++s){
	completeinfo[s] = SumNewInfo[j][k][s] / (double)samples;
	observedinfo[s] = completeinfo[s] + scoresq[s] - SumNewScoreSq[j][k][s] / (double)samples;
      }

      //score = SumNewScore[j][k] / samples;
      //completeinfo = SumNewInfo[j][k] / samples;
      //observedinfo = completeinfo + score * score.Transpose() - SumNewScoreSq[j][k] / samples;
      double det1 = determinant(completeinfo, NumberOfStates-1);
      double det2 = determinant(observedinfo, NumberOfStates-1);

      allelefreqscorestream2 << (*Loci)(j)->GetLabel(0) << ",";
      allelefreqscorestream2 << "\""<<PopLabels[k] << "\",";
      allelefreqscorestream2 << double2R(det1) << ",";
      allelefreqscorestream2 << double2R(det2) << ",";
      allelefreqscorestream2 << double2R(100*det2 / det1) << ",";
      double chisq = GaussianConditionalQuadraticForm(NumberOfStates-1, score, observedinfo, NumberOfStates-1);      
      //observedinfo.InvertUsingLUDecomposition();
      allelefreqscorestream2 << double2R(chisq)/*double2R((score.Transpose() * observedinfo * score)(0,0))*/ << "," << endl;
    }
    delete[] score;
    delete[] completeinfo;
    delete[] observedinfo;
    delete[] scoresq;
  }//end comploci loop
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

void MisSpecAlleleFreqTest::R_output3DarrayDimensions(ofstream* stream, const vector<int> dim, const vector<string> labels)
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
