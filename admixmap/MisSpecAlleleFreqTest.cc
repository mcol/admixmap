/** 
 *   ADMIXMAP
 *   MisSpecAlleleFreqTest.cc 
 *   Implements score tests for mis-specified allele frequencies
 *   Only used with fixed allele frequencies and only for monitoring. 

 * There are currently two tests:
 * 1. scalar test carried out only at composite loci with a single locus
 * 2. a generalized test , not fully implemented yet
 *
 * Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "MisSpecAlleleFreqTest.h"
#include "gsl/gsl_cdf.h"//for pvalues
#include "linalg.h"

MisSpecAlleleFreqTest::MisSpecAlleleFreqTest(){
  doTest1 = false;
  doTest2 = false;
}

MisSpecAlleleFreqTest::~MisSpecAlleleFreqTest(){
}


void MisSpecAlleleFreqTest::Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log )
{
 if( strlen( options->getAlleleFreqFilename() ) || options->getFixedAlleleFreqs()){

   doTest1 = options->getTestForMisspecifiedAlleleFreqs();
   doTest2 = options->getTestForMisspecifiedAlleleFreqs2();

   
   if(doTest1){
     Test1.Initialise(options, Loci, Log);
   }

   if(doTest2){
     Test2.Initialise(options, Loci, Log);
   }
 }
 else{
   Log << "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n"
       << "This option will be ignored\n";
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecAlleleFreqTest::Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci){
  if(doTest1)
    Test1.Update(individuals, A, Loci);

  if( doTest2 )    
    Test2.Update(individuals, A, Loci);
}

void MisSpecAlleleFreqTest::Output(int samples, const Genome* const Loci, const Vector_s& PopLabels){
    if( doTest1){
      Test1.Output(samples, Loci, PopLabels);
    }
    if( doTest2 ){
      Test2.Output(samples, Loci, PopLabels);
    }
}

//////////////////////////////////////////////////////////////////////////////
//                    scalar test for misspecified allele freqs (SNPs only)
MisSpecifiedAlleleFreqTest::MisSpecifiedAlleleFreqTest(){
  Score = 0;
  SumScore = 0;
  Info = 0;
  SumInfo = 0;
  SumScoreSq = 0;

  NumTestLoci = 0;
  NumCompLoci = 0;
  test = false;
  Populations = 1;
}

MisSpecifiedAlleleFreqTest::~MisSpecifiedAlleleFreqTest(){
  if(test){
    for(int i = NumCompLoci; i > 0; --i){
      delete[] Score[i-1];
      delete[] Info[i-1];
      delete[] SumScore[i-1];
      delete[] SumInfo[i-1];
      delete[] SumScoreSq[i-1];  
    }
    delete[] Score;
    delete[] SumScore;
    delete[] Info;
    delete[] SumInfo;
    delete[] SumScoreSq;
  }
}


void MisSpecifiedAlleleFreqTest::Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log )
{
 if( strlen( options->getAlleleFreqFilename() ) || options->getFixedAlleleFreqs()){

   test = options->getTestForMisspecifiedAlleleFreqs();
   Populations = options->getPopulations();
   NumCompLoci = Loci->GetNumberOfCompositeLoci();
   
   if(test){
     NumTestLoci = 0;
     for(int i = 0; i < NumCompLoci; ++i)
       if((*Loci)(i)->GetNumberOfLoci() == 1) NumTestLoci++;

     Score = new double*[NumCompLoci];
     SumScore = new double*[NumCompLoci];
     Info = new double*[NumCompLoci];
     SumInfo = new double*[NumCompLoci];
     SumScoreSq = new double*[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       Score[i] = new double[ Populations ];
       Info[i] = new double[ Populations * Populations ]; 
       SumScore[i] = new double[ Populations ];
       fill(SumScore[i], SumScore[i]+Populations, 0.0);
       SumInfo[i] = new double[ Populations * Populations ];
       fill(SumInfo[i], SumInfo[i]+ Populations*Populations, 0.0);
       SumScoreSq[i] = new double[ Populations * Populations ];
       fill(SumScoreSq[i], SumScoreSq[i]+Populations*Populations, 0.0);
    }
     outputfile.open( options->getAlleleFreqScoreFilename() );
     if(!outputfile.is_open()){
       string error_string = "ERROR: could not open ";
       error_string.append(options->getAlleleFreqScoreFilename());
       throw(error_string);
     }
     Log << "Writing score tests for mis-specified allele frequencies(1) to " << options->getAlleleFreqScoreFilename() << "\n";
     outputfile << "Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tsPValue\tChiSquared\tPValue"<< endl;
   }

 }
 else{
   Log << "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n"
       << "This option will be ignored\n";
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecifiedAlleleFreqTest::Reset(){
  if(test)
    for(int j = 0; j < NumCompLoci; ++j){
      fill(Score[j], Score[j]+Populations, 0.0);
      fill(Info[j], Info[j]+Populations*Populations, 0.0);
    }
}

void MisSpecifiedAlleleFreqTest::Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci){

  Reset();
  if( test ) {
    double** phi = alloc2D_d(Populations, Populations);
    for( int i = 0; i < individuals->getSize(); i++ ){
      const Individual* ind = individuals->getIndividual(i);
      
      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  phi[k][ kk ] = ind->getAdmixtureProps()[k] * ind->getAdmixtureProps()[kk];
      
      for(int j = 0; j < NumCompLoci; j++ ){
	if( !(ind->GenotypeIsMissing(j)) && 
	    (*Loci)(j)->GetNumberOfLoci() == 1  && !(A->IsRandom()) ){//CHECK: do only for SNPs?
	  int NumStates = Loci->GetNumberOfStates(j);
	  const vector<int> NumCopiesAllele1 = (*Loci)(j)->getAlleleCounts(1, ind->getSampledHapPair(j));
	  UpdateLocus( j, phi, NumCopiesAllele1[0], A->GetAlleleFreqs(j), NumStates );
	}
      }
    }

    //accumulate over iterations
    for(int j = 0; j < NumCompLoci; j++ ){
      //SumScore[j] += Score[j];
      transform(Score[j], Score[j]+Populations, SumScore[j], SumScore[j], std::plus<double>());
      //SumInfo[j] += Info[j];
      transform(Info[j], Info[j]+Populations*Populations, SumInfo[j], SumInfo[j], std::plus<double>());
      //SumScoreSq[j] += Score[j] * Score[j].Transpose();
      double* ScoreSq = new double[Populations*Populations];
      matrix_product(Score[j], ScoreSq, Populations, 1);
      transform(ScoreSq, ScoreSq+Populations*Populations, SumScoreSq[j], SumScoreSq[j], std::plus<double>());
      delete[] ScoreSq;
    }
    free_matrix(phi, Populations);
  }

}

/**
 * N.B. This only works for a single SNiP.
 * Updates what's required for the score tests. Only used with fixed
 * allele frequencies. This function is only used for monitoring.
 * Here, the number of states is 2
 */
void MisSpecifiedAlleleFreqTest::UpdateLocus(int j, const double* const* phi, int NumCopiesAllele1,
					     const double* const AlleleFreqs, int NumStates)
{
  vector<double> score( Populations );
  double Pi[3] = {0.0, 0.0, 0.0};
  
  for( int k = 0; k < Populations; k++ ){
    for( int kk = 0; kk < Populations; kk++ ){
      Pi[0] += AlleleFreqs[ NumStates*k ] * AlleleFreqs[ NumStates*kk ] * phi[k][kk];
      Pi[1] += ( ( 1 - AlleleFreqs[ NumStates*k ] ) * AlleleFreqs[ NumStates*kk ] + ( 1 - AlleleFreqs[ NumStates*kk ] ) * AlleleFreqs[ NumStates*k ] ) *phi[k][kk];
      Pi[2] += ( 1 - AlleleFreqs[ NumStates*k ] ) * ( 1 - AlleleFreqs[ NumStates*kk ] ) * phi[k][kk];
    }
  }
  
  if( NumCopiesAllele1 == 2 ){
    for( int k = 0; k < Populations; k++ ){
      score[k] = 2 * AlleleFreqs[ NumStates*k ] * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  score[k] += AlleleFreqs[ NumStates*kk ] * (phi[k][kk] + phi[kk][k]);
      score[k] /= Pi[0];
      Score[j][ k ] += score[k];
      Info[j][ k*Populations + k ] += score[k] * score[k] - 2 * phi[k][k] / Pi[0];
    }
    for( int k = 0; k < Populations; k++ )
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  Info[j][ k*Populations + kk ] += score[k] * score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[0];}
  
  else if( NumCopiesAllele1 == 1 ){
    for( int k = 0; k < Populations; k++ ){
      score[k] = 2 * ( 1 - 2 * AlleleFreqs[ NumStates*k ] ) * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  score[k] += ( 1 - 2 * AlleleFreqs[ NumStates*kk ] ) * (phi[k][kk] + phi[kk][k]);
         score[k] /= Pi[1];
         Score[j][ k ] += score[k];
         Info[j][ k*Populations + k ] += score[k] * score[k] + 4 * phi[k][k] / Pi[1];}
       for( int k = 0; k < Populations; k++ )
	 for( int kk = 0; kk < Populations; kk++ )
	   if( k != kk )
	     Info[j][ k*Populations + kk ] += score[k] * score[kk] + 2*(phi[k][kk] + phi[kk][k]) / Pi[1];}
  
  else if( NumCopiesAllele1 == 0 ){
    for( int k = 0; k < Populations; k++ ){
      score[k] = -2 * ( 1 - AlleleFreqs[ NumStates*k ] ) * phi[k][k];
      for( int kk = 0; kk < Populations; kk++ )
	if( k != kk )
	  score[k] -= ( 1 - AlleleFreqs[ NumStates*kk ] ) * (phi[k][kk] + phi[kk][k]);
      score[k] /= Pi[2];
      Score[j][ k ] += score[k];
      Info[j][ k*Populations + k ] += score[k] * score[k] - 2 * phi[k][k] / Pi[2];}
      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  if( k != kk )
	    Info[j][ k*Populations + kk ] += score[k] * score[kk] - (phi[k][kk] + phi[kk][k]) / Pi[2];}
}

void MisSpecifiedAlleleFreqTest::Output( int samples, const Genome* const Loci, const Vector_s& PopLabels)
{
  if( test){
    double* ObservedMatrix = new double[Populations*Populations];
    double* ScoreSq = new double[Populations*Populations];
    outputfile << setiosflags(ios::fixed) << setprecision(2);//output 2 decimal places
    for(int j = 0; j < NumCompLoci; j++ ){
      if( (*Loci)(j)->GetNumberOfLoci() == 1 ){
	
	for(int k = 0; k < Populations; ++k)SumScore[j][k] /= (double) samples;
	matrix_product(SumScore[j], ScoreSq, Populations, 1);//ScoreSq = Score * t(Score)
	//CompleteMatrix = SumInfo[j] / samples;
	for(int k = 0; k < Populations*Populations; ++k){
	  SumInfo[j][k] /= (double) samples;//SumInfo[j] = CompleteMatrix
	  ObservedMatrix[k] = SumInfo[j][k] + ScoreSq[k] - SumScoreSq[j][k] / (double)samples;
	}
	
	//ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumScoreSq[j] / samples;
	
	for( int k = 0; k < Populations; k++ ){
	  double observedinfo = ObservedMatrix[ k*Populations + k ];
	  double completeinfo = SumInfo[j][ k*Populations + k ];
	  // Test for mis-specification within each continental-population.
	  outputfile << (*Loci)(j)->GetLabel(0) << "\t"
		     << PopLabels[k] << "\t"
		     << SumScore[j][ k ]  << "\t" //score
		     << completeinfo  << "\t"//complete info
		     << observedinfo  << "\t";//observed info
	  if(completeinfo > 0.0)
	    outputfile << 100*observedinfo / completeinfo  << "\t";
	  else outputfile << "NA\t";
	  
	  if(observedinfo > 0.0){
	    double zscore = SumScore[j][ k ] / sqrt( observedinfo );
	    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
	    outputfile << zscore << "\t" << pvalue  << "\t";
	  }
	  else outputfile << "NA\tNA\t";
	  
	  if( k < Populations - 1 )
	    outputfile  << "NA\tNA" << endl;//output NA in chisq column
	}
	// Summary chi-sq test for mis-specification in all continental-populations.
	try{
	  double chisq = GaussianMarginalQuadraticForm(Populations, SumScore[j], ObservedMatrix, Populations);
	  double pvalue = 1.0 - gsl_cdf_chisq_P (chisq, 2.0);
	  outputfile << chisq<< "\t" << pvalue << endl;
	}
	catch(...){
	  outputfile << "NA\tNA" << endl;
	}
      }
    }//end locus loop
    delete[] ObservedMatrix;
    delete[] ScoreSq; 
  }//end if test
}

////////////////////////////////////////////////////////////////////////////////
//                       Vector test for misspecified allele freqs
MisSpecifiedAlleleFreqTest2::MisSpecifiedAlleleFreqTest2(){
  SumScore = 0;
  SumInfo = 0;
  SumScoreSq = 0;
  NumCompLoci = 0;
  test = false;
  Populations = 1;
}

MisSpecifiedAlleleFreqTest2::~MisSpecifiedAlleleFreqTest2(){
  if(test){
    for(int i = NumCompLoci; i > 0; --i){
      for(int k = Populations; k >0; --k){
	delete[] SumScore[i-1][k-1];
 	delete[] SumInfo[i-1][k-1];
	delete[] SumScoreSq[i-1][k-1];
      }

      delete[] SumScore[i-1];
      delete[] SumInfo[i-1];
      delete[] SumScoreSq[i-1];
    }
    delete[] SumScore;
    delete[] SumInfo;
    delete[] SumScoreSq;
  }
}


void MisSpecifiedAlleleFreqTest2::Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log )
{
 if( strlen( options->getAlleleFreqFilename() ) || options->getFixedAlleleFreqs()){

   test = options->getTestForMisspecifiedAlleleFreqs2();
   Populations = options->getPopulations();
   NumCompLoci = Loci->GetNumberOfCompositeLoci();
   
   if(test){
     SumScore = new double**[NumCompLoci];
     SumInfo = new double**[NumCompLoci];
     SumScoreSq = new double**[NumCompLoci];
     
     for(int i = 0; i < NumCompLoci; ++i){
       SumScore[i] = new double*[Populations];
       SumInfo[i] = new double*[Populations];
       SumScoreSq[i] = new double*[Populations];

       int NumberOfStates = (*Loci)(i)->GetNumberOfStates();
       for(int k = 0; k < Populations; ++k){
	 SumScore[i][k] = new double[NumberOfStates - 1];
	 fill(SumScore[i][k], SumScore[i][k] + (NumberOfStates-1), 0.0);
	 SumInfo[i][k] = new double[ (NumberOfStates - 1)*( NumberOfStates - 1 )];
	 fill(SumInfo[i][k], SumInfo[i][k] + (NumberOfStates-1)*(NumberOfStates-1), 0.0);
	 SumScoreSq[i][k] = new double[(NumberOfStates - 1) * (NumberOfStates - 1 )];
	 fill(SumScoreSq[i][k], SumScoreSq[i][k] + (NumberOfStates-1)*(NumberOfStates-1), 0.0);
       }
     }
     outputfile.open( options->getAlleleFreqScoreFilename2() );
     if(!outputfile.is_open()){
       string error_string = "ERROR: could not open ";
       error_string.append(options->getAlleleFreqScoreFilename2());
       throw(error_string);
     }
     Log << "Writing score tests for mis-specified allele frequencies(2) to " << options->getAlleleFreqScoreFilename2() << "\n";
     outputfile << "Locus\tPopulation\tCompleteInfo\tObservedInfo\tPercentInfo\tChiSquared" << endl;
   }
 }
 else{
   Log << "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n"
       << "This option will be ignored\n";
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecifiedAlleleFreqTest2::Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci){
  if( test )    
    for( int j = 0; j < NumCompLoci; j++ ){
      for( int i = 0; i < individuals->getSize(); i++ ){
	UpdateScoreForMisSpecOfAlleleFreqs2(j, (*Loci)(j)->GetNumberOfStates(), A->GetAlleleFreqs(j), individuals);
      }
    }
}

// presumably this calculates score test for mis-spec allele freqs at multi-allelic loci
void MisSpecifiedAlleleFreqTest2::UpdateScoreForMisSpecOfAlleleFreqs2(const int locus, const int NumberOfStates, 
								const double* const AlleleFreqs, 
								const IndividualCollection* const individuals)
{
   double rn, r, pj, pi, q;
   double* NewScore = new double[ NumberOfStates - 1];
   double* NewInfo = new double[ (NumberOfStates - 1) * (NumberOfStates - 1 )];
   for( int k = 0; k < Populations; k++ ){
     vector<int> AlleleCounts = individuals->getAlleleCounts(locus, k, NumberOfStates);
     rn = (double)( AlleleCounts[ NumberOfStates - 1] );
     q = 1.0;
     for( int j = 0; j < NumberOfStates - 1; j++ )
       q -= AlleleFreqs[j + k*NumberOfStates];
      for( int j = 0; j < NumberOfStates - 1; j++ ){
         r = AlleleCounts[ j ];
         pj = AlleleFreqs[ j + k*NumberOfStates ];
         NewScore[ j ] = ( r / pj - rn / q ) * pj * ( 1 - pj );
         NewInfo[ j*(NumberOfStates-1) + j ] = pj * ( 1 - pj )
            * ( r - ( rn / q ) * ( 2*pj - 1.0 - pj / q + pj * pj / q ) );
         for( int i = j+1; i < NumberOfStates - 1; i++ ){
            pi = AlleleFreqs[i + k*NumberOfStates ];
            NewInfo[ j*(NumberOfStates-1) + i ] = rn * pj * ( 1 - pj ) * pi * ( 1 - pi ) / ( q * q );
            NewInfo[ i*(NumberOfStates-1) + j ] = NewInfo[ j*(NumberOfStates-1) + i ];
         }
      }
      //SumScore[locus][k] += NewScore;
      transform(NewScore, NewScore+(NumberOfStates-1), SumScore[locus][k], SumScore[locus][k], std::plus<double>());
      //SumInfo[locus][k] += NewInfo;
      transform(NewInfo, NewInfo+(NumberOfStates-1)*(NumberOfStates-1), SumInfo[locus][k], SumInfo[locus][k], std::plus<double>());
      //SumScoreSq[locus][k] += NewScore * NewScore.Transpose();
      double* NewScoreSq = new double[(NumberOfStates - 1) * (NumberOfStates - 1 )];
      matrix_product(NewScore, NewScoreSq, NumberOfStates-1, 1);
      transform(NewScoreSq, NewScoreSq+(NumberOfStates-1)*(NumberOfStates-1), SumScoreSq[locus][k], SumScoreSq[locus][k], std::plus<double>());
      delete[] NewScoreSq;
   }
   delete[] NewScore;
   delete[] NewInfo;
}

void MisSpecifiedAlleleFreqTest2::Output( int samples, const Genome* const Loci, const Vector_s& PopLabels)
{
  if( test ){
    //allelefreqscorestream2 << setiosflags(ios::fixed) << setprecision(2);//output 2 decimal places
    for(int j = 0; j < NumCompLoci; j++ ){
      int NumberOfStates = (*Loci)(j)->GetNumberOfStates();
      double* score = new double[NumberOfStates-1];
      double* completeinfo = new double[(NumberOfStates-1)*(NumberOfStates-1)];
      double* observedinfo = new double[(NumberOfStates-1)*(NumberOfStates-1)];
      double* scoresq = new double[(NumberOfStates-1)*(NumberOfStates-1)];
    
      //CompositeLocus *locus = (CompositeLocus*)(*Lociptr)(j);
      for( int k = 0; k < Populations; k++ ){
	for(int s = 0; s < NumberOfStates-1; ++s){
	  score[s] = SumScore[j][k][s] / (double)samples;
	}
	
	matrix_product(score, scoresq, NumberOfStates-1, 1);//scoresq = score * t(score)
	for(int s = 0; s < (NumberOfStates-1)*(NumberOfStates-1); ++s){
	  completeinfo[s] = SumInfo[j][k][s] / (double)samples;
	  observedinfo[s] = completeinfo[s] + scoresq[s] - SumScoreSq[j][k][s] / (double)samples;
	}
	
	//score = SumScore[j][k] / samples;
	//completeinfo = SumInfo[j][k] / samples;
	//observedinfo = completeinfo + score * score.Transpose() - SumScoreSq[j][k] / samples;
	outputfile << (*Loci)(j)->GetLabel(0) << "\t"
		   << PopLabels[k] << "\t";
	try{
	  double det1 = determinant(completeinfo, NumberOfStates-1);
	  double det2 = determinant(observedinfo, NumberOfStates-1);
	  outputfile << det1 << "\t"
		     << det2 << "\t"
		     << 100*det2 / det1 << "\t";
	}
	catch(...){
	  outputfile << "NA\tNA\tNA\t";
	}
	try{
	  double chisq = GaussianMarginalQuadraticForm(NumberOfStates-1, score, observedinfo, NumberOfStates-1);      
	  //observedinfo.InvertUsingLUDecomposition();
	  outputfile << chisq/*double2R((score.Transpose() * observedinfo * score)(0,0))*/  << endl;
	}
	catch(...){
	  outputfile <<"NA" << endl;
	}
	
      }
      delete[] score;
      delete[] completeinfo;
      delete[] observedinfo;
      delete[] scoresq;
    }//end comploci loop
  }//end if test
}




