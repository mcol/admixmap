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
#include "AdmixFilenames.h"
#include "bclib/LogWriter.h"
#include "bclib/DelimitedFileWriter.h"
#include "AlleleFreqs.h"
#include "Genome.h"
#include "gsl/gsl_cdf.h"//for pvalues
#include "bclib/linalg.h"

MisSpecAlleleFreqTest::MisSpecAlleleFreqTest(){
  doTest1 = false;
  doTest2 = false;
}

MisSpecAlleleFreqTest::~MisSpecAlleleFreqTest(){
}


void MisSpecAlleleFreqTest::Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter &Log )
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

void MisSpecAlleleFreqTest::Update(const IndividualCollection* const individuals,
				   const AlleleFreqs* const A, const Genome* const Loci){
  if(doTest1)
    Test1.Update(individuals, A, Loci);

  if( doTest2 )
    Test2.Update(individuals, A, Loci);
}

void MisSpecAlleleFreqTest::Output(const string& ResultsDir, const Genome* const Loci,
				   const Vector_s& PopLabels, bclib::LogWriter& Log){
    if( doTest1){
      const string filename = ResultsDir + "/" + MISSPECALLELEFREQTEST_1;
      Test1.Output(filename.c_str(), Loci, PopLabels, Log);
    }
    if( doTest2 ){
      const string filename = ResultsDir + "/" + MISSPECALLELEFREQTEST_2;
      Test2.Output(filename.c_str(), Loci, PopLabels, Log);
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
  numUpdates = 0;
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


void MisSpecifiedAlleleFreqTest::Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter& Log)
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
    double** phi = bclib::alloc2D_d(Populations, Populations);
    for( int i = 0; i < individuals->getSize(); i++ ){
      const PedBase & ind = individuals->getElement(i);

      for( int k = 0; k < Populations; k++ )
	for( int kk = 0; kk < Populations; kk++ )
	  phi[k][ kk ] = ind.getAdmixtureProps()[k] * ind.getAdmixtureProps()[kk];

      for(int j = 0; j < NumCompLoci; j++ ){
	if( !(ind.GenotypeIsMissing(j)) &&
	    (*Loci)(j)->GetNumberOfLoci() == 1  && !(A->IsRandom()) ){//CHECK: do only for SNPs?
	  int NumStates = Loci->GetNumberOfStates(j);
	  const vector<int> NumCopiesAllele1 = (*Loci)(j)->getAlleleCounts(1, ind.getSampledHapPair(j));
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
      bclib::matrix_product(Score[j], ScoreSq, Populations, 1);
      transform(ScoreSq, ScoreSq+Populations*Populations, SumScoreSq[j], SumScoreSq[j], std::plus<double>());
      delete[] ScoreSq;
    }
    bclib::free_matrix(phi, Populations);
  }
  ++numUpdates;
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

void MisSpecifiedAlleleFreqTest::Output( const char* filename, const Genome* const Loci,
					 const Vector_s& PopLabels, bclib::LogWriter& Log)
{
  if( test){

     Log << "Writing score tests for mis-specified allele frequencies(1) to "
	 << filename << "\n";
     bclib::DelimitedFileWriter outputfile(filename);
     outputfile << "Locus" << "Population" << "Score" << "CompleteInfo"
		<< "ObservedInfo" << "PercentInfo" << "StdNormal" << "PValue"
		<< "ChiSquared" << "PValue"<< bclib::newline;

    double* ObservedMatrix = new double[Populations*Populations];
    double* ScoreSq = new double[Populations*Populations];
    outputfile.setDecimalPrecision(2);//output 2 decimal places
    for(int j = 0; j < NumCompLoci; j++ ){
      if( (*Loci)(j)->GetNumberOfLoci() == 1 ){

	for(int k = 0; k < Populations; ++k)SumScore[j][k] /= (double) numUpdates;
	bclib::matrix_product(SumScore[j], ScoreSq, Populations, 1);//ScoreSq = Score * t(Score)
	//CompleteMatrix = SumInfo[j] / samples;
	for(int k = 0; k < Populations*Populations; ++k){
	  SumInfo[j][k] /= (double) numUpdates;//SumInfo[j] = CompleteMatrix
	  ObservedMatrix[k] = SumInfo[j][k] + ScoreSq[k] - SumScoreSq[j][k] / (double)numUpdates;
	}

	//ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumScoreSq[j] / numUpdates;

	for( int k = 0; k < Populations; k++ ){
	  double observedinfo = ObservedMatrix[ k*Populations + k ];
	  double completeinfo = SumInfo[j][ k*Populations + k ];
	  // Test for mis-specification within each continental-population.
	  outputfile << (*Loci)(j)->GetLabel(0)
		     << PopLabels[k]
		     << SumScore[j][ k ] //score
		     << completeinfo     //complete info
		     << observedinfo ;   //observed info
	  if(completeinfo > 0.0)
	    outputfile << 100*observedinfo / completeinfo ;
	  else outputfile << "NA";

	  if(observedinfo > 0.0){
	    double zscore = SumScore[j][ k ] / sqrt( observedinfo );
	    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
	    outputfile << zscore << pvalue ;
	  }
	  else outputfile << "NA" << "NA";

	  if( k < Populations - 1 )
	    outputfile  << "NA" << "NA" << bclib::newline;//output NA in chisq column
	}
	// Summary chi-sq test for mis-specification in all continental-populations.
	try{
	  double chisq = bclib::GaussianMarginalQuadraticForm(Populations, SumScore[j], ObservedMatrix, Populations);
	  double pvalue = 1.0 - gsl_cdf_chisq_P (chisq, 2.0);
	  outputfile << chisq << pvalue << bclib::newline;
	}
	catch(...){
	  outputfile << "NA" << "NA" << bclib::newline;
	}
      }
    }//end locus loop
    delete[] ObservedMatrix;
    delete[] ScoreSq;
    outputfile.close();
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
  numUpdates = 0;
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


void MisSpecifiedAlleleFreqTest2::Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter &Log )
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
   }
 }
 else{
   Log << "ERROR: cannot test for mis-specification of allele frequencies unless allelefrequencies are fixed\n"
       << "This option will be ignored\n";
   //TODO: set indicators in AdmixOptions to false, remove keys from option map.
 }
}

void MisSpecifiedAlleleFreqTest2::Update(const IndividualCollection* const individuals,
					 const AlleleFreqs* const A, const Genome* const Loci){
  if( test ){
    for( int j = 0; j < NumCompLoci; j++ ){
      for( int i = 0; i < individuals->getSize(); i++ ){
	UpdateScoreForMisSpecOfAlleleFreqs2(j, (*Loci)(j)->GetNumberOfStates(), A->GetAlleleFreqs(j), individuals);
      }
    }
    ++numUpdates;
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
      bclib::matrix_product(NewScore, NewScoreSq, NumberOfStates-1, 1);
      transform(NewScoreSq, NewScoreSq+(NumberOfStates-1)*(NumberOfStates-1), SumScoreSq[locus][k], SumScoreSq[locus][k], std::plus<double>());
      delete[] NewScoreSq;
   }
   delete[] NewScore;
   delete[] NewInfo;
}

void MisSpecifiedAlleleFreqTest2::Output(const char* filename, const Genome* const Loci,
					 const Vector_s& PopLabels, bclib::LogWriter& Log)
{
  if( test ){
     Log << "Writing score tests for mis-specified allele frequencies(2) to "
	 << filename << "\n";
     bclib::DelimitedFileWriter outputfile(filename);
     outputfile << "Locus" << "Population" << "CompleteInfo" << "ObservedInfo" << "PercentInfo" << "ChiSquared" << bclib::newline;

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
	  score[s] = SumScore[j][k][s] / (double)numUpdates;
	}

	bclib::matrix_product(score, scoresq, NumberOfStates-1, 1);//scoresq = score * t(score)
	for(int s = 0; s < (NumberOfStates-1)*(NumberOfStates-1); ++s){
	  completeinfo[s] = SumInfo[j][k][s] / (double)numUpdates;
	  observedinfo[s] = completeinfo[s] + scoresq[s] - SumScoreSq[j][k][s] / (double)numUpdates;
	}

	//score = SumScore[j][k] / numUpdates;
	//completeinfo = SumInfo[j][k] / numUpdates;
	//observedinfo = completeinfo + score * score.Transpose() - SumScoreSq[j][k] / numUpdates;
	outputfile << (*Loci)(j)->GetLabel(0)
		   << PopLabels[k] ;
	try{
	  double det1 = bclib::determinant(completeinfo, NumberOfStates-1);
	  double det2 = bclib::determinant(observedinfo, NumberOfStates-1);
	  outputfile << det1
		     << det2
		     << 100*det2 / det1 ;
	}
	catch(...){
	  outputfile << "NA" << "NA" << "NA";
	}
	try{
	  double chisq = bclib::GaussianMarginalQuadraticForm(NumberOfStates-1, score, observedinfo, NumberOfStates-1);
	  //observedinfo.InvertUsingLUDecomposition();
	  outputfile << chisq/*double2R((score.Transpose() * observedinfo * score)(0,0))*/  << bclib::newline;
	}
	catch(...){
	  outputfile << "NA" << bclib::newline;
	}

      }
      delete[] score;
      delete[] completeinfo;
      delete[] observedinfo;
      delete[] scoresq;
    }//end comploci loop
    outputfile.close();
  }//end if test
}
