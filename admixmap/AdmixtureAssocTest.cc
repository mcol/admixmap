/** 
 *   ADMIXMAP
 *   AdmixtureAssocTest.cc 
 *   Class implements the score test for admixture association (admixturescoretest)
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AdmixtureAssocTest.h"
#include <math.h>//for sqrt
#ifdef PARALLEL
#include "Comms.h"
#endif

AdmixtureAssocTest::AdmixtureAssocTest(){
  Score = 0; 
  Info = 0; 
  SumScore = 0; 
  SumScore2 = 0;
  SumInfo = 0;
  NumPopulations = 0;
  NumOutcomeVars = 0;
}
AdmixtureAssocTest::~AdmixtureAssocTest(){
  //delete arrays for admixture assoc score test
  delete[] Score;
  delete[] Info;
  delete[] SumScore;
  delete[] SumInfo;
  delete[] SumScore2;
}

void AdmixtureAssocTest::Initialise(const unsigned K, const unsigned NumOutcomes, const char* filename,
				    const Vector_s& PLabels, LogWriter &Log){
  NumPopulations = K;
  NumOutcomeVars = NumOutcomes;
  test = true;
  /*----------------------
    | admixture association |
    -----------------------*/
  //TODO check conditions on this test
  if( test ){
    if ( strlen( filename ) ){
      //can't use this yet as this test doesn't write an R object
      //OpenFile(Log, &assocscorestream, options->getAssocScoreFilename(), "Tests for admixture association");
      outputfile.open( filename, ios::out );
      if( !outputfile ){
	Log << On <<"ERROR: Couldn't open admixturescorefile\n";
	exit( 1 );}
      else {
	Log << Quiet << "Writing tests for admixture association to: " << filename << "\n";
	outputfile << setiosflags( ios::fixed );
	
	Score = new double[NumPopulations * NumOutcomeVars];
	SumScore = new double[NumPopulations * NumOutcomeVars];
	SumScore2 = new double[NumPopulations * NumOutcomeVars];
	Info = new double[NumPopulations * NumOutcomeVars];
	SumInfo = new double[NumPopulations * NumOutcomeVars];
	fill(Score, Score + NumPopulations * NumOutcomeVars, 0.0);
	fill(Info, Info + NumPopulations * NumOutcomeVars, 0.0);
	fill(SumScore, SumScore + NumPopulations * NumOutcomeVars, 0.0);
	fill(SumScore2, SumScore2 + NumPopulations * NumOutcomeVars, 0.0);
	fill(SumInfo, SumInfo + NumPopulations * NumOutcomeVars, 0.0);
	InitialiseAssocScoreFile(PLabels);
#ifdef PARALLEL
	Comms::SetDoubleWorkspace(NumPopulations * NumOutcomeVars, false);
#endif

      }
    }
    else{
      Log << On << "No admixturescorefile given\n";
      exit(1);}
  }
}
void AdmixtureAssocTest::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.
  if( test ){
    fill(Score, Score + NumPopulations * NumOutcomeVars, 0.0);
    fill(Info, Info + NumPopulations * NumOutcomeVars, 0.0);
  }
}    

//Initialise ergodic average score file
void AdmixtureAssocTest::InitialiseAssocScoreFile(const Vector_s& PLabels){
  if( test ){
    //PopLabels = PLabels;
    outputfile << "Ergodic averages of score statistic for populations:\n";
    for( unsigned i = 0; i < NumPopulations; i++ ){
      outputfile << PLabels[i] << " ";
      if( !i )
	outputfile << " ";
    }
    outputfile << "\ncomplete  missing   statistic  ";
    for( unsigned i = 1; i < NumPopulations; i++ )
      outputfile << "complete  missing   statistic ";
    outputfile << endl;
  }
}

//Updates score and info for a single individual
void AdmixtureAssocTest::UpdateIndividualScore( const double* const Theta, double YMinusEY, double phi, double DInvLink, bool RandomMatingModel)
{
  double x;
  for( unsigned k = 0; k < NumPopulations; k++ ){
    if( RandomMatingModel )
      x = 0.5 * ( Theta[k] + Theta[ NumPopulations + k ] );
    else
      x = Theta[k];
    Score[ k*NumOutcomeVars] += phi * x * YMinusEY;// only for 1st outcomevar 
    Info[ k*NumOutcomeVars] += phi * x * x *DInvLink;
  }
}

void AdmixtureAssocTest::Update(){
  if( test ){
#ifdef PARALLEL
  //sum scores across processes
    Comms::ReduceAdmixtureAssocScores(Score, Info, NumPopulations * NumOutcomeVars);
#endif

    //SumScore += Score;
    transform(Score, Score + NumPopulations*NumOutcomeVars, SumScore, SumScore, std::plus<double>());
    //SumAdmixtureInfo += AdmixtureInfo;
    transform(Info, Info + NumPopulations*NumOutcomeVars, SumInfo, SumScore, std::plus<double>());
    
    for( unsigned k = 0; k < NumPopulations; k++ )
      for( unsigned j = 0; j < NumOutcomeVars; j++ )
	SumScore2[ k*NumOutcomeVars + j ] += Score[ k*NumOutcomeVars + j ] * Score[ k*NumOutcomeVars + j ];
  }
}
void AdmixtureAssocTest::Output(int iterations)
{
  for( unsigned k = 0; k < NumPopulations; k++ ){
    for( unsigned j = 0; j < NumOutcomeVars; j++ ){
      double EU = SumScore[ k*NumOutcomeVars + j ] / ( iterations );
      double complete = SumInfo[ k*NumOutcomeVars + j ] / ( iterations );
      double missing = SumScore2[ k*NumOutcomeVars + j ] / ( iterations ) - EU * EU;
      outputfile.width(9);
      outputfile << setprecision(6) << double2R(complete) << " ";
      outputfile.width(9);
      outputfile << setprecision(6) << double2R(missing) << " ";
      outputfile.width(9);
      outputfile << setprecision(6) << double2R(EU / sqrt( complete - missing )) << " ";
    }
  }
  outputfile << endl;
}
