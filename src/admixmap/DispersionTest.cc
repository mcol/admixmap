//=============================================================================
//
// Copyright (C) 2005-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file DispersionTest.cc
/// Implementation of the DispersionTest class.
//=============================================================================

#include "DispersionTest.h"
#include "AdmixFilenames.h"
#include "bclib/dist.h"
#include "bclib/linalg.h"
#include "bclib/LogWriter.h"
#include "bclib/rand.h"
#include <cmath>
#include <numeric> // for accumulate()

using namespace ::std;

DispersionTest::DispersionTest(){
  divergentallelefreqstest = 0;
  NumberOfCompositeLoci = -1;
  NumberOfPopulations = 0;
}

void DispersionTest::Initialise(const string& resultsDir, bclib::LogWriter &Log, int NumLoci, int NumPopulations){
  NumberOfCompositeLoci = NumLoci;
  NumberOfPopulations = NumPopulations;

  divergentallelefreqstest = bclib::alloc2D_i(NumberOfCompositeLoci + 1, NumberOfPopulations );
  
  Log.setDisplayMode(bclib::Quiet);
  const string filename = resultsDir + "/" + DISPERSION_TEST_FILE;
  Log << "Writing dispersion test results to " << filename << "\n";
  dispersionoutputstream.open( filename.c_str(), ios::out );
  if( !dispersionoutputstream.is_open() ){
    Log.setDisplayMode(bclib::On);
    Log << "ERROR: Couldn't open dispersiontestfile\n";
    exit( 1 );
  }

}

DispersionTest::~DispersionTest(){
  bclib::free_matrix(divergentallelefreqstest, NumberOfCompositeLoci + 1);
}

void DispersionTest::TestForDivergentAlleleFrequencies(const AlleleFreqs* const A, const IndividualCollection* const IC)
{
  vector<int> rep;
  vector<double> popfreqs;

  vector<int> AlleleCount; 
  vector<double> sum(NumberOfPopulations, 0.0);
  vector<double> repsum(NumberOfPopulations, 0.0);

  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    for( int k = 0; k < NumberOfPopulations; k++ ){
      A->GetAlleleFreqs(popfreqs, j, k);
      AlleleCount.resize(popfreqs.size());

      // Generate replicate data conditional on locus ancestry
      IC->getAlleleCounts(AlleleCount, j, k);
      int sumcounts = accumulate(AlleleCount.begin(), AlleleCount.end(), 0, plus<int>());
      bclib::Rand::genmultinomial(sumcounts, popfreqs, rep);

      // Calculate likelihood of observed and replicate data
      double LogLikelihood = log(bclib::MultinomialPDF(AlleleCount, popfreqs));
      double RepLogLikelihood =	log(bclib::MultinomialPDF(rep, popfreqs));
      if ( !(LogLikelihood < RepLogLikelihood) )
	divergentallelefreqstest[j][k] ++;
      sum[k] += LogLikelihood;
      repsum[k] += RepLogLikelihood;
    }
  }

  //test statistic for population
  for( int k = 0; k < NumberOfPopulations; k++ ){
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
