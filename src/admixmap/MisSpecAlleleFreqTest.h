//=============================================================================
//
// Copyright (C) 2005, 2006  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file MisSpecAlleleFreqTest.h
/// Definition of the MisSpecAlleleFreqTest class.
//=============================================================================

#ifndef MISSPECALLELEFREQTEST_H
#define MISSPECALLELEFREQTEST_H 1

#include "IndividualCollection.h"
#include "ScoreTestBase.h"

class AdmixOptions;
class AlleleFreqs;
class Genome;


/** \addtogroup admixmap
 * @{ */


///a scalar test for SNPs only
class MisSpecifiedAlleleFreqTest : public ScoreTestBase{
public:
  MisSpecifiedAlleleFreqTest();
  ~MisSpecifiedAlleleFreqTest();

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(const char* filename, const Genome* const Loci, const Vector_s& PopLabels, bclib::LogWriter& Log);
  void Reset();
  
private:
  double **Score;
  double **Info;
  double **SumScore;
  double **SumScoreSq;
  double **SumInfo;
  int NumTestLoci;     //number of comp loci with a single locus 
                       //ie those used in scalar score test for misspecified allelefreqs
  int NumCompLoci; //number of composite loci
  int Populations;  //number of populations
  
  void UpdateLocus(int j, const double* const* phi, int NumCopiesAllele1,
		   const double* const AlleleFreqs, int NumStates);
  
};

///a vector test, not fully tested
class MisSpecifiedAlleleFreqTest2: public ScoreTestBase{
public:
  MisSpecifiedAlleleFreqTest2();
  ~MisSpecifiedAlleleFreqTest2();
  
  void Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(const char* filename, const Genome* const Loci, const Vector_s& PopLabels, bclib::LogWriter& Log);
  void Reset(){};

private:
  int NumCompLoci; //number of composite loci
  int Populations;           //number of populations

  double ***SumScore;
  double ***SumScoreSq;
  double ***SumInfo;

  void UpdateScoreForMisSpecOfAlleleFreqs2(const int j, const int NumberOfStates, const double* const AlleleFreqs, 
					   const IndividualCollection* const individuals);

  void OutputTestsForMisSpecifiedAlleleFreqs2(const Genome* const Loci, const std::string* const PopLabels);

};

/// a container class for the tests for misspecified allele frequencies
class MisSpecAlleleFreqTest{

public:
  MisSpecAlleleFreqTest();
  ~MisSpecAlleleFreqTest();

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, bclib::LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(const string& ResultsDir, const Genome* const Loci, 
	      const Vector_s& PopLabels, bclib::LogWriter& Log);

private:
  bool doTest1, doTest2;//indicators for the two tests

  MisSpecifiedAlleleFreqTest Test1;
  MisSpecifiedAlleleFreqTest2 Test2;

};


/** @} */


#endif

