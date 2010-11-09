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
/// \file StratificationTest.h
/// Definition of the StratificationTest class.
//=============================================================================

#ifndef STRATIFICATIONTEST
#define STRATIFICATIONTEST 1

#include <fstream>
#include <vector>
#include "IndividualCollection.h"
#include "Genome.h"

class AdmixOptions;
class FreqArray;

namespace bclib{
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


///Class to implement a test for residual population stratification
class StratificationTest
{
public:
  StratificationTest();
  
  void Initialize( AdmixOptions* const options, const Genome &Loci,  
		   const IndividualCollection* const IC, bclib::LogWriter &Log);
  void OpenOutputFile( const std::string& , bclib::LogWriter &);

  void calculate( const IndividualCollection* const individuals, const FreqArray& AlleleFreqs,
		  const std::vector<std::vector<int> > ChrmAndLocus, int Populations );

  void Output(bclib::LogWriter &);

private:
  int T;
  int count;
  int NumberOfTestLoci;
  std::vector<unsigned int> TestLoci;
  bool ModelIndicator;
  std::ofstream outputstream;

  std::vector<double>
  GenerateExpectedGenotype( const genepi::PedBase * const, const double*, const int  );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAdmixture(const std::vector<double>& ProbAllele1);

  std::vector<unsigned short>
  SimGenotypeConditionalOnAncestry( const double*, const int ancestry[2] );

  std::vector<unsigned short> SampleHeterozygotePhase( const double*, const int ancestry[2] );
};


/** @} */


#endif /* !defined STRATIFICATIONTEST_H */
