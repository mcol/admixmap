//=============================================================================
//
// Copyright (C) 2006, 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file ScoreTestBase.h
/// Definition of the ScoreTestBase class.
//=============================================================================

#ifndef SCORETESTBASE_H
#define SCORETESTBASE_H 1


#include <string>
#include <vector>

namespace bclib{
  class DelimitedFileWriter;
}


/** \addtogroup base
 * @{ */


/// Abstract Base Class for score tests
class ScoreTestBase{

public:
  //virtual interface
  ScoreTestBase();
  virtual ~ScoreTestBase();
  //virtual void Initialise() = 0;
  virtual void Reset() = 0;
  //virtual void Update() = 0;
  //virtual void Output() = 0;

  static std::string double2R(double);
  static std::string double2R(double, int);

protected:
  bool test;
  unsigned numPrintedIterations;///< number of times output is written to file (for R object dimension)
  unsigned numUpdates;///< number of times Update function has been called


  //utility functions
  void OutputScalarScoreTest( int iterations, bclib::DelimitedFileWriter& outputstream, std::string label,
				     const double score, const double scoresq, const double info, bool final);
  void OutputScoreTest( int iterations, bclib::DelimitedFileWriter&, unsigned dim, std::vector<std::string> labels,
			       const double* score, const double* scoresq, const double* info, bool final, unsigned dim2);
  void OutputRaoBlackwellizedScoreTest( bclib::DelimitedFileWriter&, std::string label,
					const double score, const double scoresq, const double varscore, 
					const double info, bool final );

};

// #include <exception>
// ///exception class to throw when trying to output test results when numUpdates = 0
// class ScoreTestNotUpdatedException : public std::exception{
// public:
//   ~ScoreTestNotUpdatedException() throw(){};
//   const char* what() const throw(){
//     return "Unable to output scoretest as no updates have been made";
//   }
// };


/** @} */


#endif
