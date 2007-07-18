// *-*-C++-*-*
/* 
 *   ScoreTestBase.h 
 *   Abstract Base Class for score tests
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef SCORETESTBASE_H
#define SCORETESTBASE_H 1

#include "common.h"
namespace bclib{
  class FileWriter;
}

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
  void OutputScalarScoreTest( int iterations, bclib::FileWriter& outputstream, std::string label,
				     const double score, const double scoresq, const double info, bool final);
  void OutputScoreTest( int iterations, bclib::FileWriter&, unsigned dim, std::vector<std::string> labels,
			       const double* score, const double* scoresq, const double* info, bool final, unsigned dim2);
  void OutputRaoBlackwellizedScoreTest( bclib::FileWriter&, std::string label,
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
#endif
