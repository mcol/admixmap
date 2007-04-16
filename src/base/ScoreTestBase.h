// *-*-C++-*-*
#ifndef SCORETESTBASE_H
#define SCORETESTBASE_H 1

//#include <sstream>
//#include "IndividualCollection.h"
#include "bcppcl/LogWriter.h"
#include "common.h"

/// Abstract Base Class for score tests
class ScoreTestBase{

public:
  //virtual interface
  virtual ~ScoreTestBase();
  //virtual void Initialise() = 0;
  virtual void Reset() = 0;
  //virtual void Update() = 0;
  //virtual void Output() = 0;

  static std::string double2R(double);
  static std::string double2R(double, int);
  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);

protected:
  bool test;
  std::ofstream outputfile;
  unsigned numPrintedIterations;///< number of times output is written to file (for R object dimension)
  unsigned numUpdates;///< number of times Update function has been called
  bool onFirstLine;///< to indicate whether to use a separator before newline

  //utility functions
  void OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname, bool Robj);
  void OutputScalarScoreTest( int iterations, std::ofstream* outputstream, std::string label,
				     const double score, const double scoresq, const double info, bool final);
  void OutputScoreTest( int iterations, std::ofstream* outputstream, unsigned dim, std::vector<std::string> labels,
			       const double* score, const double* scoresq, const double* info, bool final, unsigned dim2);
  void OutputRaoBlackwellizedScoreTest( std::ofstream* outputstream, std::string label,
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
