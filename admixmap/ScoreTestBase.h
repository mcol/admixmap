// *-*-C++-*-*
#ifndef SCORETESTBASE_H
#define SCORETESTBASE_H 1

//#include <sstream>
//#include "IndividualCollection.h"
#include "utils/LogWriter.h"
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

protected:
  bool test;
  std::ofstream outputfile;
  //utility functions
  static void OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname, bool Robj);
  static void OutputScalarScoreTest( int iterations, std::ofstream* outputstream, std::string label,
				     const double score, const double scoresq, const double info, bool final);
  static void OutputScoreTest( int iterations, std::ofstream* outputstream, unsigned dim, std::vector<std::string> labels,
			       const double* score, const double* scoresq, const double* info, bool final, unsigned dim2);
  static void OutputRaoBlackwellizedScoreTest( int iterations, std::ofstream* outputstream, std::string label,
					       const double score, const double scoresq, const double varscore, 
					       const double info, bool final );

  static std::string double2R(double);
  static std::string double2R(double, int);
  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);

};


#endif
