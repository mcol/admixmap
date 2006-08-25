// *-*-C++-*-*
/*
  DataReader.h
  Utility class for reading data from file into a matrix of strings or a DataMatrix
*/
#ifndef DATAREADER_H
#define DATAREADER_H 1
#include <string>
#include "DataMatrix.h"
#include "StringSplitter.h"
#include "LogWriter.h"

///  Utility class for reading data from file into a matrix of strings or a DataMatrix
class DataReader{
public:    
  DataReader();
  ~DataReader();
  static void ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, LogWriter& Log );
  static void ReadData(const char* filename, DataMatrix& DMatrix, LogWriter& Log );
  static void ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, DataMatrix& DMatrix, LogWriter& Log );

  static void convertMatrix(const std::vector<std::vector<std::string> >& data, DataMatrix& m, 
			    size_t row0 = 0, size_t col0 = 0, size_t ncols = 0);
private:
  static StringSplitter splitter;

  static void readFile(const char *fname, Matrix_s& data, LogWriter &Log);
};

#endif /* !defined DATAREADER_H */
