// *-*-C++-*-*
/*
  DataReader.h
  Utility class for reading data from file into a matrix of strings or a DataMatrix
*/
#ifndef DATAREADER_H
#define DATAREADER_H 1
#include <string>
#include "bclib/bclib.h"
#include "bclib/DataMatrix.h"
#include "bclib/StringSplitter.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */


class LogWriter;


///  Utility class for reading data from file into a matrix of strings or a DataMatrix
class DataReader{
public:    
  DataReader();
  ~DataReader();
  static void ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, LogWriter& Log, bool requireHeader=true );
  static void ReadData(const char* filename, DataMatrix& DMatrix, LogWriter& Log, size_t row0 = 0, size_t col0 = 0, size_t ncols = 0, 
		       bool requireHeader=true );
  static void ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, DataMatrix& DMatrix, 
		       LogWriter& Log, bool requireHeader=true );

  static void convertMatrix(const std::vector<std::vector<std::string> >& data, DataMatrix& m, 
			    size_t row0 = 0, size_t col0 = 0, size_t ncols = 0);
  static void ReadHeader(const char* filename, std::vector<std::string>& labels, bool skipfirstcol = true);
private:
  static StringSplitter splitter;

  static void readFile(const char *fname, Matrix_s& data, LogWriter &Log, bool requireHeader=true);
};

/** @} */

END_BCLIB_NAMESPACE

#endif /* !defined DATAREADER_H */
