#include "bcppcl/DataReader.h"
#include "bcppcl/StringConvertor.h"
#include <stdexcept>
#include <ctype.h>
#include <algorithm>

StringSplitter DataReader::splitter;

DataReader::DataReader(){

}
DataReader::~DataReader(){

}

void DataReader::ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, LogWriter& Log ){
  readFile(filename, SMatrix, Log);
}
void DataReader::ReadData(const char* filename, DataMatrix& DMatrix, LogWriter& Log,size_t row0, size_t col0, size_t ncols ){
  std::vector< std::vector<std::string> > SMatrix;
  readFile(filename, SMatrix, Log);
  convertMatrix(SMatrix, DMatrix, row0, col0, ncols);
}
void DataReader::ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, DataMatrix& DMatrix, LogWriter& Log ){
  readFile(filename, SMatrix, Log);
  convertMatrix(SMatrix, DMatrix, 0, 0, 0);
}

///read a file into a string matrix
void DataReader::readFile(const char *fname, std::vector<std::vector<std::string> >& data, LogWriter &Log)
{
    if (0 == fname || 0 == strlen(fname)) return;

    std::ifstream in(fname);
    if (!in.is_open()) {
      std::string msg = "Cannot open file for reading: \"";
        msg += fname;
        msg += "\"";
        throw std::runtime_error(msg);
    }
    else {
      Log << "Loading " << fname << ".\n";
    }

    data.clear();
    try {
      std::string line; 

      //check for header
      getline(in, line);
      if( find_if(line.begin(), line.end(), isalpha) == line.end() ){
	std::string errstring = "ERROR: No header found in file ";
	errstring.append(fname);
	throw(errstring);
      }

      while (!in.eof()) {
	//skip blank lines
	if (!StringConvertor::isWhiteLine(line.c_str())) {
	  //tokenise, splitting on whitespace
	  data.push_back(splitter.split(line.c_str()));

	  //throw exception if line has length different from header
	  if(data.size()>1 && data[data.size()-1].size() != data[0].size()){
	    std::string errstring = "Inconsistent row lengths in file ";
	    errstring.append(fname);
	    throw errstring;
	  }
	}
	//read next line
	getline(in, line);
      }
    } catch (...) {
      in.close();
      throw;
    }
    
}

/**
 *  Auxilary function that converts a submatrix (starting at (row0, col0))of Matrix_s to DataMatrix
 */
void DataReader::convertMatrix(const std::vector<std::vector<std::string> >& data, DataMatrix& m, 
			       size_t row0, size_t col0, size_t ncols)
{       
    // If there are no rows, return empty matrix.
    if (0 == data.size()) return;

    const size_t numRows = data.size() - row0;

    // Verify that all rows have same length.
    const size_t totalnumCols = data[0].size();
    const size_t numCols = (ncols>0) ? ncols : totalnumCols - col0;
    for (size_t i = 1; i < numRows; ++i) {
        if (totalnumCols != data[i].size()) {
	  throw std::runtime_error("Invalid row length");
        }
    }
    
    // Form matrix.
    m.setDimensions(numRows, numCols);
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (StringConvertor::isMissingValue(data[i+row0][j+col0])) {
	      m.isMissing(i, j, true);
            } else {
	      m.set(i, j, StringConvertor::toFloat(data[i+row0][j+col0]));
            }
        }
    }
}

void DataReader::ReadHeader(const char* filename, std::vector<std::string>& labels, bool skipfirstcol){
  std::ifstream file(filename);
  std::string header;
  if(skipfirstcol)file >> header;//skip first column
  getline(file, header);
  file.close();
  StringSplitter::Tokenize(header, labels, " \t\"");
}
