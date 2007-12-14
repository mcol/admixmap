
#include "bclib/DataReader.h"
#include "bclib/StringConvertor.h"
#include <stdexcept>
#include <ctype.h>
#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>

BEGIN_BCLIB_NAMESPACE

StringSplitter DataReader::splitter;

DataReader::DataReader(){

}
DataReader::~DataReader(){

}

void DataReader::ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, LogWriter& Log, bool requireHeader ){
  readFile(filename, SMatrix, Log, requireHeader);
}
void DataReader::ReadData(const char* filename, DataMatrix& DMatrix, LogWriter& Log,size_t row0, size_t col0, size_t ncols, bool requireHeader ){
  std::vector< std::vector<std::string> > SMatrix;
  readFile(filename, SMatrix, Log, requireHeader);
  convertMatrix(SMatrix, DMatrix, row0, col0, ncols);
}
void DataReader::ReadData(const char* filename, std::vector< std::vector<std::string> >& SMatrix, DataMatrix& DMatrix, LogWriter& Log, bool requireHeader ){
  readFile(filename, SMatrix, Log, requireHeader);
  convertMatrix(SMatrix, DMatrix, 0, 0, 0);
}

///read a file into a string matrix
void DataReader::readFile(const char *fname, std::vector<std::vector<std::string> >& data, LogWriter &Log, bool requireHeader)
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

      //check for header by searching for alphabetic characters
      getline(in, line);
      if( requireHeader && find_if(line.begin(), line.end(), isalpha) == line.end() ){
	std::string errstring = "ERROR: No header found in file ";
	errstring.append(fname);
	throw(errstring);
      }

      do{
	//skip blank lines
	if (!StringConvertor::isWhiteLine(line.c_str())) {
	  //tokenise, splitting on whitespace
	  data.push_back(splitter.split(line.c_str()));

	  //throw exception if line has length different from header
	  if(data.size()>1 && data[data.size()-1].size() != data[0].size()){
	    std::stringstream errstring;
	    errstring << "Inconsistent row lengths in file "
		      << fname
		      << ". " << data[0].size() << " fields in header but "
		      << data[data.size()-1].size() << " fields in row " << data.size()-1 << std::endl;
	    throw errstring.str();
	  }
	}
	//read next line
      }while (getline(in, line)); 
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

void DataReader::ReadHeader(const char* filename, std::vector<std::string>& labels, bool skipfirstcol)
{
  using namespace boost;
  std::ifstream file(filename);
  std::string header;
  if(skipfirstcol) file >> header; //skip first column
  // now read rest of line
  getline(file, header);  // returns string with trailing carriage return from dos format
  if(header.find_last_of(13)==header.size()) {
    header.resize(header.size() - 1); // strip trailing carriage return
  }
  
  file.close();
  boost::trim(header); // remove leading whitespace left after skipping first column
  StringSplitter::Tokenize(header, labels, " \t");
}
END_BCLIB_NAMESPACE
