#include "RObjectWriter.h"
#include <sstream>
#include <numeric>
#include <iomanip>
#include <cmath>

using namespace::std;
#define COMMENT_CHAR " # "

//base type manipulators
ostream& comma(ostream& os){
  os << ", ";
  return os;
}
ostream& newline(ostream& os){
  os << "\n";
  return os;
}

//derived type manipulators
void comma(RObjectWriter& R){
  if(R.needComma)
    R << COMMA;
  //will need a comma next time
  R.needComma = true;
}

void newline(RObjectWriter& R){
  //print a comma if this is not the first line
  if(R.needNewLineComma)
    R << COMMA;
  //print a comment on previous line if there is one
  if(R.comment_cache.size()){
    R << COMMENT_CHAR << R.comment_cache;
    R.comment_cache.clear();
  }
  //start new line
  R << "\n";

  //will not need a comma before first element of next line
  R.needComma = false;
  //will need a comma before next line (if there is one)
  R.needNewLineComma = true;
}

//stream insertion friends

RObjectWriter& operator<<( RObjectWriter & R, void (*manip)(RObjectWriter& R)){
  manip(R);
  return R;
}

//class members

RObjectWriter::RObjectWriter() {
  needNewLineComma = false;
  needComma = false;
}

RObjectWriter::RObjectWriter(const char* filename){
  needNewLineComma = false;
  needComma = false;
  open(filename);
}

RObjectWriter::RObjectWriter(const string& filename){
  needNewLineComma = false;
  needComma = false;
  open(filename.c_str());
}

void RObjectWriter::open(const char* filename){
  ((ofstream*)this)->open(filename);
  if(!this->is_open()){
    string error_string = "ERROR: could not open ";
    error_string.append(filename);
    throw(error_string);
  }

  //start writing R object
  *this << "structure(.Data=c(" ;//<< endl;
}

void RObjectWriter::open(const string& filename){
  open(filename.c_str());
}

void RObjectWriter::close(const vector<int>& dim, const vector<vector<string> >& DimNames){
  WriteDimensions(dim, DimNames);
  ((ofstream*)this)->close();
}

//private function to prevent closing without writing dimensions
void RObjectWriter::close(){
  ((ofstream*)this)->close();
}

RObjectWriter& RObjectWriter::comment(const char* s){
  if(needComma){//comment is not on new line
    //cache comment for output later, in case this is the last line
    //strcpy(comment_cache, s);
    comment_cache = s;
  }
  else{//on new line
    *this << COMMENT_CHAR << s ;
    needNewLineComma = false;
  }

  //no need for comma after comment
  needComma = false;

  return *this;
}

void RObjectWriter::WriteDimensions(const vector<int>& dim, const vector<vector<string> >& DimNames)
{
  //sanity checks
  if(dim.size() < DimNames.size())
    throw string("Error in RobjectWriter::WriteDimensions : dimensions do not match dimnames");

  //write comment for last line if there is one
  if(comment_cache.size()){
    *this << COMMENT_CHAR << comment_cache;
    comment_cache.clear();
  }

    //finish defining data 
  *this << endl << ")," << endl

  //define dimensions
	<< ".Dim = c(";

  for(unsigned int i = 0; i < dim.size(); i++){
    *this << dim[i];
    if(i != dim.size() - 1){
      *this << ",";
    }
  }
  *this << ")," << endl
  
  //define dimnames
	<< ".Dimnames=list(";
  
  unsigned d = 0;
  for(; d < DimNames.size(); ++d){
    if(dim[d] != (int)DimNames[d].size())
      throw string("Error in RobjectWriter::WriteDimensions : dimnames do not match dimension");

    if(d > 0)//seperate dimnames with comma
      *this << ", ";    

    //write dimnames for labeled dimensions
    *this << "c(";
    for(unsigned int i = 0; i < DimNames[d].size(); i++){
      *this << "\"" << DimNames[d][i] << "\"";
      if(i != DimNames[d].size() - 1){
	*this << ",";
      }
    }
    *this << ")";
  }

  //write 'character(0)' for remaining dimensions
  for( ; d < dim.size(); ++d){
    if(d > 0)//seperate dimnames with comma
      *this << ", ";    

    *this << "character(0)";
  }
  //...and finally
  *this << "))" << endl;

}



