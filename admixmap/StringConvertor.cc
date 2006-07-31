#include "StringConvertor.h"

using namespace std;

/// removes quotes from strings. 
string StringConvertor::dequote(const string& str)
{
  if (str.size() >= 2 && (str[0] == '"' || str[0]== '\'') && (str[str.size() - 1] == '"' || (str[str.size() - 1] == '\''))) {
        return str.substr(1, str.size() - 2);
      }

    return str;
}

/// converts a string to int
int StringConvertor::toInt(const string& str)
{
    return atoi(dequote(str).c_str());
}
  /// converts a string to float
double StringConvertor::toFloat(const string& str)
{
    return atof(dequote(str).c_str());
}

/** 
    Converts a string to a vector of up to unsigned short ints, splitting on separator sep.
    Works like Tokenize but splits into only two substrings. 
*/
void StringConvertor::toIntPair(std::vector<unsigned short>* a, const std::string& str, const char* sep)
{
  //strip quotes from string
  const std::string tmp = dequote(str);
  //look for a separator
  string::size_type i = tmp.find_first_of(sep);
  //extract part of string up to first instance of separator
  a->push_back(atoi(tmp.substr(0,i).c_str()));
  if( i != string::npos){//separator found
    //extract as int the part of the string after separator
    a->push_back(atoi(tmp.substr(i+1,tmp.length()-i).c_str()));
  }                      
}

/// determines if a char array consists only of whitespace.
bool StringConvertor::isWhiteLine(const char *p)
{
    while (*p) {
        if (!isspace(*p++)) {
            return false;
        }
    }

    return true;
}

/** determines if a string is a missing value.
    Missing values are denoted by #, ., or NA.
*/
bool StringConvertor::isMissingValue(const std::string& str)
{
    static string MISSING[] = {"#", ".", "NA"};
    return MISSING + 3 != find(MISSING, MISSING + 3, str);
}

/**
   Converts a string to a vector of floats. 
   Opening and closing quotes are removed and vector elements are separated by 
   commas or spaces.
*/
void StringConvertor::StringToVec(const string s, std::vector<float>& vec)
{
  string str = s;
  string::size_type start = str.find_first_not_of("\" "), size = 0;
  str = str.substr(start, str.find_first_of("\"", start)-1);
  start = 0;
  size += str.find_first_of(" ,", start);

  vec.clear();
  while(start != string::npos){
    // Convert elements to doubles and fill resulting vector.
    vec.push_back( toFloat(str.substr(start, size)) );
    start = str.find_first_not_of(" ,", size);
    size = str.find_first_of(" ,", start);
  }
}

/**
   Converts a atring to a vector of floats. 
   As previous function but converts to doubles rather than floats.
*/
void StringConvertor::StringToVec(const string s, std::vector<double>& vec)
{
  string str = s;
  string::size_type start = str.find_first_not_of("\" "), size = 0;
  str = str.substr(start, str.find_first_of("\"", start)-1);
  start = 0;
  size += str.find_first_of(" ,", start);

  vec.clear();
  while(start != string::npos){
    // Convert elements to doubles and fill resulting vector.
    vec.push_back( toFloat(str.substr(start, size)) );
    start = str.find_first_not_of(" ,", size);
    size = str.find_first_of(" ,", start);
  }
}
