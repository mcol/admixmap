#include "bclib/StringConvertor.h"
#include <cctype>

using namespace std;

BEGIN_BCLIB_NAMESPACE

/// removes quotes from strings. 
string StringConvertor::dequote(const string& str)
{
  const string::size_type first = str.find_first_not_of("\"\'", 0);
  if(first != string::npos) {
    const string::size_type last = str.find_last_not_of("\"\'")+1;
    return str.substr(first, last-first);
  }
  else return str;
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

char LowerCase(char c){
  return (char)(tolower(c));
}

void StringConvertor::toLower(string& s){
  transform(s.begin(), s.end(), s.begin(), LowerCase);
}

void StringConvertor::RemoveCharFromString(string& s, const char c){

  string::size_type u = 0, t = 0;
  do{
     u = s.find(c, t);
     if(u != string::npos)
       s.erase(u, 1);

     t = u;
  }
  while(t != string::npos);

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

//TODO: use function templates to avoid duplication of these 
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
    // Convert elements to appropriate type and fill resulting vector.
    vec.push_back( (float)toFloat(str.substr(start, size)) );
    start = str.find_first_not_of(" ,", size);
    size = str.find_first_of(" ,", start);
  }
}

/**
   Converts a string to a vector of floats. 
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
void StringConvertor::StringToVec(const string s, std::vector<unsigned>& vec)
{
  string str = s;
  string::size_type start = str.find_first_not_of("\" "), size = 0;
  str = str.substr(start, str.find_first_of("\"", start)-1);
  start = 0;
  size += str.find_first_of(" ,", start);

  vec.clear();
  while(start != string::npos){
    // Convert elements to doubles and fill resulting vector.
    vec.push_back( (unsigned)toInt(str.substr(start, size)) );
    start = str.find_first_not_of(" ,", size);
    size = str.find_first_of(" ,", start);
  }
}

//small function to determine if a given string is in a list of strings
bool StringConvertor::isListedString(const std::string s, const std::vector<std::string>list){
  bool b = false;
  if(find(list.begin(), list.end(), s) < list.end())b = true;
  return b;
}
END_BCLIB_NAMESPACE
