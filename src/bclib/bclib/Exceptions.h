// *-*-C++-*-*
#ifndef BCLIB_EXCEPTIONS_H
#define BCLIB_EXCEPTIONS_H
#include "bclib/bclib.h"
#include <exception>
#include <string>
#include <sstream>

BEGIN_BCLIB_NAMESPACE

class InfinityException : public std::exception{
public:
  InfinityException();
  InfinityException(const std::string& fun, const std::string f) :
    function(fun), file(f){};
  ~InfinityException()throw(){};
  const char* what() const throw(){
    return (("Infinity in function \"" + function + "\" in " + file).c_str());
  }
private:
  std::string function, file;
};

class InfiniteGradient : public std::exception{
public:
  InfiniteGradient(){};
  InfiniteGradient(const std::string& fun, const std::string f) :
    function(fun), file(f){};
  ~InfiniteGradient()throw(){};
  const char* what() const throw(){
    return (("Infinite gradient in function " + function + " in " + file).c_str());
  }
private:
  std::string function, file;
};

class nanException : public std::exception{
public:
  nanException(std::string d, const char* f, int l):
    description(d), file(f), line(l){};

  ~nanException()throw(){};
  const char* what() const throw (){
  std::stringstream message;

  message << file << ": " << line << " -- " << description;
  return message.str().c_str();
  }
private:
  std::string description, file;
  int line;
  nanException();
};

/**
   exception class to throw when data from file are out of range.
   e.g. when you expect positive values and get a negative one.
*/
class DataOutOfRangeException : public std::exception{
public:
  /**
     Constructor.
     \param d name of whatever is being read
     \param r description of expected range, eg "between 0 and 1"
     \param f file being read
  */
  DataOutOfRangeException(std::string d, std::string r, std::string f) :
    dataName(d), range(r), file(f){};
  ~DataOutOfRangeException()throw(){};
  const char* what()const throw(){
    return (("ERROR: Bad values for " + dataName + " in " + file + ". Should be " + range).c_str());
  }

private:
  std::string dataName, range, file;
  //int line;

  DataOutOfRangeException();
};
END_BCLIB_NAMESPACE
#endif
