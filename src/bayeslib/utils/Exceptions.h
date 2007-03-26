// *-*-C++-*-*
#ifndef BAYESLIB_EXCEPTIONS_H
#define BAYESLIB_EXCEPTIONS_H
#include <exception>
#include <string>
#include <sstream>

class InfiniteGradient : public std::exception{
 public:
  std::string function, file;
  InfiniteGradient(){};
  InfiniteGradient(const std::string& fun, const std::string f) :
    function(fun), file(f){};
  ~InfiniteGradient()throw(){};
  const char* what() const throw(){
    return (("Infinite gradient in function " + function + " in " + file).c_str());
  }

};

class nanException : public std::exception{
public:
  std::string description, file;
  int line;

  nanException(std::string d, const char* f, int l):
    description(d), file(f), line(l){};

  ~nanException()throw(){};
  const char* what() const throw (){
  std::stringstream message;

  message << file << ": " << line << " -- " << description;
  return message.str().c_str();
  }

};

#endif
