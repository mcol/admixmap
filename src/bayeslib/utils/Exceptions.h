// *-*-C++-*-*
#ifndef BAYESLIB_EXCEPTIONS_H
#define BAYESLIB_EXCEPTIONS_H
#include <exception>
#include <string>

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

#endif
