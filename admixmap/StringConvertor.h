// *-*-C++-*-*
#ifndef STRING_CONVERTOR_H
#define STRING_CONVERTOR_H 1

#include "common.h"

///  Utility class used to convert string to other types.
 
class StringConvertor
{
public:
  /// converts a string to int
  static int toInt(const std::string& str);
  /// converts a string to float
  static double toFloat(const std::string& str);
  
  //static std::pair<int, int> toIntPair(const std::string& str);
  /// converts a string to two unsigned ints. 
  static void toIntPair(unsigned int *,const std::string& str);
  /// converts a string to two unsigned short ints. 
  static void toIntPair(unsigned short *,const std::string& str);
  /// converts a string to a vector of two unsigned short ints. 
  static void toIntPair(std::vector<unsigned short>*, const std::string& str);
  /// determines if a string is a missing value.
  static bool isMissingValue(const std::string& str);
  /// determines if a char array consists only of whitespace.
  static bool isWhiteLine(const char *p);
  /// removes quotes from strings. 
  static std::string dequote(const std::string& str);
  /// Converts a C-style char array to a vector of doubles. 
  static void CstrToVec(const char* str, std::vector<double>&);
  /// Converts a C-style char array to a vector of floats. 
  static void CstrToVec(const char* str, std::vector<float>&);

  /// a dummy function
  static bool dummy();
};

#endif /* !STRING_CONVERTOR_H */
