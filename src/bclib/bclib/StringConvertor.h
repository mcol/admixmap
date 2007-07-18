// *-*-C++-*-*
#ifndef STRING_CONVERTOR_H
#define STRING_CONVERTOR_H 1
#include "bclib/bclib.h"
#include <string>
#include <vector>

BEGIN_BCLIB_NAMESPACE
///  Utility class used to convert string to other types.
class StringConvertor
{
public:
  /// converts a string to int
  static int toInt(const std::string& str);
  /// converts a string to float
  static double toFloat(const std::string& str);

  /// converts a string to a vector of two unsigned short ints. 
  static void toIntPair(std::vector<unsigned short>*, const std::string& str, const char* sep);
  /// determines if a string is a missing value.
  static bool isMissingValue(const std::string& str);
  /// determines if a char array consists only of whitespace.
  static bool isWhiteLine(const char *p);
  /// removes quotes from strings. 
  static std::string dequote(const std::string& str);

  /// Converts a string to a vector of doubles. 
  static void StringToVec(const std::string str, std::vector<double>&);
  /// Converts a string to a vector of floats. 
  static void StringToVec(const std::string str, std::vector<float>&);
  /// Converts a string to a vector of unsigned ints. 
  static void StringToVec(const std::string str, std::vector<unsigned>&);

  ///determines if a given string is in a list of strings
  static bool isListedString(const std::string s, const std::vector<std::string>list);
  ///converts a string to lower case
  static void toLower(std::string& s);
  ///removes all instances of character c from string s
  static void RemoveCharFromString(std::string& s, const char c);

  /// a dummy function
  static bool dummy();
};
END_BCLIB_NAMESPACE
#endif /* !STRING_CONVERTOR_H */
