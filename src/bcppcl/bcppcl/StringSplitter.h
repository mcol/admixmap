// *-*-C++-*-*
#ifndef STRING_SPLITTER_H
#define STRING_SPLITTER_H 1

#include <vector>
#include <string>

typedef std::vector<std::string> Vector_s; //std vector of strings 
typedef std::vector<Vector_s>    Matrix_s; // std vector of std vectors

/**
 *  @class StringSplitter.
 *  Used to split string to separate words. Words can be 
 *      enclosed or not enclosed in quotes. Words enclosed in quotes can contain
 *      space symbols.
 *  Important:
 *      Enclosing quotes are removed from resulting strings.
 *  Implementation: 
 *      Class uses "State" pattern (GoF).
 */
class StringSplitter
{
public:
  StringSplitter();
  ~StringSplitter();
  
  /// = Splitting (see comments for class above).

  /// splits a char array into a vector of strings, splitting on whitespace ro supplied delimiter.  
  const Vector_s& split(const char *str, char delim = ' ');
  /// splits a string into a vector of strings, splitting on whitespace or supplied delmiter. 
  const Vector_s& split(const std::string& str, char delim = ' ');
  ///   Splits a string into tokens, splitting on any character(s).
  static void Tokenize(const std::string& str,
		       std::vector<std::string>& tokens,
		       const std::string& delimiters);
  static void QuickTokenize(const std::string& text, std::vector<std::string>& tokens, const std::string& delim);

private:
    struct StringSplitterState;
    struct StateWhite;
    struct StateInWord;
    struct StateInString;

    friend struct StateWhite;
    friend struct StateInWord;
    friend struct StateInString;

    static StateWhite    STATE_WHITE;
    static StateInWord   STATE_IN_WORD;
    static StateInString STATE_IN_STRING;

    void addChar(char c);
    void wordComplete();

    std::string current_;
    Vector_s result_;

private:
    // To avoid undesired copying.
    StringSplitter(const StringSplitter&);
    void operator=(const StringSplitter&);
};

#endif /* STRING_SPLITTER_H */
