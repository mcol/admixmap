// *-*-C++-*-*
#ifndef STRING_SPLITTER_H
#define STRING_SPLITTER_H 1

#include "common.h"

/**
 *  @class StringSplitter.
 *  Purpose: used to split string to separate words. Words can be 
 *      enclosed or not enclosed in quotes. Words enclosed in quotes can contain
 *      space symbols.
 *  Important:
 *      Enclosing quotes are removed from resulting strings.
 *  Implementaion: 
 *      Class uses "State" pattern (GoF).
 */
class StringSplitter
{
public:
    /// = Ctr/dtr.

    StringSplitter();
    ~StringSplitter();


    /// = Splitting (see comments for class above).

    const Vector_s& split(const char *str);
    const Vector_s& split(const std::string& str);

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
