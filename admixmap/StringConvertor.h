// *-*-C++-*-*
#ifndef STRING_CONVERTOR_H
#define STRING_CONVERTOR_H 1

#include "common.h"

/**
 *  Utility class used to convert string to other types.
 *  String can be enclosed to doouble quotes.
 */
class StringConvertor
{
public:
    static int toInt(const std::string& str);

    static double toFloat(const std::string& str);

    static std::pair<int, int> toIntPair(const std::string& str);

    static bool isMissingValue(const std::string& str);

private:
    static std::string dequote(const std::string& str);
};

#endif /* !STRING_CONVERTOR_H */
