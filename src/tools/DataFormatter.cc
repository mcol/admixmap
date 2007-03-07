#include "DataFormatter.h"

/**
 * An abstract base class for data formatters. 
 * 
 * To write a new formatter, create a class with DataFormatter as a base class,
 * Implement int getNoFiles() which should return the number
 * of files to be written and then implement void renderFile(ostream&, int)
 * which should write file with given number, to the given output stream.
 */

DataFormatter::DataFormatter()
{
}

DataFormatter::DataFormatter(UserData *pUserData, Legend *pLegend)
  : vpUserData(pUserData)
  , vpLegend(pLegend)
{
}

DataFormatter::~DataFormatter()
{
}
