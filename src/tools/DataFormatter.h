#ifndef DATAFORMATTER_H_
#define DATAFORMATTER_H_

#include <ostream>
#include <string>

class UserData;
class Legend;

using std::ostream;
using std::string;

class DataFormatter
{
public:
  DataFormatter();
  DataFormatter(UserData *, Legend *);
  virtual ~DataFormatter();
  virtual int getNoFiles() = 0;
  virtual void renderFile(ostream&, int) = 0;
protected:
  UserData *vpUserData;
  Legend *vpLegend;
};

#endif /*DATAFORMATTER_H_*/
