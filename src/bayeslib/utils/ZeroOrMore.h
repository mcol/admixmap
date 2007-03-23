#ifndef ZEROORMORE_H_
#define ZEROORMORE_H_

template<class T>
class ZeroOrMore
{
public:
	ZeroOrMore() {};
	virtual ~ZeroOrMore() {};
  inline bool operator () (const T& o) const
  {
    return (o >= 0.0);
  }
};

#endif /*ZEROORMORE_H_*/
