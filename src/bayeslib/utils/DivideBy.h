#ifndef DIVIDEBY_H_
#define DIVIDEBY_H_

/// Provides a function object which performs division
template <class T>
class DivideBy
{
private:
  T divisor;
  T reciprocal;
public:
  DivideBy();
  DivideBy(const T&);
	virtual ~DivideBy();
  inline const T operator () (const T& o) const
  {
    // Multiplication by reciprocal (as opposed to dividing)
    // makes normalize() method run two times faster.
    return o * reciprocal;
  }
  void setDivisor(const T&);
};

#endif /*DIVIDEBY_H_*/
