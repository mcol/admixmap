#ifndef IFREQARRAY_H_
#define IFREQARRAY_H_

class IFreqArray
{
public:
	// IFreqArray();
	virtual ~IFreqArray();
  virtual void dealloc(int) = 0;
  virtual double* operator[](unsigned) = 0;
  virtual const double* operator[](unsigned i)const = 0;
};

#endif /*IFREQARRAY_H_*/
