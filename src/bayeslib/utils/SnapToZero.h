#ifndef SNAPTOZERO_H_
#define SNAPTOZERO_H_

template <class T>
class SnapToZero
{
private:
  T threshold;
public:
  SnapToZero() {}
  SnapToZero(const T& t) : threshold(t) {}
	virtual ~SnapToZero() {};
  void setThreshold(const T& t) {
    threshold = t;
  }
  inline T operator () (const T& v) const
  {
    if (v < threshold) { return 0; }
    else               { return v; }
  }
};

#endif /*SNAPTOZERO_H_*/
