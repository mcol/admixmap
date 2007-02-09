// *-*-C++-*-*

#ifndef COLUMN_ITERATOR_H
#define COLUMN_ITERATOR_H

///class to iterate over the columns of a 2D array/matrix
class ColumnIterator{

public:

  ColumnIterator(){
    p = 0;
    stride = 0;
  }

  ColumnIterator(const double* const x, unsigned n = 1){
    p = x;
    stride = n;
  };

  // A = B is the same as A.operator=(B)
  ColumnIterator& operator=(const ColumnIterator& rhs){
    stride = rhs.getStride();
    p = rhs.getPointer();
    return *this;
  }

  void setStride(unsigned n){
    stride = n;
  };

  double operator[](unsigned i)const{
    if(!p)throw ("pointer error in ColumnIterator");
    return *(p + i*stride);
  };

  unsigned getStride()const{
    return stride;
  }
  const double* getPointer()const{
    return p;
  }

private:
  const double* p;
  unsigned stride;

};

#endif
