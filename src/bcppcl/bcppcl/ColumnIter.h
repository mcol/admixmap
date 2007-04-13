// *-*-C++-*-*

#ifndef COLUMN_ITERATOR_H
#define COLUMN_ITERATOR_H

///class to iterate over the columns of a 2D array/matrix
class ColumnIterator{

public:

  ColumnIterator(){
    p = 0;
    stride = 0;
    offset = 0;
  }

  ColumnIterator(const double* const x, unsigned n = 1){
    assign(x, n);
  };
  void assign(const double* const x, unsigned n, unsigned t = 0){
    p = x;
    stride = n;
    offset = t;
  }
  ///assign pointer without changing stride or offset
  void assign(const double* const x){
    p = x;
  }
  void setStride(unsigned n){
    stride = n;
  };
  void setOffset(unsigned t){
    offset = t;
  }

  // A = B is the same as A.operator=(B)
  ColumnIterator& operator=(const ColumnIterator& rhs){
    stride = rhs.getStride();
    p = rhs.getPointer();
    return *this;
  }

  double operator[](unsigned i)const{
    if(!p)throw ("pointer error in ColumnIterator");
    return *(p + offset + i*stride);
  };

  bool isNull()const{
    return (bool)(!p);
  }

  //the remaining functions are to facilitate the assignment operator
  unsigned getStride()const{
    return stride;
  }
  const double* getPointer()const{
    return p;
  }
  unsigned getOffset()const{
    return offset;
  }

private:
  const double* p;
  unsigned stride;
  unsigned offset;

  ColumnIterator(const ColumnIterator& );

};

#endif
