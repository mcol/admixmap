// *-*-C++-*-*

#ifndef GPI_H
#define GPI_H

#include "FreqArrays.h"
#include <vector>

class GPIBase{
 public:

 virtual ~GPIBase(){};
 virtual const ColumnIterator& operator[](unsigned i) = 0;
 unsigned getStride()const{
  return stride;
 }
 virtual bool isNull()const = 0;

 protected:
  ColumnIterator C;
  //unsigned K;
  unsigned stride;
} ;


class HapMixGenotypeProbIterator : public GPIBase{

public:

  HapMixGenotypeProbIterator(){
    p = 0;
    stride = 0;
    offset = 0;
    //K = 0;
  }
  void assign(const FreqArray* const x, const std::vector<unsigned short>::const_iterator geno,
                             const unsigned n = 1, const unsigned t = 0){
    p = x;
    g = geno;
    //K = k;
    stride = n;
    offset = t;
  }

  HapMixGenotypeProbIterator(const FreqArray* const x, const std::vector<unsigned short>::const_iterator geno,
                             const unsigned n = 1, const unsigned t = 0){
    assign (x, geno, n, t);
  };

  // A = B is the same as A.operator=(B)
  void operator=(const HapMixGenotypeProbIterator& rhs){
    stride = rhs.getStride();
    p = rhs.getPointer();
    //K = rhs.getNumStrata();
  }

  const ColumnIterator& operator[](unsigned i){
    if(!p)throw ("pointer error in GenotypeProbIterator");
    C.assign( (*p)[i+offset] + *(g+i+offset)-1, stride);
    //C.setOffset(g[i]);
    return C;
  };

//  unsigned getNumStrata()const{
//   return K;
//  }
  const FreqArray* getPointer()const{
    return p;
  }
  bool isNull()const{
   return (bool)(p==0);
  }

private:
  const DoubleArray* p;
  std::vector<unsigned short>::const_iterator g;
  unsigned offset;


  HapMixGenotypeProbIterator(const HapMixGenotypeProbIterator& );

};

class AdmixGenotypeProbIterator : public GPIBase{

public:

  AdmixGenotypeProbIterator(){
   p = 0;
   stride = 1;
  }

  AdmixGenotypeProbIterator(const double* x, unsigned n){
   assign(x, n);
  };

  void assign(const double* x, unsigned n = 1){
    p = x;
    stride = n;
  }
  // A = B is the same as A.operator=(B)
  void operator=(const AdmixGenotypeProbIterator& rhs){
    p = rhs.getPointer();
    stride = rhs.getStride();
  }

  const ColumnIterator& operator[](unsigned i){
    if(!p)throw ("pointer error in GenotypeProbIterator");
    C.assign( p+i*stride , 1);
    return C;
  };

  const double* getPointer()const{
    return p;
  }
  bool isNull()const{
   return (bool) (p==0);
  }

private:
  const double* p;

  AdmixGenotypeProbIterator(const AdmixGenotypeProbIterator& );

};

class GenotypeProbIterator {

GPIBase* GPI;
 public:

  GenotypeProbIterator(){
    GPI = 0;
  }

  GenotypeProbIterator(GPIBase* p){
    assign(p);
  }

  const ColumnIterator& operator[](unsigned i){
   return GPI->operator[](i);
  }

  void assign(GPIBase* p){
   GPI = p;
   }

  void operator=(const GenotypeProbIterator& rhs){
    this->GPI = rhs.GPI;
  }
  bool isNull()const{
   return (!GPI || GPI->isNull());
  }
} ;

#endif



