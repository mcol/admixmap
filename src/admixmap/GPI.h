// *-*-C++-*-*

#ifndef GPI_H
#define GPI_H

#include "FreqArrays.h"
#include <vector>

class GenotypeProbIterator{

public:

  GenotypeProbIterator(){
    p = 0;
    offset = 0;
    g = 0;
  }
  GenotypeProbIterator(const FreqArray* const x, const unsigned cstride =1){
    setPointer (x, cstride);
  };

  void setPointer(const FreqArray* const x, const unsigned cstride = 1){
    p = x;
    C.setStride(cstride);
  }

  void assign(const FreqArray* const x, std::vector<unsigned short>* GI,
               const unsigned cstride, const unsigned t){
    //C.assign( (*p)[offset]);
    setPointer(x, cstride);
    Offset(t);
    g = GI;
  }
  
  void Offset(const unsigned t){
    offset = t;
    //C.setOffset(Coffset);
  }

  // A = B is the same as A.operator=(B)
  void operator=(const GenotypeProbIterator& rhs){
    p = rhs.getPointer();
    g = rhs.getGenotypePointer();
    C = rhs.getCI();
    offset = rhs.getOffset();
  }

  const ColumnIterator& operator[](unsigned i){
    if(!p)
      throw ("pointer error in GenotypeProbIterator");

    if(g)//g is pointing to something
      C.setOffset( (*g)[i+offset] - 1);//-1 because genotypes count from 1
    C.assign( (*p)[i+offset]);
 
   //C.setOffset(g[i]);
    //C.setOffset(i+offset);
    return C;
  };

  bool isNull()const{
    return (bool)(p==0/* || C.isNull()*/);
  }

  //the remaining functions are to facilitate the assignment operator
  const FreqArray* getPointer()const{
    return p;
  }

  unsigned getOffset()const{
    return offset;
  }
  const ColumnIterator& getCI()const{
    return C;
  };
  const std::vector<unsigned short>* getGenotypePointer()const{
    return g;
  }

private:
  ColumnIterator C;
  const FreqArray* p;
  const std::vector<unsigned short>* g;//genotype pointer to indicate columns of the Probs array to use
                                //using pointer becasue iterators have no null state
  unsigned offset;

  GenotypeProbIterator(const GenotypeProbIterator& );

};

#endif



