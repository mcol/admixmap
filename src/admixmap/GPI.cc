/** 
 *   GPI.cc 
 *   (G)enotype(P)robability(I)terator class
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "GPI.h"
#include "utils/ColumnIter.h"

///default c'tor ,initialises to null values
GenotypeProbIterator::GenotypeProbIterator(){
  p = 0;
  offset = 0;
  g = 0;
}

///c'tor with pointer assignment. Equivalent to default c'tor + setPointer
GenotypeProbIterator::GenotypeProbIterator(const FreqArray* const x, const unsigned cstride){
  setPointer (x, cstride);
}

///sets pointer to GenotypeProbs and initial stride in ColumnIterator member
void GenotypeProbIterator::setPointer(const FreqArray* const x, const unsigned cstride){
  p = x;
  C.setStride(cstride);
}

///assign the GenotypeProb pointer, genotype pointer (may be null), ColumnStride and Offset(number of columns)
void GenotypeProbIterator::assign(const FreqArray* const x, std::vector<unsigned short>* GI,
				  const unsigned cstride, const unsigned t){
  //C.assign( (*p)[offset]);
  setPointer(x, cstride);
  Offset(t);
  g = GI;
}
  
///set offset (number of columns in GenotypeProb array
void GenotypeProbIterator::Offset(const unsigned t){
  offset = t;
  //C.setOffset(Coffset);
}

///Bracket operator, returns a reference to a ColumnIterator object, which can then be used to access individual genotype probs.
// const ColumnIterator& GenotypeProbIterator::operator[](unsigned i){
//   if(!p)
//     throw ("Error in GenotypeProbIterator: unassigned pointer");
//   if(g){//g is pointing to something
//     if((*g)[i+offset] ==0 ){//missing genotype
//       throw ("Error in GenotypeProbIterator: Missing genotype");
//     }
    
//     C.setOffset( (*g)[i+offset] - 1);//-1 because genotypes count from 1
//   }
//   C.assign( (*p)[i+offset]);
  
//   //C.setOffset(g[i]);
//   //C.setOffset(i+offset);
//   return C;
// }

double GenotypeProbIterator::operator()(unsigned i, unsigned j){
  if(!p)
    throw ("Error in GenotypeProbIterator: unassigned pointer");
  if(g){//g is pointing to something
    if((*g)[i+offset] ==0 ){//missing genotype
      return 1.0;
    }
    
    C.setOffset( (*g)[i+offset] - 1);//-1 because genotypes count from 1
  }
  C.assign( (*p)[i+offset]);
  
  //C.setOffset(g[i]);
  //C.setOffset(i+offset);
  return C[j];

}


///Checks if the genotype prob pointer has been assigned
bool GenotypeProbIterator::isNull()const{
  return (bool)(p==0/* || C.isNull()*/);
}

