// *-*-C++-*-*
/** 
 *   GPI.h 
 *   (G)enotype(P)robability(I)terator class
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef GPI_H
#define GPI_H

#include "FreqArrays.h"
#include <vector>

/**
   Class to provide an interface to GenotypeProbs stored in Individual objects.
   Holds a pointer to Genotype probs, stored as a FreqArray object, and a pointer to a vector of genotypes (optional).
   The [] operator is used to access the required probs via a ColumnIterator object. The row number is offset by the number of columns and the column number is offset by the genotype (if provided).
   For example, GPI[locus][state]


   There are 4 main usages:
   (1) 1 GPI per individual: Use default c'tor then setPointer once to set the FreqArray pointer and Offset before passing as a function argument(see AdmixedIndividual)
   (2) 1 GPI shared between individuals: Use default c'tor then assign to effectively reset every time (see HapMixIndividual) 
   (3) As temporary object
   (3)Use assignment operator(=) to copy (see HMM)
*/
class GenotypeProbIterator{

public:

  ///default c'tor ,initialises to null values
  GenotypeProbIterator(){
    p = 0;
    offset = 0;
    g = 0;
  }
  ///c'tor with pointer assignment. Equivalent to default c'tor + setPointer
  GenotypeProbIterator(const FreqArray* const x, const unsigned cstride =1){
    setPointer (x, cstride);
  };

  ///sets pointer to GenotypeProbs and initial stride in ColumnIterator member
  void setPointer(const FreqArray* const x, const unsigned cstride = 1){
    p = x;
    C.setStride(cstride);
  }

  ///assign the GenotypeProb pointer, genotype pointer (may be null), ColumnStride and Offset(number of columns)
  void assign(const FreqArray* const x, std::vector<unsigned short>* GI,
               const unsigned cstride, const unsigned t){
    //C.assign( (*p)[offset]);
    setPointer(x, cstride);
    Offset(t);
    g = GI;
  }
  
  ///set offset (number of columns in GenotypeProb array
  void Offset(const unsigned t){
    offset = t;
    //C.setOffset(Coffset);
  }

  // Assignment operator. A = B is the same as A.operator=(B)
  void operator=(const GenotypeProbIterator& rhs){
    p = rhs.getPointer();
    g = rhs.getGenotypePointer();
    C = rhs.getCI();
    offset = rhs.getOffset();
  }

  ///Bracket operator, returns a reference to a ColumnIterator object, which can then be used to access individual genotype probs.
  const ColumnIterator& operator[](unsigned i){
    if(!p)
      throw ("Error in GenotypeProbIterator: unassigned pointer");
    if(g){//g is pointing to something
      if((*g)[i+offset] ==0 ){//missing genotype
        throw("Error in GenotypeProbIterator: missing genotype");
      }

      C.setOffset( (*g)[i+offset] - 1);//-1 because genotypes count from 1
    }
    C.assign( (*p)[i+offset]);
 
   //C.setOffset(g[i]);
    //C.setOffset(i+offset);
    return C;
  };

  ///Checks if the genotype prob pointer has been assigned
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



