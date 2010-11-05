//=============================================================================
//
// Copyright (C) 2007  David O'Donnell
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file GPI.cc
/// Implementation of the GenotypeProbIterator class.
//=============================================================================

#include "GPI.h"
#include "bclib/ColumnIter.h"

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

