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
/// \file GPI.h
/// Definition of the GenotypeProbIterator class.
//=============================================================================

#ifndef GPI_H
#define GPI_H

#include "FreqArrays.h"
#include "bclib/ColumnIter.h"
#include <vector>


/** \addtogroup base
 * @{ */


/**
   \brief Class to provide an interface to GenotypeProbs stored in Individual objects.

   Holds a pointer to Genotype probs, stored as a FreqArray object, and a pointer to a vector of genotypes (optional).
   The () operator is used to access the required probs via a ColumnIterator object. 
   The row number is offset by the number of columns and the column number is offset by the genotype (if provided).
   For example, GPI(locus, state). 
   Returns 1.0 if the genotype is being used and is missing.

   There are 4 main usages:
   (1) 1 GPI per individual: Use default c'tor then setPointer once to set the FreqArray pointer and Offset before passing as a function argument(see AdmixedIndividual)
   (2) 1 GPI shared between individuals: Use default c'tor then assign to effectively reset every time (see HapMixIndividual) 
   (3) As temporary object
   (3)Use assignment operator(=) to copy (see HMM)
*/
class GenotypeProbIterator{

public:

  ///default c'tor ,initialises to null values
  GenotypeProbIterator();

  ///c'tor with pointer assignment. Equivalent to default c'tor + setPointer
  GenotypeProbIterator(const FreqArray* const x, const unsigned cstride = 1);

  ///sets pointer to GenotypeProbs and initial stride in ColumnIterator member
  void setPointer(const FreqArray* const x, const unsigned cstride = 1);

  ///assign the GenotypeProb pointer, genotype pointer (may be null), ColumnStride and Offset(number of columns)
  void assign(const FreqArray* const x, std::vector<unsigned short>* GI,
	      const unsigned cstride, const unsigned t);
  
  ///set offset (number of columns in GenotypeProb array
  void Offset(const unsigned t);

  ///Bracket operator, returns a reference to a ColumnIterator object, which can then be used to access individual genotype probs.
  //const ColumnIterator& operator[](unsigned i);

  ///Alternative bracket operator, returns genotype prob for given row, column
  double operator()(unsigned i, unsigned j);

  ///Checks if the genotype prob pointer has been assigned
  bool isNull()const;

private:
  ColumnIterator C;
  const FreqArray* p;
  const std::vector<unsigned short>* g; ///genotype pointer to indicate columns of the Probs array to use
					///using pointer becasue iterators have no null state
  unsigned offset;

  GenotypeProbIterator(const GenotypeProbIterator& );

};



/** @} */



#endif
