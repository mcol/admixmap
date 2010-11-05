//=============================================================================
//
// Copyright (C) 2006  David O'Donnell
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
/// \file GenotypeProbOutputter.h
/// Definition of the GenotypeProbOutputter class.
//=============================================================================

#ifndef GPOUTPUTTER_H
#define GPOUTPUTTER_H


#include "bclib/pvector.h"
#include "bclib/RObjectWriter.h"
#include <string>
#include <vector>


/** \addtogroup base
 * @{ */


class GenotypeProbOutputter{
public:

  void Initialise(unsigned Nindivs, unsigned Nloci);
  void Update(unsigned i, unsigned j, const std::vector<std::vector<double> >&);
  void Output(const char* filename, const std::vector<std::string>& LocusLabels);

private:
  unsigned NumMaskedIndivs;
  unsigned NumMaskedLoci;
  unsigned NumIterations;
  bclib::pvector<double> Probs;
  //std::ofstream outfile;
  bclib::RObjectWriter outfile;

  /** Three-dimensional vector for probabilities.
   * Dimensions are:
   * 1. Individual
   * 2. Locus
   * 3. Genotype
   */
  std::vector<std::vector<std::vector<double> > > sumProbs;
};




/** @} */


#endif
