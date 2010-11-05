//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AncestryAssocTest.h
/// Definition of the AncestryAssocTest class.
//=============================================================================

#ifndef ANCESTRYASSOCTEST_H
#define ANCESTRYASSOCTEST_H

#include "CopyNumberAssocTest.h"
#include "common.h"   // for Vector_s

class Genome;
namespace bclib{
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


class AncestryAssocTest : public CopyNumberAssocTest{

public:
  AncestryAssocTest();
  ~AncestryAssocTest();
  void Initialise(const char* filename, const int NumPopulations, const int NumLoci);
  void Output(const Vector_s& PopLabels, const Genome& Loci);
  void WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
		       const Genome& Loci, bclib::LogWriter& Log);
private:
  unsigned firstpoplabel;

};


/** @} */


#endif
