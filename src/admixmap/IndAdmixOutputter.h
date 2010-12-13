//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file IndAdmixOutputter.h
/// Definition of the IndAdmixOutputter class.
//=============================================================================

#ifndef IND_ADMIX_OUTPUTTER
#define IND_ADMIX_OUTPUTTER 1

#include "AdmixOptions.h"
#include "Genome.h"
#include "PedBase.h"
#include "bclib/RObjectWriter.h"
#include <vector>

class AdmixIndividualCollection;


/** \addtogroup admixmap
 * @{ */


///Class to output individual admixture proportions and sumintensities to file
class IndAdmixOutputter
{
public:
  IndAdmixOutputter(const AdmixOptions& , const Genome&, const Vector_s& PopLabels);
  virtual ~IndAdmixOutputter();
  void visitIndividual(const genepi::PedBase&, const std::vector<int>&);
  void visitIndividualCollection(const AdmixIndividualCollection&);


private:
  bclib::RObjectWriter _out;

  const AdmixOptions&  _options;
  const Genome& _Loci;
  const Vector_s&  _PopulationLabels;
  const bool _RandomMatingModelIndicator;

  int _iterations;
  int _totalIndividuals;
  int _currentIndividual;

  // UNIMPLEMENTED
  // to avoid use
  IndAdmixOutputter();
  IndAdmixOutputter(const IndAdmixOutputter&);
  IndAdmixOutputter& operator=(const IndAdmixOutputter&);

};


/** @} */


#endif /* !defined IND_ADMIX_OUTPUTTER */
