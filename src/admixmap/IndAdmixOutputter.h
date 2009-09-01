// *-*-C++-*-*
/**
 *   ADMIXMAP
 *   IndAdmixOutputter.h
 *   Class to output individual admix parameters
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */
#ifndef IND_ADMIX_OUTPUTTER
#define IND_ADMIX_OUTPUTTER 1

#include "AdmixOptions.h"
#include "Genome.h"
#include "PedBase.h"
#include "bclib/RObjectWriter.h"
#include <vector>
#include <iostream>

class AdmixIndividualCollection;


///Class to output individual admixture proportions and sumintensities to file
class IndAdmixOutputter
{
public:
  IndAdmixOutputter(const AdmixOptions& , const Genome&, const Vector_s& PopLabels);
  virtual ~IndAdmixOutputter();
  void visitIndividual(const genepi::PedBase &, const std::vector<int>);
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

#endif /* !defined IND_ADMIX_OUTPUTTER */
