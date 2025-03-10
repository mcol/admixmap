//=============================================================================
//
// Copyright (C) 2007  David O'Donnell and Paul McKeigue
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
/// \file InputAdmixData.h
/// Definition of the InputAdmixData class.
//=============================================================================

#ifndef INPUTADMIXDATA_H
#define INPUTADMIXDATA_H

#include "InputData.h"

class AdmixOptions;


/** \addtogroup admixmap
 * @{ */


///class to read ADMIXMAP data
class InputAdmixData : public InputData{

public:
  /** Constructor
  \param options pointer to AdmixOptions, to access filenames
  \param Log LogWriter, for writing messages
  */
  InputAdmixData( AdmixOptions & options, bclib::LogWriter & Log );

  /// The object is not constructed and cannot be used until this is called, but
  /// it cannot be done from the constructor (see note in source file for why).
  /// Ugly and brittle, this is due to the circular dependency between
  /// AdmixOptions and InputAdmixData.
  void finishConstructing( const AdmixOptions & options );

  ///Destructor
  ~InputAdmixData(){};

  void Delete();
  /*
   *  Getters to retrieve data (in string form).
   */
  ///returns contents of allelefreqfile
  const Matrix_s& getAlleleFreqData() const;
  ///returns contents of historicallelefreqfile
  const Matrix_s& getHistoricalAlleleFreqData() const;
  ///returns contents of etapriorfile
  const Matrix_s& getEtaPriorData() const;
  ///returns contents of reportedancesryfile
  const Matrix_s& getReportedAncestryData() const;

  /*
   *  Getters to retrieve data (converted to DataMatrix).
   */    
  //const bclib::DataMatrix& getAlleleFreqMatrix() const;
  //const bclib::DataMatrix& getHistoricalAlleleFreqMatrix() const;
  //const bclib::DataMatrix& getPriorAlleleFreqMatrix() const;
  const bclib::DataMatrix& getEtaPriorMatrix() const;
  const bclib::DataMatrix& getReportedAncestryMatrix() const;

  /** 
      Retrieves an Individual's genotypes.
      \param i the Individual number
      \param Loci Genome object, to read locus sizes
      \param genotypes pointer to Individual's genotypes
      \param Missing Individual's array of missing genotype indicators
      \see GenotypeLoader::GetGenotype
  */
  void GetGenotype(int i, const Genome &Loci, std::vector<genotype>* genotypes, bool **Missing)const;

  ///get population labels
  const Vector_s& GetPopLabels() const;
private:
  Matrix_s alleleFreqData_;
  Matrix_s historicalAlleleFreqData_;
  Matrix_s etaPriorData_;
  Matrix_s reportedAncestryData_;

  bclib::DataMatrix alleleFreqMatrix_;
  bclib::DataMatrix historicalAlleleFreqMatrix_;
  bclib::DataMatrix etaPriorMatrix_;
  bclib::DataMatrix reportedAncestryMatrix_;

  /**
     Check data after it has been read in.
     \param options pointer to options
     \param Log LogWriter, for messages
  */
  void CheckData(AdmixOptions *options, bclib::LogWriter &Log);

  /**
     Reads population labels.
     Reads from header of (prior/historic)allelefreqfile, if specified.
     Otherwise, assigns default names - Pop1, Pop2 etc.
  */
  void ReadPopulationLabels(AdmixOptions *options);
  /**
     checks contents of (prior/historic)allelefreqfile.
     checks consistency of supplied allelefreqs with locusfile
     and determines number of populations and population labels.
     \param options pointer to options object
     \param Log LogWriter, for error messages
  */
  void CheckAlleleFreqs(AdmixOptions *options, bclib::LogWriter &Log);
  /**
     Check contents of reportedancestryfile.
     Checks number of rows is the same as in genotypesfile and number of
     columns is the same as the number of populations.
     \param populations number of populations in model
     \param Log LogWriter, for error messages
  */
  void CheckRepAncestryFile(int populations, bclib::LogWriter &Log)const;
  /*
   *  UNIMPLEMENTED: to avoid undesired copying.
   */
  /// Copy ctor, not implemented.
  InputAdmixData(const InputAdmixData&);
  /// Assignment operator, not implemented.
  void operator=(const InputAdmixData&);
};


/** @} */


#endif
