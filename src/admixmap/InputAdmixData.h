// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   InputAdmixData.h 
 *   class to read ADMIXMAP data
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef INPUTADMIXDATA_H
#define INPUTADMIXDATA_H

#include "InputData.h"

class AdmixOptions;

///class to read ADMIXMAP data
class InputAdmixData : public InputData{

public:
  /** Constructor
  \param options pointer to AdmixOptions, to access filenames
  \param Log LogWriter, for writing messages
  */
  InputAdmixData(AdmixOptions *options, bclib::LogWriter &Log);

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
      Retrives an Individual's genotypes.
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
  ///default ctor, not implemented    
  InputAdmixData();
  ///copy ctor, not implemented    
  InputAdmixData(const InputAdmixData&);
  ///assignment operator, not implemented
  void operator=(const InputAdmixData&);
};
#endif
