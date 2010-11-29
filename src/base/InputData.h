//=============================================================================
//
// Copyright (C) 2005-2007  David O'Donnell and Paul McKeigue
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
/// \file InputData.cc
/// Implementationtion of the InputData class.
//=============================================================================

#ifndef INPUT_DATA_H
#define INPUT_DATA_H 1


#include "bclib/DataMatrix.h"
#include "GeneticDistanceUnit.h"
#include "Genotype.h"
#include "SimpleLocusArray.h"
#include "Pedigree.h"
#include <vector>

typedef genepi::GenotypeArray genotype;
using genepi::SimpleLocusArray;

/**
 *  Forward declarations.
 */
class Options;
class Genome;
class Chromosome;

namespace bclib {
  class LogWriter;
}
namespace genepi {
  class GenotypeParser;
}



/** \addtogroup base
 * @{ */


///Class to read and check all input data files
class InputData{

private:

  SimpleLocusArray simpleLoci;

public:

  /**
   *  Constructor.
   */
  InputData();

  /**
   *  Destructor.
   */
  virtual ~InputData();

  ///erases (nearly) all memory used by object
  virtual void Delete();

  /*
   *  Getters to retrieve data (in string form).
   */

  ///return contents of covariatesfile
  const Matrix_s& getCovariatesData()	 const;
  ///return contents of outcomevarfile
  const Matrix_s& getOutcomeVarData()	const;
  ///return contents of priorallelefreqfile
  const Matrix_s& getPriorAlleleFreqData() const { return priorAlleleFreqData_; }

  /*
   *  Getters to retrieve data (converted to DataMatrix).
   */
  const SimpleLocusArray & getSimpleLoci() const { return simpleLoci; }
  const bclib::DataMatrix& getOutcomeVarMatrix() const;
  const bclib::DataMatrix& getCoxOutcomeVarMatrix() const;
  const bclib::DataMatrix& getCovariatesMatrix() const;

  void            getOutcomeTypes      (DataType*) const; ///<
  DataType        getOutcomeType       (unsigned ) const; ///< returns the DataType of an outcome variable
  const Vector_s& GetHiddenStateLabels (         ) const; ///< returns population/hidden state labels
  Vector_s        getOutcomeLabels     (         ) const; ///< returns outcome variable labels
  const Vector_s  getCovariateLabels   (         ) const; ///< returns covariate labels

  const Vector_s & getLocusLabels() const { return LocusLabels; } ///< (from locusfile)


  /**
     Determines if an individual is female.
     Determined from sex column of genotypesfile, if there is one.
     \param i Individual number (count from 1)
     \return true if female, false if male or unknown
  */
  bool isFemale(int i)const;

  /// Returns number of individuals (lines in the genotype/pedigree file)
  size_t getNumberOfIndividuals() const;
  size_t getNumberOfSimpleLoci () const; ///< returns total number of loci in locusfile

  /// Returns number of composite loci
  unsigned getNumberOfCompositeLoci() const
	    { return simpleLoci.getNComposite(); }
  ///returns unit of distance used in locusfile (defaults to Morgans)
  GeneticDistanceUnit getUnitOfDistance()const;
  ///returns a string representing the unit of distance in locusfile
  const char * getUnitOfDistanceAsString() const { return gduAsString(getUnitOfDistance()); }
  float getLocusDistanceThreshold(bool hapmixmodelindicator)const;

protected:
  /** Extracts population labels from header line of allelefreq file.
      \param data header line of file.
      \param Populations number of populations/hidden states
      \param labels vector to fill with labels
  */
  void getPopLabels(const Vector_s& data, size_t Populations, Vector_s& labels);

  void SetLocusLabels();

  /// Determine number of composite loci from locusfile, but mysteriously don't
  /// store it in the corresponding data member
  unsigned determineNumberOfCompositeLoci() const;

  //virtual void CheckAlleleFreqs(Options *options, bclib::LogWriter &Log);
  ///check contents of outcomevarfile
  void CheckOutcomeVarFile(unsigned N, Options * const options, bclib::LogWriter &Log);
  ///check contsnts of coxoutcomevarfile
  void CheckCoxOutcomeVarFile(bclib::LogWriter &log)const;
  ///check contents of covariatesfile
  void CheckCovariatesFile(unsigned NumIndividuals, Options* const options, bclib::LogWriter &Log);
  /**
   *  Read input data and store in internal structures.
   */
  void ReadData(Options *options, bclib::LogWriter &log);

  Matrix_s covariatesData_;
  Matrix_s outcomeVarData_;
  Matrix_s coxOutcomeVarData_;
  private:
    Matrix_s priorAlleleFreqData_;
  protected:

  bclib::DataMatrix covariatesMatrix_;
  bclib::DataMatrix outcomeVarMatrix_;
  bclib::DataMatrix coxOutcomeVarMatrix_;
  private:
    bclib::DataMatrix priorAlleleFreqMatrix_;
  protected:

  genepi::GenotypeParser *genotypeParser;

 private:
  genepi::cvector<genepi::Pedigree> peds;
 protected:
  genepi::cvector<genepi::Pedigree> & getPeds()		{ return peds	   ; }
  genepi::Pedigree &		  getPed( size_t idx )	{ return peds[idx] ; }
 public:
  const genepi::cvector<genepi::Pedigree> & getPeds() const { return peds; }
  const genepi::Pedigree &			getPed( size_t idx ) const { return peds[idx] ; }
  bool isPedFile() const { return genotypeParser->isPedFile(); }

 protected:

  Vector_s HiddenStateLabels;
  std::vector<DataType> OutcomeType;
  Vector_s CovariateLabels;
  Vector_s OutcomeLabels;


private:

  Vector_s LocusLabels;

  //void CheckData(Options *options, bclib::LogWriter &Log);

  /*
   *  UNIMPLEMENTED: to avoid undesired copying.
   */
  ///copy constructor, not implemented
  InputData(const InputData&);
  //assignment operator, not implemented
  void operator=(const InputData&);

protected:
  void generatePedigrees( const Options & options );
};

//void getLabels(const Vector_s& data, std::string *labels);


inline size_t InputData::getNumberOfSimpleLoci() const {
  return getSimpleLoci().size();
}



/** @} */


#endif /* !defined INPUT_DATA_H */
