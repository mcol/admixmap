// *-*-C++-*-*
/**
 *   InputData.h
 *   header file for InputData class
 */

/*
 *   Copyright (c) 2005 - 2007 David O'Donnell and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */
#ifndef INPUT_DATA_H
#define INPUT_DATA_H 1


#include "common.h"
#include "bclib/DataMatrix.h"
#include "GeneticDistanceUnit.h"
#include "config.h" // USE_GENOTYPE_PARSER

#if USE_GENOTYPE_PARSER

    namespace genepi { class GenotypeParser; }
    #define GenotypeLoader genepi::GenotypeParser

    #include "Genotype.h"
    typedef genepi::GenotypeArray genotype;

    #include "SimpleLocusArray.h"
    using genepi::SimpleLocusArray;

    #include "Pedigree.h"

#else
    #include "GenotypeLoader.h"
    #define GenotypeArray std::vector<genotype>
#endif


#include <vector>

/**
 *  Forward declarations.
 */
class Options;
class Genome;
class Chromosome;


namespace bclib{
  class LogWriter;
}



/** \addtogroup base
 * @{ */


///Class to read and check all input data files
class InputData{

private:

  #if USE_GENOTYPE_PARSER
    SimpleLocusArray simpleLoci;
  #endif


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
  #if ! USE_GENOTYPE_PARSER
      ///return contents of locusfile
      const Matrix_s& getLocusData() const;
  #endif
  ///return contents of covariatesfile
  const Matrix_s& getCovariatesData()	 const;
  ///return contents of outcomevarfile
  const Matrix_s& getOutcomeVarData()	const;
  ///return contents of priorallelefreqfile
  const Matrix_s& getPriorAlleleFreqData() const { return priorAlleleFreqData_; }

  /*
   *  Getters to retrieve data (converted to DataMatrix).
   */
  #if USE_GENOTYPE_PARSER
    const SimpleLocusArray & getSimpleLoci() const { return simpleLoci; }
  #else
    const bclib::DataMatrix& getLocusMatrix() const;
  #endif
  const bclib::DataMatrix& getOutcomeVarMatrix() const;
  const bclib::DataMatrix& getCoxOutcomeVarMatrix() const;
  const bclib::DataMatrix& getCovariatesMatrix() const;

  void            getOutcomeTypes      (DataType*) const; ///<
  DataType        getOutcomeType       (unsigned ) const; ///< returns the DataType of an outcome variable
  const Vector_s& GetHiddenStateLabels (         ) const; ///< returns population/hidden state labels
  Vector_s        getOutcomeLabels     (         ) const; ///< returns outcome variable labels
  const Vector_s  getCovariateLabels   (         ) const; ///< returns covariete labels

  #if 1 || ! USE_GENOTYPE_PARSER
    const Vector_s & getLocusLabels() const { return LocusLabels; } ///< (from locusfile)
  #endif


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
	#if USE_GENOTYPE_PARSER
	    { return simpleLoci.getNComposite(); }
	#else
	    { return NumCompositeLoci; }
	#endif
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

  /**
     Checks contents of locusfile and extracts locus labels.
     If check = true, checks that loci labels in locusfile are unique and that they match the names in the genotypes file.
     \param sexColumn the index of the sex column in genotypesfile if there is one; 0 is not.
     \param threshold maximum allowable distance between loci on a chromosome. If a distance exceeds this, the chromosome is broken up and a warning is printed.
     \param check determines whether to check locus labels
  */
  #if ! USE_GENOTYPE_PARSER
    bool checkLocusFile(Options *options, bclib::LogWriter& Log);
  #endif
  void SetLocusLabels();

  /// Determine number of composite loci from locusfile, but mysteriously don't
  /// store it in the corresponding data member
  unsigned determineNumberOfCompositeLoci() const;
  #if ! USE_GENOTYPE_PARSER
    GeneticDistanceUnit DetermineUnitOfDistance(); ///< determine unit of distance in locusfile
  #endif

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

  #if ! USE_GENOTYPE_PARSER
    Matrix_s locusData_;
  #endif
  Matrix_s covariatesData_;
  Matrix_s outcomeVarData_;
  Matrix_s coxOutcomeVarData_;
  private:
    Matrix_s priorAlleleFreqData_;
  protected:

  #if ! USE_GENOTYPE_PARSER
    bclib::DataMatrix locusMatrix_;
  #endif
  bclib::DataMatrix covariatesMatrix_;
  bclib::DataMatrix outcomeVarMatrix_;
  bclib::DataMatrix coxOutcomeVarMatrix_;
  private:
    bclib::DataMatrix priorAlleleFreqMatrix_;
  protected:

  GenotypeLoader * genotypeLoader;
  #if USE_GENOTYPE_PARSER
    private:
      genepi::cvector<genepi::Pedigree> peds;
    protected:
      genepi::cvector<genepi::Pedigree> & getPeds()		{ return peds	   ; }
      genepi::Pedigree &		  getPed( size_t idx )	{ return peds[idx] ; }
    public:
      const genepi::cvector<genepi::Pedigree> & getPeds() const { return peds; }
  const genepi::Pedigree &			getPed( size_t idx ) const { return peds[idx] ; }
      bool isPedFile() const { return genotypeLoader->isPedFile(); }
    protected:
  #else
    size_t		NumSimpleLoci	;
    GeneticDistanceUnit distanceUnit	;
    unsigned		NumCompositeLoci;
  #endif
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


#if USE_GENOTYPE_PARSER
    inline size_t InputData::getNumberOfSimpleLoci() const { return getSimpleLoci().size(); }
#else
    inline size_t InputData::getNumberOfSimpleLoci() const { return (locusData_.size() - 1); }
#endif




/** @} */


#endif /* !defined INPUT_DATA_H */
