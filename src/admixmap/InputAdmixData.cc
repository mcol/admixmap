/** 
 *   ADMIXMAP
 *   InputAdmixData.cc 
 *   class to read ADMIXMAP data
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *   Copyright (C) 2009  David D. Favro  gpl@meta-dynamic.com
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "InputAdmixData.h"

#include "bclib/DataReader.h"
#include "config.h" // USE_GENOTYPE_PARSER
#include "estr.h"

#include <cstring>  // strlen()
#include <typeinfo> // typeid


#include "AdmixOptions.h"

#if USE_GENOTYPE_PARSER
    #include "AdmixmapGenotypeConverter.h"
#endif

using bclib::LogWriter;

InputAdmixData::InputAdmixData(AdmixOptions *options, LogWriter &Log){
  using bclib::DataReader;

  #if ! USE_GENOTYPE_PARSER
    genotypeLoader = new GenotypeLoader;
  #endif

  Log.setDisplayMode(bclib::Quiet);
  // Read all input files.
  try {
   //read base data files
    ReadData(options, Log);

    DataReader::ReadData(options->getAlleleFreqFilename(), alleleFreqData_, Log);
    DataReader::ReadData(options->getHistoricalAlleleFreqFilename(), historicalAlleleFreqData_, Log);            
    DataReader::ReadData(options->getEtaPriorFilename(), etaPriorData_,etaPriorMatrix_,  Log, false);//no header
    DataReader::ReadData(options->getReportedAncestryFilename(), reportedAncestryData_, reportedAncestryMatrix_, Log);
    
    Log << "\n";
  }

#if USE_GENOTYPE_PARSER
  // Catch data-validation errors and re-throw, just so that they are not caught
  // by the other handlers here:
  catch ( genepi::DataValidError & e )
    {
    #if 0
	cerr << options->getProgramName() << ": input data validation error: " << e.what() << endl;
	exit(1);
    #else
	throw;
    #endif
    }
#endif

  catch (const exception& e) {
    cerr << "\nException (" << typeid(e).name() << ") occured during parsing of input file:\n" << e.what() << endl;
    exit(1);
  }
  catch(string s){
    cerr << "\nException (strin) occured during parsing of input file:\n" << s << endl;;
    exit(1);
  }
 
  CheckData(options, Log);
}



//-----------------------------------------------------------------------------
// CheckData()
//-----------------------------------------------------------------------------

void InputAdmixData::CheckData(AdmixOptions *options, LogWriter &Log){

#if ! USE_GENOTYPE_PARSER
  NumSimpleLoci = getNumberOfSimpleLoci();
  distanceUnit = DetermineUnitOfDistance();
  NumCompositeLoci = determineNumberOfCompositeLoci();
#endif

  Log.setDisplayMode(bclib::Quiet);
 
  bool badData = false;
  if(options->CheckData())
    badData = !checkLocusFile(options, Log);

  if(badData)
    exit(1);

  SetLocusLabels();

  CheckAlleleFreqs(options, Log);
  ReadPopulationLabels(options);

  //detects regression model
  if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
    const unsigned N = (genotypeLoader->getNumberOfIndividuals() - options->getTestOneIndivIndicator());
    if ( strlen( options->getOutcomeVarFilename() ) != 0 )
      CheckOutcomeVarFile( N, options, Log);
    if ( strlen( options->getCoxOutcomeVarFilename() ) != 0 ){
      OutcomeType.push_back( CoxData );
      if(options->CheckData())
	CheckCoxOutcomeVarFile( Log);
    }
    if ( strlen( options->getCovariatesFilename() ) != 0 )
      CheckCovariatesFile(N, options, Log);
    //append population labels to covariate labels
    if(!options->getTestForAdmixtureAssociation()){
      for( vector<string>::const_iterator i = HiddenStateLabels.begin()+1; i !=HiddenStateLabels.end(); ++i ){
	cout << "HiddenStateLabels " << *i << endl;
	CovariateLabels.push_back("slope." + *i); 
      }
    }
  }
  
  if ( strlen( options->getReportedAncestryFilename() ) != 0 )
    CheckRepAncestryFile(options->getPopulations(), Log);

}



//-----------------------------------------------------------------------------
// GetGenotype()
//-----------------------------------------------------------------------------

void InputAdmixData::GetGenotype(int i, const Genome &Loci, 
			   std::vector<genotype>* genotypes, bool **Missing) const {
  #if USE_GENOTYPE_PARSER
    convert( (*genotypeLoader)[i-1], Loci, *genotypes, Missing );
  #else
    genotypeLoader->GetGenotype(i, Loci, genotypes, Missing);
  #endif

  #if 0 // ********* DEBUG *********

      fprintf( stderr, "\n===== GENOTYPES[%d] (size: %lu) =====\n\n", i, genotypes->size() );
      for ( size_t idx1 = 0; idx1 < genotypes->size(); ++idx1 )
	{
	const genotype & gt = (*genotypes)[ idx1 ];
	fprintf( stderr, "  %lu[%lu]: ", idx1, gt.size() );
	for ( size_t idx2 = 0 ; idx2 < gt.size(); ++idx2 )
	    #if USE_GENOTYPE_PARSER
		fprintf( stderr, " %s", gt[idx2].desc().c_str() );
	    #else
		{
		const vector<unsigned short> & x = gt[idx2];
		if ( x.size() == 0 )
		    fprintf( stderr, " |" );
		else if ( x.size() == 1 )
		    fprintf( stderr, " %hu", x[0] );
		else if ( x.size() == 2 )
		    {
		    if ( x[0] == 0 )
			fprintf( stderr, " " );
		    else
			fprintf( stderr, " %hu", x[0] );
		    if ( x[1] == 0 )
			fprintf( stderr, "," );
		    else
			fprintf( stderr, ",%hu", x[1] );
		    }
		else
		    fprintf( stderr, " WTF[%lu]", x.size() );
		}
	    #endif
	fprintf( stderr, "\n" );
	}

      fprintf( stderr, "\n===== MISSING[%d] =====\n\n", i );
      for ( size_t chrm = 0; chrm < Loci.GetNumberOfChromosomes(); ++chrm )
	{
	fprintf( stderr, "  %lu:", chrm );
	for ( size_t loc = 0 ; loc < Loci.GetSizeOfChromosome(chrm) ; ++loc )
	    fprintf( stderr, " %s", Missing[chrm][loc] ? "M" : "-" );
	fprintf( stderr, "\n" );
	}

  #endif // ********* DEBUG *********
}



//-----------------------------------------------------------------------------
// ReadPopulationLabels()
//-----------------------------------------------------------------------------

void InputAdmixData::ReadPopulationLabels(AdmixOptions *options){
  using bclib::DataReader;
  //  if(strlen(options->getAlleleFreqFilename()) || strlen(options->getPriorAlleleFreqFilename()) || strlen(options->getHistoricalAlleleFreqFilename())){
  if(strlen(options->getAlleleFreqFilename()))
    DataReader::ReadHeader(options->getAlleleFreqFilename(), HiddenStateLabels);
  else if ( strlen(options->getPriorAlleleFreqFilename()) != 0 )
    DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), HiddenStateLabels);
  else if(strlen(options->getHistoricalAlleleFreqFilename()))
    DataReader::ReadHeader(options->getHistoricalAlleleFreqFilename(), HiddenStateLabels);
  else{
    //read labels from 'poplabels' option
    bclib::StringSplitter::Tokenize(options->getPopLabelString(), HiddenStateLabels, " ,");
    if(HiddenStateLabels.size() != (unsigned)options->getPopulations()){
      HiddenStateLabels.clear();

      //set default pop labels
      for( int j = 0; j < options->getPopulations(); j++ )
	HiddenStateLabels.push_back( genepi::estr("Pop") + (j+1) );
    }
  }
}
 

void InputAdmixData::CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log){
  string freqtype = "";
  bool infile = false;//indicates whether either of the three allelefreq files are specified
  int nrows=0, expectednrows=0;
  int Populations = options->getPopulations();
  int NumberOfStates = 0;

  #if USE_GENOTYPE_PARSER

    int nStatesThisCompLoc = 0;
    for ( SimpleLocusArray::ConstIter slIter = getSimpleLoci().begin() ;
			slIter != getSimpleLoci().end() ; ++slIter )
	if ( slIter->isCompositeWithPrevious() )
	    nStatesThisCompLoc *= slIter->getNumAlleles();
	else
	    {
	    NumberOfStates += nStatesThisCompLoc;
	    nStatesThisCompLoc = slIter->getNumAlleles();
	    }
    NumberOfStates += nStatesThisCompLoc;

  #else // ! USE_GENOTYPE_PARSER:

    unsigned index = 0;
    for ( unsigned i = 0; i < NumCompositeLoci; ++i ) {
      int states = 1;
      do{
	states *= (int)locusMatrix_.get( index, 0 );
	index++;
	}
	while( index < locusMatrix_.nRows() && !locusMatrix_.isMissing(index, 1) && locusMatrix_.get( index, 1 ) == 0 );
      NumberOfStates += states;
    }

  #endif // ! USE_GENOTYPE_PARSER


  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() )){
    freqtype = "";
    infile = true;
    nrows = alleleFreqData_.size()-1;
    expectednrows = NumberOfStates-getNumberOfCompositeLoci();
    Populations = alleleFreqData_[0].size() - 1;// -1 for ids in first col
    //getPopLabels(alleleFreqData_[0], Populations, PopulationLabels);
  }
  
  //Historic allelefreqs
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    freqtype = "historic";
    infile = true;
    nrows = historicalAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = historicalAlleleFreqData_[0].size() - 1;
    //getPopLabels(historicalAlleleFreqData_[0], Populations, PopulationLabels);

  }
  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
    freqtype = "prior";
    infile = true;
    nrows = getPriorAlleleFreqData().size();
    expectednrows = NumberOfStates+1;
    Populations = getPriorAlleleFreqData()[0].size() - 1;
    //getPopLabels(getPriorAlleleFreqData()[0], Populations, PopulationLabels);
  }
  if(infile){
    if(nrows != expectednrows){
      Log << "Incorrect number of rows in " << freqtype << "allelefreqfile.\n" 
	  << "Expecting " << expectednrows << " rows, but there are " << nrows << " rows.\n";
      exit(1);
    }
    options->setPopulations(Populations);
  }
  else{//'populations' option
    if(Populations < 1){
      Log << "ERROR: populations = " << options->getPopulations() << "\n";
      exit(1);
    }

    //     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    //       if(Loci->GetNumberOfStates(i) < 2){
    // 	Log << "ERROR: The number of alleles at a locus is " << Loci->GetNumberOfStates(i) << "\n";
    // 	exit(1);
    //       }
    //     }
  }
}

void InputAdmixData::CheckRepAncestryFile(int populations, LogWriter &Log)const{
  if( reportedAncestryMatrix_.nRows() != 2 * genotypeLoader->getNumberOfIndividuals() ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nRows() << " rows\n"
	    "Genotypesfile has " << genotypeLoader->getNumberOfIndividuals() << " rows\n";
    exit(1);}
  if( (int)reportedAncestryMatrix_.nCols() != populations ){
    Log << "ERROR: " << "ReportedAncestry file has " << reportedAncestryMatrix_.nCols() << " cols\n"
	<< "AlleleFreq file has "<< populations << " cols\n";
    exit(1);
  }
}

const Matrix_s& InputAdmixData::getAlleleFreqData() const{
  return alleleFreqData_;
}

const Matrix_s& InputAdmixData::getHistoricalAlleleFreqData() const{
  return historicalAlleleFreqData_;
}

const Matrix_s& InputAdmixData::getEtaPriorData() const{
  return etaPriorData_;
}

const Matrix_s& InputAdmixData::getReportedAncestryData() const{
  return reportedAncestryData_;
}

const bclib::DataMatrix& InputAdmixData::getEtaPriorMatrix() const{
  return etaPriorMatrix_;
}

// const bclib::DataMatrix& InputAdmixData::getAlleleFreqMatrix() const {
//     return alleleFreqMatrix_;
// }

// const bclib::DataMatrix& InputAdmixData::getHistoricalAlleleFreqMatrix() const {
//     return historicalAlleleFreqMatrix_;
// }

const bclib::DataMatrix& InputAdmixData::getReportedAncestryMatrix() const{
  return reportedAncestryMatrix_;
}
const Vector_s& InputAdmixData::GetPopLabels() const{
  return HiddenStateLabels;
}

void InputAdmixData::Delete(){
  InputData::Delete();

  //erase string matrices
  for(unsigned i = 0; i < alleleFreqData_.size(); ++i)
    alleleFreqData_[i].clear();
  alleleFreqData_.clear();
   for(unsigned i = 0; i < historicalAlleleFreqData_.size(); ++i)
    historicalAlleleFreqData_[i].clear();
  historicalAlleleFreqData_.clear();
  for(unsigned i = 0; i < etaPriorData_.size(); ++i)
    etaPriorData_[i].clear();
  etaPriorData_.clear();
  for(unsigned i = 0; i < reportedAncestryData_.size(); ++i)
    reportedAncestryData_[i].clear();
  reportedAncestryData_.clear();

  //erase data matrices 
  //alleleFreqMatrix_.clear();
  //historicalAlleleFreqMatrix_.clear();
  etaPriorMatrix_.clear();
  reportedAncestryMatrix_.clear();
}
