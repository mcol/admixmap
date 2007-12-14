/** 
 *   ADMIXMAP
 *   InputAdmixData.cc 
 *   class to read ADMIXMAP data
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "InputAdmixData.h"
#include "AdmixOptions.h"
#include "bclib/DataReader.h"
#include <sstream>

using bclib::LogWriter;

InputAdmixData::InputAdmixData(AdmixOptions *options, LogWriter &Log){
  using bclib::DataReader;
  genotypeLoader = new GenotypeLoader;

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
    
  } catch (const exception& e) {
    cerr << "\nException occured during parsing of input file: \n" << e.what() << endl;
    exit(1);
  }
  catch(string s){
    cerr << "\nException occured during parsing of input file: \n" << s << endl;;
    exit(1);
  }
 
  CheckData(options, Log);
}

void InputAdmixData::CheckData(AdmixOptions *options, LogWriter &Log){
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  distanceUnit = DetermineUnitOfDistance();

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

void InputAdmixData::GetGenotype(int i, const Genome &Loci, 
			   std::vector<genotype>* genotypes, bool **Missing)const{
  genotypeLoader->GetGenotype(i, Loci, genotypes, Missing);
}

void InputAdmixData::ReadPopulationLabels(AdmixOptions *options){
  using bclib::DataReader;
  //  if(strlen(options->getAlleleFreqFilename()) || strlen(options->getPriorAlleleFreqFilename()) || strlen(options->getHistoricalAlleleFreqFilename())){
  if(strlen(options->getAlleleFreqFilename()))
    DataReader::ReadHeader(options->getAlleleFreqFilename(), HiddenStateLabels);
  else if(strlen(options->getPriorAlleleFreqFilename()))
    DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), HiddenStateLabels);
  else if(strlen(options->getHistoricalAlleleFreqFilename()))
    DataReader::ReadHeader(options->getHistoricalAlleleFreqFilename(), HiddenStateLabels);
  else{
    //read labels from 'poplabels' option
    bclib::StringSplitter::Tokenize(options->getPopLabelString(), HiddenStateLabels, " ,");
    if(HiddenStateLabels.size() != (unsigned)options->getPopulations()){
      HiddenStateLabels.clear();

      //set default pop labels
      for( int j = 0; j < options->getPopulations(); j++ ){
	stringstream poplabel;
	poplabel << "Pop" << j+1;
	HiddenStateLabels.push_back(poplabel.str());
      }
    }
  }
}
 
void InputAdmixData::CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log){
  string freqtype = "";
  bool infile = false;//indicates whether either of the three allelefreq files are specified
  int nrows=0, expectednrows=0;
  int Populations = options->getPopulations();
  int NumberOfStates = 0;

  unsigned index = 0;
  for(unsigned i = 0; i < NumCompositeLoci; ++i){
    int states = 1;

    do{
      states *= (int)locusMatrix_.get( index, 0 );
      index++;
    }
    while( index < locusMatrix_.nRows() && !locusMatrix_.isMissing(index, 1) && locusMatrix_.get( index, 1 ) == 0 );
    NumberOfStates += states;
  }


  //fixed allele freqs
  if( strlen( options->getAlleleFreqFilename() )){
    freqtype = "";
    infile = true;
    nrows = alleleFreqData_.size()-1;
    expectednrows = NumberOfStates-NumCompositeLoci;
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
    nrows = priorAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = priorAlleleFreqData_[0].size() - 1;
    //getPopLabels(priorAlleleFreqData_[0], Populations, PopulationLabels);
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
	<<"Genotypesfile has " << genotypeLoader->getNumberOfIndividuals() << " rows\n";
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


