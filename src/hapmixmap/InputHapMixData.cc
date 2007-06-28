/**
 *   HAPMIXMAP 
 *   InputHapMixData.cc 
 *   Class to read HAPMIXMAP data
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "InputHapMixData.h"
#include "HapMixOptions.h"
#include "bcppcl/DataReader.h"
#include "bcppcl/StringConvertor.h"
#include <sstream>

InputHapMixData::InputHapMixData(HapMixOptions *options, LogWriter &Log){
  hGenotypeLoader = new HapMixGenotypeLoader;
  //assign pointer in base class to same address
  genotypeLoader = (GenotypeLoader*) hGenotypeLoader;

  Log.setDisplayMode(Quiet);
  try
    {
      ReadData(options, Log);
      hGenotypeLoader->ReadCaseControlGenotypes(options->getCCGenotypesFilename(), Log);

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

void InputHapMixData::CheckData(HapMixOptions *options, LogWriter &Log){
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  distanceUnit = DetermineUnitOfDistance();

  Log.setDisplayMode(Quiet);
  DetermineSexColumn();
  if(options->CheckData())
    if(!genotypeLoader->CheckForUnobservedAlleles(locusMatrix_, genotypesSexColumn, Log))
      exit(1);
 
  bool badData = false;
  if(options->CheckData())
    badData = !checkLocusFile(Log);

  if(badData)
    exit(1);

  SetLocusLabels();
 
  CheckAlleleFreqs(options, Log);
  ReadBlockStateLabels(options);

  //detects regression model
  if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
    unsigned N =  hGenotypeLoader->getNumberOfCaseControlIndividuals();
    if(N == 0) N = genotypeLoader->getNumberOfIndividuals();

    if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
      CheckOutcomeVarFile(N, options, Log);
    }
    if ( strlen( options->getCoxOutcomeVarFilename() ) != 0 ){
      OutcomeType.push_back( CoxData );
	if(options->CheckData())
	  CheckCoxOutcomeVarFile( Log);
    }
    if ( strlen( options->getCovariatesFilename() ) != 0 )
      CheckCovariatesFile(N, options, Log);
  }

}

///tells if a given locus is typed, in a hapmix case-control analysis
bool InputHapMixData::isTypedLocus(unsigned locus)const{
  return hGenotypeLoader->isTypedLocus(locus);
}

//returns the number of typed loci in a hapmix case-control analysis
unsigned InputHapMixData::getNumTypedLoci()const{
  return hGenotypeLoader->getNumTypedLoci();
}
int InputHapMixData::getNumberOfCaseControlIndividuals()const {
  return hGenotypeLoader->getNumberOfCaseControlIndividuals();
}

bool InputHapMixData::IsCaseControl(int i)const{
  return hGenotypeLoader->IsCaseControl(i);
}

void InputHapMixData::GetGenotype(int i, const Genome &Loci, 
			   std::vector<genotype>* genotypes, bool **Missing)const{
  hGenotypeLoader->GetGenotype(i, genotypesSexColumn, Loci, genotypes, Missing);
}
bool InputHapMixData::GetHapMixGenotype(int i, const Genome &Loci, 
				  std::vector<unsigned short>* genotypes, bool** Missing){
  return hGenotypeLoader->GetHapMixGenotype(i, genotypesSexColumn, Loci, genotypes, Missing);
}

void InputHapMixData::ReadBlockStateLabels(HapMixOptions *options){
  
  if(strlen(options->getPriorAlleleFreqFilename()))
    DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), HiddenStateLabels);

  else{
    //set default pop labels
    for( int j = 0; j < options->getPopulations(); j++ ){
      stringstream poplabel;
      poplabel << "BlockState" << j+1;
      HiddenStateLabels.push_back(poplabel.str());
    }
  }
}

////checks consistency of supplied allelefreqs with locusfile
///and determines number of populations and population labels.
void InputHapMixData::CheckAlleleFreqs(HapMixOptions *options, LogWriter &Log){
  bool infile = false;//indicates whether a priorallelefreqfile has been specified
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

  //prior allelefreqs
  if( strlen( options->getPriorAlleleFreqFilename() )) {
     infile = true;
    nrows = priorAlleleFreqData_.size();
    expectednrows = NumberOfStates+1;
    Populations = priorAlleleFreqData_[0].size() - 1;
  }
  if(infile){
    if(nrows != expectednrows){
      Log << "Incorrect number of rows in priorallelefreqfile.\n" 
	  << "Expecting " << expectednrows << " rows, but there are " << nrows << " rows.\n";
      exit(1);
    }
    options->setPopulations(Populations);
  }
  //  else{//'states' option
    //     for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    //       if(Loci->GetNumberOfStates(i) < 2){
    // 	Log << "ERROR: The number of alleles at a locus is " << Loci->GetNumberOfStates(i) << "\n";
    // 	exit(1);
    //       }
    //     }
  //  }
}
 
