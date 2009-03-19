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
#include "bclib/DataReader.h"
#include "bclib/StringConvertor.h"
#include <sstream>
#include <string.h>

InputHapMixData::InputHapMixData(HapMixOptions *options, bclib::LogWriter &Log){
  hGenotypeLoader = new HapMixGenotypeLoader;
  //assign pointer in base class to same address
  genotypeLoader = (GenotypeLoader*) hGenotypeLoader;

  Log.setDisplayMode(bclib::Quiet);
  try
    {
      ReadData(options, Log);
      hGenotypeLoader->ReadTestGenotypes(options->getTestGenotypesFilename(), Log);

    } catch (const exception& e) {
    cerr << "\nException occured during parsing of input file: \n" << e.what() << endl;
    exit(1);
  }
  catch(string s){
    cerr << "\nException occured during parsing of input file: \n" << s << endl;;
    exit(1);
  }
  NumSimpleLoci = getNumberOfSimpleLoci();
  NumCompositeLoci = determineNumberOfCompositeLoci();
  distanceUnit = DetermineUnitOfDistance();
  SetLocusLabels();
  ReadBlockStateLabels(options);

  CheckRegressionData(options, Log);
}

bool InputHapMixData::CheckData(HapMixOptions *options, bclib::LogWriter &Log){

  Log.setDisplayMode(bclib::Quiet);

  bool DataOK = true;

  try{
    DataOK = DataOK & genotypeLoader->CheckForUnobservedAlleles(locusMatrix_, Log);
  }catch(string s){
    Log << bclib::On << s << bclib::Quiet;
    DataOK = false;
  }
  DataOK = DataOK & checkLocusFile(options, Log);
  
  if( strlen( options->getPriorAlleleFreqFilename() )) {
    DataOK = DataOK & CheckAlleleFreqs(options, Log);
  }

  return DataOK;
}

/**
   checks outcome variables, covariates and extracts required variables, 
   detecting which type of regression (linear/logistic) is required.
   Exits if any problames encountered
*/
void InputHapMixData::CheckRegressionData(HapMixOptions *options, bclib::LogWriter &Log){

  //detects regression model
  if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
    unsigned N =  hGenotypeLoader->getNumberOfTestIndividuals();
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

///tells if a given locus is typed, in a testing analysis
bool InputHapMixData::isTypedLocus(unsigned locus)const{
  return hGenotypeLoader->isTypedLocus(locus);
}

//returns the number of typed loci in a testing analysis
unsigned InputHapMixData::getNumTypedLoci()const{
  return hGenotypeLoader->getNumTypedLoci();
}
int InputHapMixData::getNumberOfTestIndividuals()const {
  return hGenotypeLoader->getNumberOfTestIndividuals();
}

bool InputHapMixData::IsTestIndividual(int i)const{
  return hGenotypeLoader->IsTestIndividual(i);
}

void InputHapMixData::GetGenotype(int i, const Genome &Loci, 
			   std::vector<genotype>* genotypes, bool **Missing)const{
  hGenotypeLoader->GetGenotype(i, Loci, genotypes, Missing);
}
bool InputHapMixData::GetHapMixGenotype(int i, const Genome &Loci, 
				  std::vector<unsigned short>* genotypes, bool** Missing){
  return hGenotypeLoader->GetHapMixGenotype(i, Loci, genotypes, Missing);
}

void InputHapMixData::ReadBlockStateLabels(HapMixOptions *options){
  
  if(strlen(options->getPriorAlleleFreqFilename()))
    bclib::DataReader::ReadHeader(options->getPriorAlleleFreqFilename(), HiddenStateLabels);

  else{
    //set default pop labels
    for( int j = 0; j < options->getPopulations(); j++ ){
      stringstream poplabel;
      poplabel << "BlockState" << j+1;
      HiddenStateLabels.push_back(poplabel.str());
    }
  }
}

///checks the priorallelefreqfile has the correct number of rows
bool InputHapMixData::CheckAlleleFreqs(HapMixOptions *options, bclib::LogWriter &Log){
  unsigned int NumberOfStates = 2*NumCompositeLoci;
  
  //use the following block if allowing multiallelic loci
  //   unsigned index = 0;
  //   for(unsigned i = 0; i < NumCompositeLoci; ++i){
  //     int states = 1;
  
  //     do{
  //       states *= (int)locusMatrix_.get( index, 0 );
  //       index++;
  //     }
  //     while( index < locusMatrix_.nRows() && !locusMatrix_.isMissing(index, 1) && locusMatrix_.get( index, 1 ) == 0 );
  //     NumberOfStates += states;
  //   }
  
  if(getPriorAlleleFreqData().size() != NumberOfStates+1){
    Log << "Incorrect number of rows in priorallelefreqfile.\n" 
	<< "Expecting " << (int)NumberOfStates+1 << " rows, but there are " <<(int) getPriorAlleleFreqData().size() << " rows.\n";
    return false;
  }
  //set number of block states as the number of columns in priorallelefreqfile -1
  options->setPopulations(getPriorAlleleFreqData()[0].size() - 1);

  //if all is ok return true
  return true;
}
 
