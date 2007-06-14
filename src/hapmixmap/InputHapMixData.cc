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
 
  double threshold = 100.0;//if(options->getHapMixModelIndicator())threshold /= options->getRhoPriorMean();
  checkLocusFile(genotypesSexColumn, threshold, options->CheckData());
  //locusMatrix_ = locusMatrix_.SubMatrix(1, locusMatrix_.nRows() - 1, 1, 2);//remove header and first column of locus file

  CheckAlleleFreqs(options, Log);
  ReadBlockStateLabels(options);

  //detects regression model
  if(strlen( options->getOutcomeVarFilename() ) || strlen( options->getCoxOutcomeVarFilename() )){//if outcome specified
    if ( strlen( options->getOutcomeVarFilename() ) != 0 )
      CheckOutcomeVarFile( options, Log);
    if ( strlen( options->getCoxOutcomeVarFilename() ) != 0 ){
      OutcomeType.push_back( CoxData );
	if(options->CheckData())
	  CheckCoxOutcomeVarFile( Log);
    }
    if ( strlen( options->getCovariatesFilename() ) != 0 )
      CheckCovariatesFile(Log);
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

void InputHapMixData::CheckForMonomorphicLoci(LogWriter& Log)const{
  hGenotypeLoader->CheckForMonomorphicLoci(Log);
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
 
void InputHapMixData::CheckOutcomeVarFile(Options* const options, LogWriter& Log){
  unsigned N = genotypeLoader->getNumberOfIndividuals();
  const unsigned NumCCIndividuals = hGenotypeLoader->getNumberOfCaseControlIndividuals();

  //check outcomevarfile and genotypes file have the same number of rows
  if(NumCCIndividuals>0 )
    N = NumCCIndividuals;
  if( outcomeVarMatrix_.nRows() - 1 != NumCCIndividuals ){
    stringstream s;
    s << "ERROR: Case-Control Genotypes file has " << NumCCIndividuals << " observations and Outcomevar file has "
      << outcomeVarMatrix_.nRows() - 1 << " observations.\n";
    throw(s.str());
  }

  //check the number of outcomes specified is not more than the number of cols in outcomevarfile
  int Firstcol = options->getTargetIndicator();
  int NumOutcomes = options->getNumberOfOutcomes();
  if(strlen(options->getCoxOutcomeVarFilename())) --NumOutcomes;

  int numoutcomes = NumOutcomes;
  if(NumOutcomes > -1){//options 'numberofregressions' used
    if((int)outcomeVarMatrix_.nCols() - Firstcol < NumOutcomes){
      numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;//adjusts if too large
      Log << "ERROR: 'outcomes' is too large, setting to " << numoutcomes;
    }
  }
  else numoutcomes = (int)outcomeVarMatrix_.nCols() - Firstcol;
  options->setNumberOfOutcomes(numoutcomes);

  RegressionType RegType = None;  
  if(numoutcomes >0){
    //RegType = Both;
    //extract portion of outcomevarfile needed
    std::string* OutcomeVarLabels = new string[ outcomeVarMatrix_.nCols() ];
    getLabels(outcomeVarData_[0], OutcomeVarLabels);
    DataMatrix Temp = outcomeVarMatrix_.SubMatrix(1, N, Firstcol, Firstcol+numoutcomes-1);
    outcomeVarMatrix_ = Temp;
    
    //determine type of outcome - binary/continuous
    for( int j = 0; j < numoutcomes; j++ ){
      OutcomeType.push_back( Binary );
      for(unsigned i = 0; i < N; ++i)
	if(!outcomeVarMatrix_.isMissing(i, j) && !(outcomeVarMatrix_.get( i, j ) == 0 || outcomeVarMatrix_.get( i, j ) == 1) ){
	  OutcomeType[j] =  Continuous ;
	  break;
	}
      //in this way, the outcome type is set as binary only if all individuals have outcome values of 1 or 0
      //otherwise, a continuous outcome of 1.0 or 0.0 could lead to the type being wrongly set to binary.
      
      //need to check for allmissing
      //     if(i == NumIndividuals){
      //       Log << "ERROR: all outcomes missing\n";
      //       exit(1);
      //     }
      
      Log << "Regressing on ";    
      if( OutcomeType[j] == Binary ){
	Log << "Binary variable: ";
	if(numoutcomes==1)RegType = Logistic;//one logistic regression
	else {
	  if(RegType == Logistic) RegType = Mlogistic;//more than one logistic regression
	  else RegType = Multiple;//linear and logistic 
	}
      }
      else if(OutcomeType[j] == Continuous ){
	Log << "Continuous variable: ";
	if(numoutcomes==1)RegType = Linear;//one linear regression
	else {
	  if(RegType == Linear) RegType = Mlinear;//more than one linear regression
	  else RegType = Multiple;//linear and logistic 
	}
      }
      Log << outcomeVarData_[0][j+Firstcol];
      Log << ".\n";
      OutcomeLabels.push_back(outcomeVarData_[0][j+Firstcol]);
    }
    Log << "\n";
    delete[] OutcomeVarLabels;
  }
  options->setRegType(RegType);
}

