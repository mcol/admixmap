/**
 *   InputData.cc 
 *   Class to read and check all input data files
 *   Copyright (c) 2005 - 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "InputData.h"
#include "Options.h"
#include "bclib/StringConvertor.h"
#include "bclib/DataReader.h"
#include "bclib/LogWriter.h"
#include "Genome.h"
#include <string>
#include <sstream>

using namespace std;
using namespace bclib;

/// Extracts population labels from header line of allelefreq input file
void InputData::getPopLabels(const Vector_s& data, size_t Populations, Vector_s& labels){
  if(data.size() != Populations+1){cout << "Error in getPopLabels\n";exit(1);}

  for (size_t i = 0; i < Populations; ++i) {
    labels.push_back( StringConvertor::dequote(data[i+1]) );
  }
}
void getLabels(const Vector_s& data, string *labels){
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = StringConvertor::dequote(data[i]);
  }
}

InputData::InputData(){
  NumSimpleLoci = 0;
  NumCompositeLoci = 0;
  genotypeLoader = 0;
}

InputData::~InputData(){
  delete genotypeLoader;
}

void InputData::ReadData(Options *options, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  try
    {
      // Read all input files.
      DataReader::ReadData(options->getLocusFilename(), locusData_, Log);   //locusfile
      //convert to DataMatrix, dropping header and first col and use only 2 cols
      DataReader::convertMatrix(locusData_, locusMatrix_, 1, 1,2);

      //read genotype data
      genotypeLoader->Read(options->getGenotypesFilename(), locusData_.size() - 1, Log);
      
      DataReader::ReadData(options->getCovariatesFilename(), covariatesData_, covariatesMatrix_,Log);     //covariates file
      DataReader::ReadData(options->getOutcomeVarFilename(), outcomeVarData_,outcomeVarMatrix_, Log);//outcomevar file
      DataReader::ReadData(options->getCoxOutcomeVarFilename(), coxOutcomeVarData_, Log);            //coxoutcomevar file
      DataReader::convertMatrix(coxOutcomeVarData_, coxOutcomeVarMatrix_, 1, 0,0);//drop first row in conversion
      
      DataReader::ReadData(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_, Log);

      Log << "\n";
      
    } catch (const exception& e) {
    cerr << "\nException occured during parsing of input file: \n" << e.what() << endl;
    exit(1);
  }
  catch(string s){
    cerr << "\nException occured during parsing of input file: \n" << s << endl;;
    exit(1);
  }
 
}

///determine number of individuals by counting lines in genotypesfile 
int InputData::getNumberOfIndividuals()const {
  return genotypeLoader->getNumberOfIndividuals();
}

///determine number of loci by counting rows of locusfile
int InputData::getNumberOfSimpleLoci()const {
  return(locusData_.size() - 1);
}
///determines number of composite loci from locusfile
unsigned InputData::determineNumberOfCompositeLoci()const{
  unsigned NumberOfCompositeLoci = locusMatrix_.nRows();
  for( unsigned i = 0; i < locusMatrix_.nRows(); i++ )
    if( !locusMatrix_.isMissing(i,1) && locusMatrix_.get( i, 1 ) == 0.0 ) NumberOfCompositeLoci--;
  return NumberOfCompositeLoci;
}

///determine unit of distance from locus file header. Defaults to Morgans if not specified.
GeneticDistanceUnit InputData::DetermineUnitOfDistance(){
  GeneticDistanceUnit u = Morgans;//default, usual for admixture mapping
  string distance_header = locusData_[0][2];
  if(distance_header.find("cm")!=string::npos || distance_header.find("CM")!=string::npos 
     || distance_header.find("cM")!=string::npos) 
    u = centimorgans;
  else if(distance_header.find("mb")!=string::npos || distance_header.find("Mb")!=string::npos 
     || distance_header.find("MB")!=string::npos) 
    u = megabases;
  
  return u;
}

GeneticDistanceUnit InputData::getUnitOfDistance()const{
  return distanceUnit;
}
const string& InputData::getUnitOfDistanceAsString()const{
  return GeneticDistanceUnitString[getUnitOfDistance()];
}

bool InputData::checkLocusFile(LogWriter& Log){
  bool badData = false;

  const vector<string>& GenotypesFileHeader = genotypeLoader->getHeader();
  const float threshold = getLocusDistanceThreshold();

  for (unsigned i = 1; i < locusData_.size(); ++i) {//rows of locusfile

    const float distance = locusMatrix_.get(i-1,1);
    const string locusName = locusData_[i][0];

    //check number of alleles is >1

    if(locusMatrix_.get(i-1,0) <2){
      Log << On << "ERROR on line " << i << " of locusfile: number of alleles must be >1.\n";
      badData = true;
    }
    
    //check distances are not negative
    if( distance < 0.0){
      badData = true;
      Log << On << "Error: distance on line "<< i <<" of locusfile is negative.\n";
    }
    //check distances are not too large 
    if(distance >= threshold) {
      //badData = true;
      if(distance != 100 )//for backward-compatibility; no warning if 100 used to denote new chromosome      
	Log << On << "Warning: distance of " << distance << " " << GeneticDistanceUnitString[distanceUnit] << "  at locus " << i << "\n";
      locusMatrix_.isMissing(i-1,1, true);//missing value for distance denotes new chromosome
    }
    
    // Check loci names are unique    
    for (size_t j = 0; j < i-1; ++j) {   
      if (locusName == locusData_[j][0]) {
	badData = true;
	Log << On << "Error in locusfile. Two different loci have the same name: "
	    << locusName << "\n";
      }
    }
    // Compare loci names in locus file and genotypes file.
    if (StringConvertor::dequote(locusName) != StringConvertor::dequote(GenotypesFileHeader[i + genotypeLoader->getSexColumn()])) {
      Log << On << "Error. Locus names in locus file and genotypes file are not the same.\n"
	  << "Locus names causing an error are: " << locusName << " and " 
	  << GenotypesFileHeader[i + genotypeLoader->getSexColumn()] << "\n";
      return false;
    }
    
  }//end loop over loci
  return !badData;

}

///determines the distance threshold for a new chromosome
//100 Morgans or 10 Mb
float InputData::getLocusDistanceThreshold()const{
  switch(distanceUnit){
  case basepairs:{
    return (1e7);
    }
  case kilobases:{
    return(10000.0);
    }
  case megabases:{
    return(10.0);
    }
  case centimorgans:{
    return(10000.0);
  }
  case Morgans:{
    return (100.0);
  }
  default:{
    return(100.0);
   }
  }
}

void InputData::SetLocusLabels(){
  for (size_t i = 1; i < locusData_.size(); ++i) {//rows of locusfile
    LocusLabels.push_back(StringConvertor::dequote(locusData_[i][0]));
  }
}

void InputData::CheckOutcomeVarFile(unsigned N, Options* const options, LogWriter& Log){
  if( outcomeVarMatrix_.nRows() - 1 !=  N){
    stringstream s;
    s << "ERROR: Genotypes file has " << N << " observations and Outcomevar file has "
      << outcomeVarMatrix_.nRows() - 1 << " observations.\n";
    throw(s.str());
  }

  const vector<unsigned>& OutcomeVarCols = options->getOutcomeVarColumns();
  vector<unsigned> ColsInUse;
  //check all in range and unique
  if(OutcomeVarCols.size()){
    for(vector<unsigned>::const_iterator i = OutcomeVarCols.begin(); i != OutcomeVarCols.end(); ++i){
      if(*i > 0 && *i <= outcomeVarData_[0].size())
	ColsInUse.push_back(*i -1);
      else{
	Log << On << "Warning: there is no column " << (unsigned)*i << " in " << options->getOutcomeVarFilename() << "\n";
      }
    }
    // sort(ColsInUse.begin(), ColsInUse.end());
    vector<unsigned>::iterator new_end = unique(ColsInUse.begin(), ColsInUse.end());//reorder, putting duplicates at end
    ColsInUse.erase(new_end, ColsInUse.end());//remove duplicates
  }
  else{//use all if no outcomevarcols options specified
    for(unsigned i = 0; i < outcomeVarData_[0].size(); ++i)
      ColsInUse.push_back(i);
  }

  options->setNumberOfOutcomes(ColsInUse.size());

  RegressionType RegType = None;
  DataMatrix TempOutcome(N, ColsInUse.size());  
  unsigned col = 0;
  for(vector<unsigned>::const_iterator ov = ColsInUse.begin(); ov != ColsInUse.end(); ++ov, ++col){
    bool isContinuous = false;
    for(unsigned i = 0; i < N; ++i){
      TempOutcome.set(i, col, outcomeVarMatrix_.get( i+1, *ov ));
      TempOutcome.isMissing(i, col, outcomeVarMatrix_.isMissing( i+1, *ov ));
      //in this way, the outcome type is set as binary only if all individuals have outcome values of 1 or 0
      //otherwise, a continuous outcome of 1.0 or 0.0 could lead to the type being wrongly set to binary.
      if(!outcomeVarMatrix_.isMissing(i+1, *ov) && !(outcomeVarMatrix_.get( i+1, *ov ) == 0 || outcomeVarMatrix_.get( i+1, *ov ) == 1) ){
	isContinuous = true;
      }
    }
    if(isContinuous)
      OutcomeType.push_back( Continuous );
    else
      OutcomeType.push_back( Binary );
    Log << "Regressing on ";    
    if( !isContinuous ){
      Log << "Binary variable: ";
      if(ColsInUse.size()==1)RegType = Logistic;//one logistic regression
      else {
	if(RegType == Logistic) RegType = Mlogistic;//more than one logistic regression
	else RegType = Multiple;//linear and logistic 
      }
    }
    else {
      Log << "Continuous variable: ";
      if(ColsInUse.size()==1)RegType = Linear;//one linear regression
      else {
	if(RegType == Linear) RegType = Mlinear;//more than one linear regression
	else RegType = Multiple;//linear and logistic 
      }
    }
    Log << outcomeVarData_[0][*ov];
    Log << ".\n";
    OutcomeLabels.push_back(StringConvertor::dequote(outcomeVarData_[0][*ov]));
  }
  outcomeVarMatrix_.clear();
  outcomeVarMatrix_ = TempOutcome;
   
  options->setRegType(RegType);
}

void InputData::CheckCoxOutcomeVarFile(LogWriter &Log)const{
  if(coxOutcomeVarMatrix_.nCols() !=3){
    Log << "ERROR: 'coxoutcomevarfile should have 3 columns but has " << coxOutcomeVarMatrix_.nCols() << "\n";
    exit(1);
  }
  if( coxOutcomeVarMatrix_.nRows() != genotypeLoader->getNumberOfIndividuals() ){
    stringstream s;
    s << "ERROR: Genotypes file has " << genotypeLoader->getNumberOfIndividuals() 
      << " observations and coxoutcomevar file has "
      << coxOutcomeVarMatrix_.nRows() - 1 << " observations.\n";
    throw(s.str());
  }

  bool flag = false;
  for(unsigned i = 0; i < coxOutcomeVarMatrix_.nRows(); ++i)
    if(coxOutcomeVarMatrix_.get(i, 0) >= coxOutcomeVarMatrix_.get(i, 1)){
      Log << "Error in coxoutcomevarfile on line " <<i << " : finish times must be later than start times\n";
      flag = true;
    }
  if(flag)exit(1);

}

void InputData::CheckCovariatesFile(unsigned NumIndividuals, Options* const options, LogWriter &Log){
  if( NumIndividuals != covariatesMatrix_.nRows() - 1 ){
    Log << "ERROR: Genotypes file has " << NumIndividuals 
	<< " observations and Covariates file has "
	<< covariatesMatrix_.nRows() - 1 << " observations.\n";
    exit(1);
  }

  const vector<unsigned>& CovariateCols = options->getCovariateColumns();
  vector<unsigned> ColsInUse;
  //check all in range and unique
  if(CovariateCols.size()){
    for(vector<unsigned>::const_iterator i = CovariateCols.begin(); i != CovariateCols.end(); ++i){
      if(*i > 0 && *i <= covariatesData_[0].size())
	ColsInUse.push_back(*i -1);
      else{
	Log << On << "Warning: there is no column " << (unsigned)*i << " in " << options->getCovariatesFilename() << "\n";
      }
    }
    // sort(ColsInUse.begin(), ColsInUse.end());
    vector<unsigned>::iterator new_end = unique(ColsInUse.begin(), ColsInUse.end());//reorder, putting duplicates at end
    ColsInUse.erase(new_end, ColsInUse.end());//remove duplicates
  }
  else{//use all if no covariatecols options specified
    for(unsigned i = 0; i < covariatesData_[0].size(); ++i)
      ColsInUse.push_back(i);
  }

  DataMatrix Temp(NumIndividuals, ColsInUse.size());  
  unsigned col = 0;
  for(vector<unsigned>::const_iterator cv = ColsInUse.begin(); cv != ColsInUse.end(); ++cv, ++col){
    for(unsigned i = 0; i < NumIndividuals; ++i){
      Temp.set(i, col, covariatesMatrix_.get( i+1, *cv ));
      Temp.isMissing(i, col, covariatesMatrix_.isMissing( i+1, *cv ));
    }
      //set covariate labels
      CovariateLabels.push_back(StringConvertor::dequote(covariatesData_[0][*cv]));
  }
    covariatesMatrix_.clear();
    covariatesMatrix_ = Temp;
}


///determines if an individual is female
bool InputData::isFemale(int i)const{
  return genotypeLoader->isFemale(i);
}

void InputData::getOutcomeTypes(DataType* T)const{
  for(unsigned i = 0; i < OutcomeType.size(); ++i)
    T[i] = OutcomeType[i];
}
DataType InputData::getOutcomeType(unsigned i)const{
  return OutcomeType[i];
}
const Matrix_s& InputData::getLocusData() const
{
  return locusData_;
}

const Matrix_s& InputData::getCovariatesData() const
{
  return covariatesData_;
}

const Matrix_s& InputData::getOutcomeVarData() const
{
  return outcomeVarData_;
}

const Matrix_s& InputData::getPriorAlleleFreqData() const
{
  return priorAlleleFreqData_;
}

const DataMatrix& InputData::getLocusMatrix() const
{
  return locusMatrix_;
}

// const DataMatrix& InputData::getPriorAlleleFreqMatrix() const
// {
//     return priorAlleleFreqMatrix_;
// }

const DataMatrix& InputData::getOutcomeVarMatrix() const
{
  return outcomeVarMatrix_;
}
const DataMatrix& InputData::getCoxOutcomeVarMatrix() const
{
  return coxOutcomeVarMatrix_;
}

const DataMatrix& InputData::getCovariatesMatrix() const
{
  return covariatesMatrix_;
}
const Vector_s& InputData::GetHiddenStateLabels() const{
  return HiddenStateLabels;
}
Vector_s InputData::getOutcomeLabels()const{
  return OutcomeLabels;
}
const Vector_s& InputData::getLocusLabels()const{
  return LocusLabels;
}
const Vector_s InputData::getCovariateLabels()const{
  return CovariateLabels;
}

// const string& InputData::getIndivId(unsigned i)const{
//   return geneticData_[i+1][0];//+1 to skip header
// }

void InputData::Delete(){
  //erase string matrices
  genotypeLoader->clear();

  for(unsigned i = 0; i < locusData_.size(); ++i)
    locusData_[i].clear();
  locusData_.clear();
  for(unsigned i = 0; i < covariatesData_.size(); ++i)
    covariatesData_[i].clear();
  covariatesData_.clear();
  for(unsigned i = 0; i < outcomeVarData_.size(); ++i)
    outcomeVarData_[i].clear();
  outcomeVarData_.clear();
  for(unsigned i = 0; i < coxOutcomeVarData_.size(); ++i)
    coxOutcomeVarData_[i].clear();
  coxOutcomeVarData_.clear();
  for(unsigned i = 0; i < priorAlleleFreqData_.size(); ++i)
    priorAlleleFreqData_[i].clear();
  priorAlleleFreqData_.clear();


  //erase data matrices 
  locusMatrix_.clear();
  covariatesMatrix_.clear();
  outcomeVarMatrix_.clear();
  coxOutcomeVarMatrix_.clear();
  //priorAlleleFreqMatrix_.clear();
}



