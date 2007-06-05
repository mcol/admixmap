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
#include "bcppcl/StringConvertor.h"
#include "bcppcl/DataReader.h"
#include "Genome.h"
#include <string>
#include <sstream>

using namespace std;

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
      DataReader::convertMatrix(locusData_, locusMatrix_, 1, 1,2);//drop first row, first col and last col

      //read genotype data
      genotypeLoader->Read(options->getGenotypesFilename(), Log);
      
      DataReader::ReadData(options->getCovariatesFilename(), inputData_, covariatesMatrix_,Log);     //covariates file
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

void InputData::checkLocusFile(int sexColumn, double threshold, bool check){
  // if check = true, Checks that loci labels in locusfile are unique and that they match the names in the genotypes file.
  //also extracts locus labels
  bool flag = false;

  if(getUnitOfDistance()==centimorgans)threshold *= 100.0;
  for (size_t i = 1; i < locusData_.size(); ++i) {//rows of locusfile
    if(check){
      //check number of alleles is >1
      if(locusMatrix_.get(i-1,0) <2){
	cerr << "ERROR on line " << i+1 << " of locusfile: number of alleles must be >1." << endl;
	flag = true;
      }

      //check distances are not negative
      if(locusMatrix_.get(i-1,1) < 0.0){
	flag = true;
	cerr<<"Error: distance on line "<<i<<" of locusfile is negative."<<endl;
      }
      //check distances are not too large 
      if(locusMatrix_.get(i-1,1) >= threshold) {
	//flag = true;
	if(locusMatrix_.get(i-1,1) != 100 )//for backward-compatibility; no warning if 100 used to denote new chromosome      
	  cerr << "Warning: distance of " <<locusMatrix_.get(i-1,1)<< "  at locus " <<i<<endl;
	locusMatrix_.isMissing(i-1,1, true);//missing value for distance denotes new chromosome
      }
      
      // Check loci names are unique    
      for (size_t j = 0; j < i-1; ++j) {   
	if (locusData_[i][0] == locusData_[j][0]) {
	  flag = true;
	  cerr << "Error in locusfile. Two different loci have the same name: "
	       << locusData_[i][0] << endl;
	}
      }

    }//end if check
    LocusLabels.push_back(StringConvertor::dequote(locusData_[i][0]));
  }//end loop over loci
  if(flag)exit(1);
  if(check){
    const size_t numLoci = locusData_.size() - 1;//number of simple loci
    const vector<string>& GenotypesFileHeader = genotypeLoader->getHeader();

    // Compare loci names in locus file and genotypes file.
    for (size_t i = 1; i <= numLoci; ++i) {
      if (StringConvertor::dequote(locusData_[i][0]) != StringConvertor::dequote(GenotypesFileHeader[i + sexColumn])) {
	cout << "Error. Locus names in locus file and genotypes file are not the same." << endl;
	cout << "Locus names causing an error are: " << locusData_[i][0] << " and " 
	     << GenotypesFileHeader[i + sexColumn] << endl;
	//cout << options->getgenotypesSexColumn() << endl;
	exit(2);
      }
    }
  } 
}

///checks number of loci in genotypes file is the same as in locusfile, 
///and determines if there is a sex column
void InputData::DetermineSexColumn(){
  const size_t numLoci = locusData_.size() - 1; //number of loci in locus file
  int SexCol = 0;
  // Determine if "Sex" column present in genotypes file.
  if (numLoci == genotypeLoader->NumLoci() - 1) {
    SexCol = 0;//no sex col
  } else if (numLoci == genotypeLoader->NumLoci() - 2) {
    SexCol  = 1;//sex col
  } else {//too many cols
    cerr << "Error: " << numLoci << " loci in locus file but " 
	 <<  genotypeLoader->NumLoci() - 1 << " loci in genotypes file." << endl;
    exit(2);
  }
  genotypesSexColumn = SexCol;
}
 
void InputData::CheckOutcomeVarFile(Options* const options, LogWriter& Log){
  unsigned N = genotypeLoader->getNumberOfIndividuals();

  if( outcomeVarMatrix_.nRows() - 1 != (N - options->getTestOneIndivIndicator()) ){
    stringstream s;
    s << "ERROR: Genotypes file has " << N << " observations and Outcomevar file has "
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

void InputData::CheckCovariatesFile(LogWriter &Log){
  if( genotypeLoader->getNumberOfIndividuals() != covariatesMatrix_.nRows() - 1 ){
    Log << "ERROR: Genotypes file has " << genotypeLoader->getNumberOfIndividuals() 
	<< " observations and Covariates file has "
	<< covariatesMatrix_.nRows() - 1 << " observations.\n";
    exit(1);
  }
  for (size_t i = 0; i < inputData_[0].size(); ++i) {
    CovariateLabels.push_back(StringConvertor::dequote(inputData_[0][i]));
  }
}


///determines if an individual is female
bool InputData::isFemale(int i)const{
  if (genotypesSexColumn == 1)
    return genotypeLoader->isFemale(i);
  else
    return false;
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

const Matrix_s& InputData::getInputData() const
{
  return inputData_;
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
  for(unsigned i = 0; i < inputData_.size(); ++i)
    inputData_[i].clear();
  inputData_.clear();
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



