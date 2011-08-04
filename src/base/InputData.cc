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
/// Definition of the InputData class.
//=============================================================================

#include "InputData.h"
#include "Options.h"
#include "bclib/StringConvertor.h"
#include "bclib/DataReader.h"
#include "bclib/LogWriter.h"
#include "DataValidError.h"
#include "GenotypeParser.h"
#include "SimpleLocusParser.h"

#include <algorithm>  // for unique()
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace bclib;
using namespace genepi;


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
  genotypeParser = 0;
}

InputData::~InputData(){
  delete genotypeParser;
}

void InputData::ReadData( Options * options, LogWriter & Log )
    {
    Log.setDisplayMode(Quiet);

    //read genotype data
    SimpleLocusParser::parse( options->getLocusFilename(), simpleLoci );
    genotypeParser = new GenotypeParser(options->getGenotypesFilename(),
                                        simpleLoci,
                                        options->getOutcomeIsBinary(),
                                        options->getMaleXHomozygWarn());

    // We used to generate the pedigrees here, but it turns out the number
    // of populations is not yet known.  This is now called from
    // InputAdmixData's constructor.  (except that turned out to also be too
    // soon, so it is called from InputAdmixData::finishConstructing()).
    #if 0
      generatePedigrees( options );
    #endif

    if (options->getWarningsAreErrors() && genotypeParser->getNWarnings()) {
      std::cerr << "Warnings treated as errors: aborting execution.\n";
      exit(1);
    }

    DataReader::ReadData(options->getCovariatesFilename(), covariatesData_, covariatesMatrix_,Log);	//covariates file
    DataReader::ReadData(options->getOutcomeVarFilename(), outcomeVarData_,outcomeVarMatrix_, Log);//outcomevar file
    DataReader::ReadData(options->getCoxOutcomeVarFilename(), coxOutcomeVarData_, Log);		   //coxoutcomevar file
    DataReader::convertMatrix(coxOutcomeVarData_, coxOutcomeVarMatrix_, 1, 0,0);//drop first row in conversion

    DataReader::ReadData(options->getPriorAlleleFreqFilename(), priorAlleleFreqData_, Log);

    Log << "\n";
  }



//-------------------------------------------------------------------------
// generatePedigrees()
//-------------------------------------------------------------------------

void InputData::generatePedigrees( const Options & options )
    {
    if ( isPedFile() || options.getUsePedForInd() )
	genepi::Pedigree::generatePedigrees( *genotypeParser, peds,
					     options.getMaxPedigreeSize() );
    }



//-------------------------------------------------------------------------
/// Determine number of individuals by counting lines in genotypesfile
//-------------------------------------------------------------------------
size_t InputData::getNumberOfIndividuals()const {
  return genotypeParser->getNumOrganisms();
}


GeneticDistanceUnit InputData::getUnitOfDistance()const{
  return simpleLoci.getGDU();
}


void InputData::SetLocusLabels(){
  for (size_t i = 0; i < simpleLoci.size(); ++i) // rows of locusfile
    LocusLabels.push_back( simpleLoci[i].getName() );
}

void InputData::CheckOutcomeVarFile(unsigned N, Options* const options, LogWriter& Log){
  if( outcomeVarMatrix_.nRows() - 1 !=	N){
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
    Log << "Regressing on ";
    if( !isContinuous ){
      OutcomeType.push_back( Binary );
      Log << "Binary variable: ";
      if(ColsInUse.size()==1)RegType = Logistic;//one logistic regression
      else {
	if(RegType == Logistic) RegType = Mlogistic;//more than one logistic regression
	else RegType = Multiple;//linear and logistic
      }
    }
    else {
      OutcomeType.push_back( Continuous );
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

  // Pedigree files contain an outcome column; in the case of genotype files, we
  // must "merge" the data from the outcome file back into the organism objects:
  if ( ! isPedFile() ) {
    // For now, just use the first column.
    const unsigned col_to_use = 0;
    if ( outcomeVarMatrix_.nCols() != 0 ) {
      gp_assert_eq(outcomeVarMatrix_.nRows(), genotypeParser->getNumOrganisms());
      for ( unsigned row = outcomeVarMatrix_.nRows() ; row-- != 0 ; )
        genotypeParser->getOrganism(row).setOutcome(
            // So ugly: cast double to integer to enumerated value.
	    static_cast<Organism::OutcomeType>(int(outcomeVarMatrix_.get(row,col_to_use))) );
    }
  }
}


void InputData::CheckCoxOutcomeVarFile(LogWriter &Log)const{
  if(coxOutcomeVarMatrix_.nCols() !=3){
    Log << "ERROR: 'coxoutcomevarfile should have 3 columns but has " << coxOutcomeVarMatrix_.nCols() << "\n";
    exit(1);
  }
  if (coxOutcomeVarMatrix_.nRows() != genotypeParser->getNumberOfIndividuals()) {
    stringstream s;
    s << "ERROR: Genotypes file has " << genotypeParser->getNumberOfIndividuals()
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
  return (*genotypeParser)[i].isFemale();
}

void InputData::getOutcomeTypes(DataType* T)const{
  for(unsigned i = 0; i < OutcomeType.size(); ++i)
    T[i] = OutcomeType[i];
}
DataType InputData::getOutcomeType(unsigned i)const{
  return OutcomeType[i];
}


const Matrix_s& InputData::getCovariatesData() const
{
  return covariatesData_;
}

const Matrix_s& InputData::getOutcomeVarData() const
{
  return outcomeVarData_;
}


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


const Vector_s InputData::getCovariateLabels()const{
  return CovariateLabels;
}

// const string& InputData::getIndivId(unsigned i)const{
//   return geneticData_[i+1][0];//+1 to skip header
// }

/// Erase string matrices
void InputData::Delete(){

  // Make this a data-member rather than a pointer when the parser and
  // container are separated:
  delete genotypeParser;
  genotypeParser = 0;

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
  covariatesMatrix_.clear();
  outcomeVarMatrix_.clear();
  coxOutcomeVarMatrix_.clear();
  //priorAlleleFreqMatrix_.clear();
}
