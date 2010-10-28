/*
 *   Options.cc 
 *   Base for classes to hold program options
 *   Copyright (c) 2006, 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

//=============================================================================
/// \file Options.cc
/// Implementation of the Options class.
//=============================================================================

#include "Options.h"
#include "bclib/StringConvertor.h"
#include "bclib/StringSplitter.h"
#include <string.h>
#include <sstream>
#include "Filenames.h"

using namespace std;
using namespace bclib;

Options::Options(){
  SetDefaultValues();
  DefineOptions();
}
bool Options::ReadUserOptions(int argc, char** argv, const char* fileargIndicator ){
  //one arg without leading '-'
  if(argc == 2 && strncmp(argv[1], "-", 1)){
    if(!ReadArgsFromFile(argv[1]))
      return false;
  }
  else if(argc > 1){
    //NOTE: command line args will not be supported soon
    //cout << "Warning: command-line arguments are deprecated" << endl;
   if(!ReadCommandLineArgs(argc, argv, fileargIndicator))
     return false;
  }
  return true;
}
void Options::SetDefaultValues(){
  // global variables store option values: variable names not necessarily same as 
  // command-line option names which are all lower-case 
  burnin = 100;
  TotalSamples = 1100;
  SampleEvery = 10;
  Seed = 1;
  displayLevel = 2; 
  checkData = true;
  thermoIndicator = false; // calculate marginal likelihood by thermodynamic method
  NumAnnealedRuns = 20; // default if option thermo not specified 
  ResultsDir = "results";
  LogFilename = "log.txt";
  fixedallelefreqs = false;
  NumberOfOutcomes = 0;
  RegType = None;
  TestForAllelicAssociation = false;
  TestForResidualAllelicAssoc = false;
  HWTest = false;
  OutputAlleleFreq = false;
  regressionPriorPrecision = 0.25;
  EYFilename = "ExpectedOutcomes.txt";
  DeleteOldResultsIndicator = true;

  useroptions["burnin"] = "100";
  useroptions["samples"] = "1100";
  useroptions["every"] = "10";
  useroptions["numannealedruns"] = "20";
  useroptions["displaylevel"] = "2";
  useroptions["logfile"] = "log.txt";
  useroptions["resultsdir"] = "results";
  useroptions["thermo"] = "0";
  useroptions["seed"] = "1";
  useroptions["checkdata"]="1";
  useroptions["fixedallelefreqs"] = "0";
  useroptions["regressionpriorprecision"] = "0.25";
}

Options::~Options()
{
}

// each option has a function to return its value
const string Options::getResultsDir() const{
  return ResultsDir;
}

long Options::getBurnIn() const
{
  return burnin;
}

long Options::getTotalSamples() const
{
  return TotalSamples;
}
int Options::getDisplayLevel() const
{
  return displayLevel;
}
long Options::getSampleEvery() const
{
  return SampleEvery;
}
long Options::getNumAnnealedRuns()const{
  return NumAnnealedRuns;
}
bool Options::getThermoIndicator()const{
  return thermoIndicator;
}
const char *Options::getLogFilename() const
{
  return LogFilename.c_str();
}

long Options::getSeed() const
{
  return Seed;
}
const char *Options::getLocusFilename() const
{
  return LocusFilename.c_str();
}

const char *Options::getGenotypesFilename() const
{
  return GenotypesFilename.c_str();
}
bool Options::CheckData()const{
    return checkData;
}
const char *Options::getAlleleFreqOutputFilename() const
{
  return AlleleFreqOutputFilename.c_str();
}
const char *Options::getErgodicAverageFilename() const
{
  return ErgodicAverageFilename.c_str();
}
bool Options::getFixedAlleleFreqs() const
{
  return (fixedallelefreqs);
}
const char *Options::getParameterFilename() const
{
  return ParameterFilename.c_str();
}

const char *Options::getRegressionOutputFilename() const
{
  return RegressionOutputFilename.c_str();
}
const vector<unsigned>& Options::getOutcomeVarColumns() const
{
  return OutcomeVarColumns;
}
const vector<unsigned>& Options::getCovariateColumns() const
{
  return CovariateColumns;
}
bool Options::getOutputAlleleFreq() const
{
  return OutputAlleleFreq;
}

int Options::getNumberOfOutcomes() const{
  return NumberOfOutcomes;
}
void Options::setNumberOfOutcomes(unsigned i){
  NumberOfOutcomes = i;
  if(i==0){
    useroptions.erase("outcomevarfile");
    useroptions.erase("covariatesfile");
    OutcomeVarFilename.clear();
    CovariatesFilename.clear();
    OutcomeVarColumns.clear();
    CovariateColumns.clear();
  }
}
void Options::setRegType(RegressionType R){
  RegType = R;
  if(R == None)
    setNumberOfOutcomes(0);
}
const char *Options::getCovariatesFilename() const
{
  return CovariatesFilename.c_str();
}

const char *Options::getPriorAlleleFreqFilename() const
{
  return PriorAlleleFreqFilename.c_str();
}
bool Options::getTestForResidualAllelicAssoc()const{
  return TestForResidualAllelicAssoc;
}

bool Options::getHWTestIndicator() const
{
  return HWTest;
}

const genepi::cvector<float>& Options::getrhoSamplerParams()const{
  return rhoSamplerParams;
}

double Options::getRegressionPriorPrecision()const{
  return regressionPriorPrecision;
}
const char *Options::getOutcomeVarFilename() const
{
  return OutcomeVarFilename.c_str();
}

const char *Options::getCoxOutcomeVarFilename() const
{
  return CoxOutcomeVarFilename.c_str();
}
bool Options::getTestForAllelicAssociation() const
{
  return TestForAllelicAssociation;
}

void Options::setTestForAllelicAssociation(bool b){
  TestForAllelicAssociation = b;
  if(!b)useroptions.erase("allelicassociationtest");
}
const char* Options::getEYFilename()const{
  return EYFilename.c_str();
}

bool Options::getDeleteOldResultsIndicator()const{
  return DeleteOldResultsIndicator;
}

void Options::DefineOptions(){
  //set up Option map

  addOption( "use-pedigree-for-individual", usePedForInd	 , false );
  addOption( "single-gamete-founder"	  , sglGameteFounder	 , true  );
  addOption( "warnings-are-errors"	  , warningsAreErrors	 , false );
  addOption( "ignore-invalid-parents"	  , ignoreInvalidParents , false );
  addOption( "outcome-is-binary"	  , outcomeIsBinary	 , false );
  addOption( "max-cpus"			  , maxCPUsToUse	 , 0	 );
  addOption( "max-pedigree-size"	  , maxPedigreeSize	 , 0	 );
  addOption( "male-x-homozyg-warn"	  , maleXHomozygWarn	 , false );
  addOption( "exclude-mendelian-errors"	  , excludeMendelError	 , true  );
  addOption( "exclude-unaffected-sibs"	  , excludeUnaffectedSibs, false );

  addFlag('h', "help");
  addFlag('o', "options");
  addFlag('v', "version");
  addFlag('c', "checkmode");
  addFlag('b', "printbuildinfo");

  addOption("samples", intOption, &TotalSamples);
  addOption("burnin", intOption, &burnin);
  addOption("every", intOption, &SampleEvery);
  addOption("locusfile", stringOption, &LocusFilename, true);
  addOption("genotypesfile", stringOption, &GenotypesFilename, true);
  addOption("priorallelefreqfile", stringOption, &PriorAlleleFreqFilename);
  //regression data files
  addOption("outcomevarfile", stringOption, &OutcomeVarFilename);
  addOption("coxoutcomevarfile", stringOption, &CoxOutcomeVarFilename);
  addOption("covariatesfile", stringOption, &CovariatesFilename);
  addOption("outcomevarcols", uivectorOption, &OutcomeVarColumns);
  addOption("covariatecols", uivectorOption, &CovariateColumns);
  //output files 
  addOption("logfile", outputfileOption, &LogFilename);
  addOption("paramfile", outputfileOption, &ParameterFilename);
  addOption("regparamfile", outputfileOption, &RegressionOutputFilename);
  addOption("allelefreqoutputfile", outputfileOption, &AlleleFreqOutputFilename);
  addOption("ergodicaveragefile", outputfileOption, &ErgodicAverageFilename);
  addOption("resultsdir", stringOption, &ResultsDir);

  addOption("regressionpriorprecision", doubleOption, &regressionPriorPrecision);
  addOption("fixedallelefreqs", boolOption, &fixedallelefreqs);
  // test options
  addOption("allelicassociationtest", boolOption, &TestForAllelicAssociation);
  addOption("residualldtest", boolOption, &TestForResidualAllelicAssoc);
  addOption("hwtest", boolOption, &HWTest);
  addOption("rhosamplerparams", rhoSamplerParams);

  // Other options
  addOption("numannealedruns", intOption, &NumAnnealedRuns);// number of coolnesses 
  addOption("displaylevel", intOption, &displayLevel);// output detail, 0 to 3
  addOption("seed", longOption, &Seed);// random number seed
  addOption("thermo", boolOption, &thermoIndicator);// Marginal likelihood by thermodynamic integration
  addOption("checkdata", boolOption, &checkData);// set to 0 to skip some data checks
  addOption("deleteoldresults", boolOption, &DeleteOldResultsIndicator);

  //null options - so argsfile can be used for input
  addOption("rparamfile", nullOption, 0);
}

bool Options::SetOptions(){

  if(!OptionReader::SetOptions())
    return false;

  //iterate over user options again, appending resultsdir to output filenames
  //this must be done here as resultsdir must be set first
  for(map<string, string>::iterator i = useroptions.begin(); i != useroptions.end(); ++i){
    if(ProgOptions[i->first].second == outputfileOption)
      *((string*)(ProgOptions[i->first].first)) = ResultsDir + "/" + i->second;
  }
  EYFilename = ResultsDir + "/" + EYFilename;
  //set indicators
  OutputAlleleFreq = (AlleleFreqOutputFilename.size()>0);
  checkData = checkData | getFlag("checkmode");

  return true;
}

bool Options::checkOptions(LogWriter& Log, int) {

  bool badOptions = false;//to indicate invalid options. Prog will exit at end of function if true.
  Log.setDisplayMode(Quiet);

  //check for burnin >= samples
  if(burnin >= TotalSamples){
    Log << "ERROR: 'samples' must be greater than 'burnin'.\n";
    badOptions = true;
  }
  //
  if(SampleEvery >= TotalSamples){
    Log << "ERROR: 'samples' must be greater than 'every'.\n";
    badOptions = true;
  }
  if(10*SampleEvery > (TotalSamples-burnin)){
    Log << "WARNING: 'every' should be less than ('samples' - 'burnin') / 10.\n"
        << "         Some output files may be empty.\n";
  }

  // **** Check whether genotypes file has been specified ****
  if ( GenotypesFilename.length() == 0 )
    {
      Log << "ERROR: Must specify genotypesfile.\n";
      badOptions = true;
    }
  // **** Check whether locus file has been specified ****
  if ( LocusFilename.length() == 0 )
    {
      Log << "ERROR: Must specify locusfile.\n";
      badOptions = true;
    }

  if( TestForAllelicAssociation ){
    if( NumberOfOutcomes < 1 ){
      Log << "WARNING: allelic association score test is not valid without a regression model."
	  << "         This option will be ignored.\n";
      setTestForAllelicAssociation(false);
    }
  }

  return badOptions;
}

///output Options table to file
void Options::PrintUserOptions(const char* filename){
  string ss = ResultsDir;
  ss.append("/");
  ss.append(filename);

  useroptions["rparamfile"] = PARAMFILE_ROBJECT;

  OptionReader::PrintUserOptions(ss.c_str());
}
