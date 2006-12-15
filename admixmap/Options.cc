/** 
 *   Options.cc 
 *   Base for classes to hold program options
 *   Copyright (c) 2006 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "Options.h"
#include "utils/StringConvertor.h"
#include <string.h>
#include <sstream>

using namespace std;

Options::Options()
{
  SetDefaultValues();
}
void Options::ReadUserOptions(int argc,  char** argv){
  if(argc == 2)
    ReadArgsFromFile(argv[1], useroptions);
  else{
    //NOTE: command line args will not be supported soon
    cout << "Warning: command-line arguments are deprecated" << endl;
    ReadCommandLineArgs(argc, argv);
  }
  //NOTE: derived classes must call SetOptions in their own constructors
  //  OptionMap ProgOptions;
  //SetOptions(ProgOptions);

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
  genotypesSexColumn = 0;
  TargetIndicator = 0;
  fixedallelefreqs = false;
  NumberOfOutcomes = -1;
  RegType = None;
  ScoreTestIndicator = false; //indicator for any of the score tests in ScoreTests class
  TestForAllelicAssociation = false;
  TestForResidualAllelicAssoc = false;
  OutputAlleleFreq = false;
  regressionPriorPrecision = 0.25;
  EYFilename = "ExpectedOutcomes.txt";

  useroptions["burnin"] = "100";
  useroptions["samples"] = "1100";
  useroptions["every"] = "5";
  useroptions["numannealedruns"] = "20";
  useroptions["displaylevel"] = "2";
  useroptions["logfile"] = "log.txt";
  useroptions["resultsdir"] = "results";
  useroptions["thermo"] = "0";
  useroptions["seed"] = "1";
  useroptions["checkdata"]="1";
  useroptions["targetindicator"] = "0";
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
unsigned int Options::getgenotypesSexColumn() const
{
  return genotypesSexColumn;
}

void Options::setgenotypesSexColumn(unsigned int i)
{
  genotypesSexColumn = i;
}
bool Options::CheckData()const{
    return checkData;
}
const char *Options::getAlleleFreqFilename() const
{
  return alleleFreqFilename.c_str();
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
  return (fixedallelefreqs || alleleFreqFilename.length());
}
const char *Options::getParameterFilename() const
{
  return ParameterFilename.c_str();
}

const char *Options::getRegressionOutputFilename() const
{
  return RegressionOutputFilename.c_str();
}
int Options::getTargetIndicator() const
{
  return TargetIndicator;
}

bool Options::getScoreTestIndicator() const
{
  return ScoreTestIndicator;
}
bool Options::getOutputAlleleFreq() const
{
  return OutputAlleleFreq;
}

int Options::getNumberOfOutcomes() const{
  return NumberOfOutcomes;
}
void Options::setNumberOfOutcomes(int i){
  NumberOfOutcomes = i;
  if(i==0){
    useroptions.erase("outcomevarfile");
    OutcomeVarFilename = "";
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
const char *Options::getAllelicAssociationScoreFilename() const
{
  return AllelicAssociationScoreFilename.c_str();
}
const char *Options::getPriorAlleleFreqFilename() const
{
  return PriorAlleleFreqFilename.c_str();
}
bool Options::getTestForResidualAllelicAssoc()const{
  return TestForResidualAllelicAssoc;
}
const char* Options::getResidualAllelicAssocScoreFilename()const{
  return ResidualAllelicAssocScoreFilename.c_str();
}
bool Options::getMHTest()const{
  return (bool)(MHTestFilename.size());
}
const char* Options::getMHTestFilename()const{
  return MHTestFilename.c_str();
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
  if(!b)useroptions.erase("allelicassociationscorefile");
}
const char* Options::getEYFilename()const{
  return EYFilename.c_str();
}

///////////////////////////////////////////////////////////////////////////////
int Options::ReadArgsFromFile(const char* filename, UserOptions& opt){
  ifstream fin(filename);
  if (0 == filename || 0 == strlen(filename)) return 1;
  if (!fin.is_open()) {
    string msg = "Cannot open file \"";
    msg += filename;
    msg += "\". Aborting";
    cerr << msg << endl;
    exit(1);
  } 

  std::string str;
  //read in line from file
  while (getline(fin,str,'\n')){// ## apparent memory leak 

    if( str.find_first_of("#") < str.length() ) str.erase( str.find_first_of("#") );//ignore #comments
    if(str.find_first_not_of(" \t\n\r") < str.length() )//skip blank lines. 
      {   
	str.erase(0, str.find_first_not_of(" \t\n\r") );//trim leading whitespace
	//trim remaining whitespace
	str.erase( str.find_last_not_of(" \t\n\r") + 1 );//trailing whitespace
	if( str.find_first_of(" \t\n\r") <= str.length() ){//check for any whitespace left
	  string::size_type eq = str.find("="), pos = str.find_first_of(" \t\n\r");
	  if( pos < eq )//check for space before '='
	    str.erase( pos, eq - pos );//remove space before '='
	  //str.erase( str.find_first_of(" \t\n\r"),str.find_last_of(" \t\n") - str.find_first_of(" \t\n\r") +1 );//after '='
	  eq = str.find("=");
	  pos = str.find_first_of(" \t\n\r", eq);//position of first space after the = 
	  str.erase( pos, str.find_first_not_of(" \t\n", pos) - pos );//remove space after '='
	}
	//add line to xargv
	string::size_type eq = str.find("=");
	opt[str.substr(0, eq)] = str.substr(eq+1 );
      }
    str.clear();
  }
  fin.close();
  return 0;
}

/// sets the value of a data member corresponding to a user option.
/// converts the value string to an appropriate type, with method indicated by opt.second
/// and assigns converted value to the data member with address opt.first
int Options::assign(OptionPair& opt, const string value){
  if(opt.second == "int")
    *((int*)opt.first) = atoi(value.c_str());
  else if(opt.second =="long")
    *((long*)opt.first) = atoi(value.c_str());
  else if(opt.second == "double")
    *((double*)opt.first) = atof(value.c_str());
  else if(opt.second =="float")
    *((float*)opt.first) = atof(value.c_str());
  else if(opt.second == "bool"){
    if (atoi(value.c_str()) ==1) *((bool*)opt.first) = true;
    else *((bool*)opt.first) = false;
  }
  else if(opt.second == "string")
    *((string*)opt.first) = value;

  else if(opt.second =="dvector"){
    StringConvertor::StringToVec(value, *((vector<double>*)opt.first));
  }
  else if(opt.second =="fvector"){
    StringConvertor::StringToVec(value, *((vector<float>*)opt.first));
  }
  else if(opt.second =="uivector"){
    StringConvertor::StringToVec(value, *((vector<unsigned>*)opt.first));
  }
  else if(opt.second =="old"){//deprecated option - return signal to erase
    return 2;
  }
  else if(opt.second != "null" && opt.second != "outputfile"){
    //skipping output file names as they are set later, prefixing resultsdir
    return 1;//unrecognised option
  }
  return 0;//success
}

void Options::SetOptions(OptionMap& ProgOptions)
{
  //set up Option map

  // Required options
  ProgOptions["samples"] = OptionPair(&TotalSamples, "int");
  ProgOptions["burnin"] = OptionPair(&burnin, "int");
  ProgOptions["every"] = OptionPair(&SampleEvery, "int");
  ProgOptions["locusfile"] = OptionPair(&LocusFilename, "string");
  ProgOptions["genotypesfile"] = OptionPair(&GenotypesFilename, "string");
  ProgOptions["priorallelefreqfile"] = OptionPair(&PriorAlleleFreqFilename, "string");
  ProgOptions["allelefreqfile"] = OptionPair(&alleleFreqFilename, "string");
  //regression data files
  ProgOptions["outcomevarfile"] = OptionPair(&OutcomeVarFilename, "string");
  ProgOptions["coxoutcomevarfile"] = OptionPair(&CoxOutcomeVarFilename, "string");
  ProgOptions["covariatesfile"] = OptionPair(&CovariatesFilename, "string");
  ProgOptions["outcomes"] = OptionPair(&NumberOfOutcomes, "int");
  ProgOptions["targetindicator"] = OptionPair(&TargetIndicator, "int");
  //standard output files (optional)
  ProgOptions["logfile"] = OptionPair(&LogFilename, "outputfile");
  ProgOptions["paramfile"] = OptionPair(&ParameterFilename, "outputfile");
  ProgOptions["regparamfile"] = OptionPair(&RegressionOutputFilename, "outputfile");
  ProgOptions["allelefreqoutputfile"] = OptionPair(&AlleleFreqOutputFilename, "outputfile");
  ProgOptions["ergodicaveragefile"] = OptionPair(&ErgodicAverageFilename, "outputfile");

  //optional results directory name option - default is 'results'
  ProgOptions["resultsdir"] = OptionPair(&ResultsDir, "string");
  ProgOptions["regressionpriorprecision"] = OptionPair(&regressionPriorPrecision, "double");
  ProgOptions["fixedallelefreqs"] = OptionPair(&fixedallelefreqs, "bool");
  // test options
  ProgOptions["allelicassociationscorefile"] = OptionPair(&AllelicAssociationScoreFilename, "outputfile");
  ProgOptions["residualallelicassocscorefile"] = OptionPair(&ResidualAllelicAssocScoreFilename, "outputfile");
  ProgOptions["mhtestfile"] = OptionPair(&MHTestFilename, "outputfile");

  // Other options
  ProgOptions["numannealedruns"] = OptionPair(&NumAnnealedRuns, "int");// number of coolnesses 
  ProgOptions["displaylevel"] = OptionPair(&displayLevel, "int");// output detail, 0 to 3
  ProgOptions["seed"] = OptionPair(&Seed, "long");// random number seed
  ProgOptions["thermo"] = OptionPair(&thermoIndicator, "bool");// Marginal likelihood by thermodynamic integration
  ProgOptions["checkdata"] = OptionPair(&checkData, "bool");// set to 0 to skip some data checks

  //parse user options
  bool badOptions = false;
  for(UserOptions::iterator i = useroptions.begin(); i != useroptions.end(); ++i){
    int status = assign(ProgOptions[i->first], i->second);
    if(status == 1){
      cerr << "Unknown option: " << i->first
	   << " with arg: " << i->second
	   << endl;
      badOptions = true;
    }
    else if (status==2) useroptions.erase(i->first);
  }
  if(badOptions)exit(1);

  //iterate over user options again, appending resultsdir to output filenames
  //this must be done here as resultsdir must be set first
  for(UserOptions::iterator i = useroptions.begin(); i != useroptions.end(); ++i){
    if(ProgOptions[i->first].second == "outputfile")
      *((string*)(ProgOptions[i->first].first)) = ResultsDir + "/" + i->second;
  }
  EYFilename = ResultsDir + "/" + EYFilename;//TODO: make a user option?
  //set indicators
  OutputAlleleFreq = (AlleleFreqOutputFilename.size()>0);
  TestForAllelicAssociation = (AllelicAssociationScoreFilename.size()>0);
  TestForResidualAllelicAssoc = (ResidualAllelicAssocScoreFilename.size()>0);
}

void Options::PrintOptions(){
  //set populations value in case it has changed
  //NB do similar for any option that can be changed outside Options
  std::ostringstream s;
  if (s << getPopulations()) // conversion worked
    {
    useroptions["populations"] = (char *)s.str().c_str();
    }
  //Now output Options table to args.txt
  string ss;
  ss = ResultsDir + "/args.txt";
  ofstream argstream(ss.c_str());

  for( UserOptions::iterator p= useroptions.begin(); p!=useroptions.end(); p++) {
    argstream << p->first << "=" << p->second <<endl;
  }
  argstream.close();
}

int Options::checkOptions(LogWriter &Log, int){
  bool badOptions = false;//to indicate invalid options. Prog will exit at end of function if true.
  Log.setDisplayMode(Quiet);


  //check for burnin >= samples
  if(burnin >= TotalSamples){
    Log << "ERROR: 'samples' must be greater than 'burnin'\n";
    badOptions = true;
  }
  //
  if(SampleEvery >= TotalSamples){
    Log << "ERROR: 'samples' must be greater than 'every'\n";
    badOptions = true;
  }
  if(10*SampleEvery > (TotalSamples-burnin)){
    Log << "WARNING: 'every' should be less than ('samples' - 'burnin') / 10. Some output files may be empty.\n";
  }

  if(badOptions) return 1;
  else return 0;
}

void Options::ReadCommandLineArgs(const int argc, char** argv){
  string name, value;
  vector<char*> args;
  char delims[] = "-=";
  for(int i = 1; i < argc; ++i){

    //tokenise argv, splitting on '-' and '='
    char *result = NULL;
    result = strtok( argv[i], delims );
    while( result != NULL ) {
      args.push_back(result);
      result = strtok( NULL, delims );
    }  
  }
  if(args.size() % 2){
    cerr << "ERROR: mismatched arguments" << endl;
    exit(1);
  }
  //TODO: allow spaces in vector args
  for( vector<char*>::iterator i = args.begin(); i != args.end(); ){
    name.assign(*i);
    ++i;
    value.assign(*i);
    ++i;
    useroptions[name] = value;
  }
}
