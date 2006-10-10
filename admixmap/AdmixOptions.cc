/** 
 *   ADMIXMAP
 *   AdmixOptions.cc 
 *   Class to hold program options
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AdmixOptions.h"
#include "utils/StringConvertor.h"
#include <string.h>
#include <sstream>
#include <numeric> // for checkInitAlpha
#include "Latent.h"// for POPADMIXSAMPLER

using namespace std;

AdmixOptions::AdmixOptions()
{
  Initialise();
}
AdmixOptions::AdmixOptions(int argc,  char** argv){
  Initialise();
  if(argc == 2)
    ReadArgsFromFile(argv[1], useroptions);
  else{
    //NOTE: command line args will not be supported soon
    cout << "Warning: command-line arguments are deprecated" << endl;
    ReadCommandLineArgs(argc, argv);
  }
  SetOptions();

}

void AdmixOptions::Initialise(){
  // global variables store option values: variable names not necessarily same as 
  // command-line option names which are all lower-case 
  burnin = 100;
  TotalSamples = 1100;
  SampleEvery = 10;
  Seed = 1;
  TargetIndicator = 0;
  Populations = 1;

  displayLevel = 2; 
  OutputFST = false;
  genotypesSexColumn = 0;
  locusForTestIndicator = false;
  LocusForTest = -1;
  fixedallelefreqs = false;
  correlatedallelefreqs = false;
  RandomMatingModel = false;
  NumberOfOutcomes = -1;
  RegType = None;
  GlobalRho = true;//corresponds to globalrho = 1;
  PopAdmixPropsAreEqual = false;
  IndAdmixHierIndicator = true; //hierarchical model on ind admixture
  HapMixModelIndicator = false; //model haplotypes with mixture model
  chibIndicator = false;//calculate marginal likelihood by Chib method
  thermoIndicator = false; // calculate marginal likelihood by thermodynamic method
  TestOneIndivIndicator = false; // evaluate marginal likelihood for single individual
  NumAnnealedRuns = 20; // default if option thermo not specified 
  ScoreTestIndicator = false; //indicator for any of the score tests in ScoreTests class
  TestForAdmixtureAssociation = false;
  StratificationTestIndicator = false;
  TestForAffectedsOnly = false;
  TestForAllelicAssociation = false;
  TestForHaplotypeAssociation = false;
  TestForDispersion = false;
  TestForLinkageWithAncestry = false;
  TestForResidualAllelicAssoc = false;
  TestForMisspecifiedAlleleFreqs = false;
  TestForMisspecifiedAlleleFreqs2 = false;
  HWTest = false;
  OutputAlleleFreq = false;
  checkData = true;


//   globalrhoPrior.push_back(3.0);//rhoalpha 
//   globalrhoPrior.push_back(0.5);//rhobeta

//   rhoPrior.push_back(6.0);//rhoalpha 
//   rhoPrior.push_back(5.0);//rhobeta shape
//   rhoPrior.push_back(4.0);//rhobeta rate

  initalpha.resize(2);//TODO: check if this is necessary
  //gamma(3, 0.01) prior on dispersion parameter
  etamean = 100.0; 
  etavar = 2500.0; 

  regressionPriorPrecision = 0.25;
  ResultsDir = "results";
  LogFilename = "log.txt";

  //prior on frequency Dirichlet prior params in hapmixmodel
  //allelefreqprior.push_back(0.01);//shape
  //allelefreqprior.push_back(0.1);//prior shape of rate
  //allelefreqprior.push_back(1.0);//prior rate of rate

  LikRatioFilename = "LikRatioFile.txt";//hardcoding for now, can change later
  EYFilename = "ExpectedOutcomes.txt";

  // option names and default option values are stored as strings in a map container 
  // these are default values
  // other specified options will be appended to this array 
  useroptions["burnin"] = "100";
  useroptions["samples"] = "1100";
  useroptions["every"] = "5";
  useroptions["numannealedruns"] = "20";
  useroptions["targetindicator"] = "0";
  useroptions["displaylevel"] = "2";
  useroptions["fixedallelefreqs"] = "0";
  useroptions["correlatedallelefreqs"] = "0";
  useroptions["logfile"] = "log.txt";
  useroptions["resultsdir"] = "results";
  useroptions["randommatingmodel"] = "0";
  useroptions["globalrho"] = "1";
  useroptions["indadmixhiermodel"] = "1";
  useroptions["hapmixmodel"] = "0";
  useroptions["chib"] = "0";
  useroptions["thermo"] = "0";
  //global rho: default gamma (3, 0.5) prior has mean 6, variance 12 
  useroptions["globalsumintensitiesprior"] = "3.0,0.5";
  // non-global rho: default gamma-gamma prior with parameters n=6, alpha=5, beta=4
  // effective prior mean is 6*4/(5-1) = 6 and effective prior variance is 6*7 / (5-2) = 14
  useroptions["sumintensitiesprior"] = "6.0,5.0,4.0";
  useroptions["hapmixrhopriormeans"] = "8, 8, 1.4";
  useroptions["hapmixrhopriorvars"] = "1000, 1000, 1000";
  useroptions["seed"] = "1";
  useroptions["regressionpriorprecision"] = "0.25";
}

AdmixOptions::~AdmixOptions()
{
}

// each option has a function to return its value
const string AdmixOptions::getResultsDir() const{
  return ResultsDir;
}
const char *AdmixOptions::getAlleleFreqFilename() const
{
  return alleleFreqFilename.c_str();
}

long AdmixOptions::getBurnIn() const
{
  return burnin;
}

const char *AdmixOptions::getStratTestFilename() const
{
  return StratTestFilename.c_str();
}

const char *AdmixOptions::getErgodicAverageFilename() const
{
  return ErgodicAverageFilename.c_str();
}

const char *AdmixOptions::getFSTOutputFilename() const
{
  return FSTOutputFilename.c_str();
}

bool AdmixOptions::getOutputFST() const
{
  return OutputFST;
}

bool AdmixOptions::getFixedAlleleFreqs() const
{
  return (fixedallelefreqs || alleleFreqFilename.length());
}
bool AdmixOptions::getCorrelatedAlleleFreqs() const
{
  return correlatedallelefreqs;
}

const char *AdmixOptions::getLocusFilename() const
{
  return LocusFilename.c_str();
}

const char *AdmixOptions::getGenotypesFilename() const
{
  return GenotypesFilename.c_str();
}

const char *AdmixOptions::getIndAdmixtureFilename() const
{
  return IndAdmixtureFilename.c_str();
}

const char *AdmixOptions::getParameterFilename() const
{
  return ParameterFilename.c_str();
}

const char *AdmixOptions::getRegressionOutputFilename() const
{
  return RegressionOutputFilename.c_str();
}

const char *AdmixOptions::getEtaOutputFilename() const
{
  return EtaOutputFilename.c_str();
}

const char *AdmixOptions::getReportedAncestryFilename() const
{
  return ReportedAncestryFilename.c_str();
}

long AdmixOptions::getTotalSamples() const
{
  return TotalSamples;
}

int AdmixOptions::getTargetIndicator() const
{
  return TargetIndicator;
}
int AdmixOptions::getDisplayLevel() const
{
  return displayLevel;
}

bool AdmixOptions::getScoreTestIndicator() const
{
  return ScoreTestIndicator;
}
const char *AdmixOptions::getAffectedsOnlyScoreFilename() const
{
  return AffectedsOnlyScoreFilename.c_str();
}

const char *AdmixOptions::getAlleleFreqScoreFilename() const
{
  return AlleleFreqScoreFilename.c_str();
}

const char *AdmixOptions::getAlleleFreqScoreFilename2() const
{
  return AlleleFreqScoreFilename2.c_str();
}

const char *AdmixOptions::getAlleleFreqOutputFilename() const
{
  return AlleleFreqOutputFilename.c_str();
}
const char *AdmixOptions::getAlleleFreqPriorOutputFilename() const
{
  return AlleleFreqPriorOutputFilename.c_str();
}

bool AdmixOptions::getOutputAlleleFreq() const
{
  return OutputAlleleFreq;
}

int AdmixOptions::getNumberOfOutcomes() const{
  return NumberOfOutcomes;
}
void AdmixOptions::setNumberOfOutcomes(int i){
  NumberOfOutcomes = i;
  if(i==0){
    useroptions.erase("outcomevarfile");
    OutcomeVarFilename = "";
  }
}
void AdmixOptions::setRegType(RegressionType R){
  RegType = R;
  if(R == None)
    setNumberOfOutcomes(0);
}
const char *AdmixOptions::getAssocScoreFilename() const
{
  return AssocScoreFilename.c_str();
}

long AdmixOptions::getSampleEvery() const
{
  return SampleEvery;
}

const char *AdmixOptions::getDispersionTestFilename() const
{
  return DispersionTestFilename.c_str();
}

bool AdmixOptions::getHWTestIndicator() const
{
  return HWTest;
}

const char *AdmixOptions::getHWTestFilename() const
{
  return HWTestFilename.c_str();
}
const char *AdmixOptions::getHistoricalAlleleFreqFilename() const
{
  return HistoricalAlleleFreqFilename.c_str();
}

const char *AdmixOptions::getCovariatesFilename() const
{
  return CovariatesFilename.c_str();
}

bool AdmixOptions::isRandomMatingModel() const
{
  return RandomMatingModel;
}

bool AdmixOptions::getIndAdmixHierIndicator() const{
  return IndAdmixHierIndicator;
}
bool AdmixOptions::getHapMixModelIndicator() const{
  return HapMixModelIndicator;
}
bool AdmixOptions::getChibIndicator()const{
  return chibIndicator;
}
bool AdmixOptions::getThermoIndicator()const{
  return thermoIndicator;
}
bool AdmixOptions::getTestOneIndivIndicator()const{
  return TestOneIndivIndicator;
}

long AdmixOptions::getNumAnnealedRuns()const{
  return NumAnnealedRuns;
}

bool AdmixOptions::isGlobalRho() const
{
  return GlobalRho;
}

bool AdmixOptions::getLocusForTestIndicator() const
{
  return locusForTestIndicator;
}

int AdmixOptions::getLocusForTest() const
{
  return LocusForTest;
}

const char *AdmixOptions::getAllelicAssociationScoreFilename() const
{
  return AllelicAssociationScoreFilename.c_str();
}

const char *AdmixOptions::getAncestryAssociationScoreFilename() const
{
  return AncestryAssociationScoreFilename.c_str();
}

const char *AdmixOptions::getLogFilename() const
{
  return LogFilename.c_str();
}

int AdmixOptions::getPopulations() const
{
  return Populations;
}

void AdmixOptions::setPopulations(int num)
{
  Populations = num;
}
const char *AdmixOptions::getPriorAlleleFreqFilename() const
{
  return PriorAlleleFreqFilename.c_str();
}

bool AdmixOptions::getTestForAdmixtureAssociation() const
{
  return TestForAdmixtureAssociation;
}

bool AdmixOptions::getTestForResidualAllelicAssoc()const{
  return TestForResidualAllelicAssoc;
}
const char* AdmixOptions::getResidualAllelicAssocScoreFilename()const{
  return ResidualAllelicAssocScoreFilename.c_str();
}
double AdmixOptions::getRhoPriorMean()const{
  if( !HapMixModelIndicator && (GlobalRho || !IndAdmixHierIndicator ) )
    return globalrhoPrior[0] / globalrhoPrior[1];
  else 
    return rhoPrior[0] * rhoPrior[2] / (rhoPrior[1] - 1.0);
}
long AdmixOptions::getSeed() const
{
  return Seed;
}

double AdmixOptions::getRhoalpha() const {
  if(!HapMixModelIndicator && (GlobalRho || !IndAdmixHierIndicator)) {
    return globalrhoPrior[0];
  } else {
    return rhoPrior[0];
  }
}

double AdmixOptions::getRhobeta() const
{
  return globalrhoPrior[1];
}
double AdmixOptions::getRhobetaShape()const{
  return rhoPrior[1];
}
double AdmixOptions::getRhobetaRate()const{
  return rhoPrior[2];
}
const std::vector<double> &AdmixOptions::getHapMixRhoPriorMeans()const{
  return hapmixrhopriormeans;
}
const std::vector<double> &AdmixOptions::getHapMixRhoPriorVars()const{
  return hapmixrhopriorvars;
}
double AdmixOptions::getEtaMean() const{
  return etamean;
}
double AdmixOptions::getEtaVar() const{
  return etavar;
}
double AdmixOptions::getRegressionPriorPrecision()const{
  return regressionPriorPrecision;
}

bool AdmixOptions::getStratificationTest() const
{
  return StratificationTestIndicator;
}

void AdmixOptions::setStratificationTest(bool b){
  StratificationTestIndicator = b;
}
const char *AdmixOptions::getOutcomeVarFilename() const
{
  return OutcomeVarFilename.c_str();
}

const char *AdmixOptions::getCoxOutcomeVarFilename() const
{
  return CoxOutcomeVarFilename.c_str();
}

bool AdmixOptions::getTestForAffectedsOnly() const
{
  return TestForAffectedsOnly;
}

void AdmixOptions::setTestForAffectedsOnly(bool b){
  TestForAffectedsOnly = b;
  if(b && AffectedsOnlyScoreFilename.length()==0){
    //set default filename
  }
  else
    useroptions.erase("affectedsonlyscorefile");
}

bool AdmixOptions::getTestForAllelicAssociation() const
{
  return TestForAllelicAssociation;
}

void AdmixOptions::setTestForAllelicAssociation(bool b){
  TestForAllelicAssociation = b;
  if(!b)useroptions.erase("allelicassociationscorefile");
}

bool AdmixOptions::getTestForDispersion() const
{
  return TestForDispersion;
}

bool AdmixOptions::getTestForLinkageWithAncestry() const
{
  return TestForLinkageWithAncestry;
}

void AdmixOptions::setTestForLinkageWithAncestry(bool b){
  TestForLinkageWithAncestry = b;
  if(b && AncestryAssociationScoreFilename.length()==0){
    //set default filename
  }
  else
    useroptions.erase("ancestryassociationscorefile");
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs() const
{
  return TestForMisspecifiedAlleleFreqs;
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs2() const
{
  return TestForMisspecifiedAlleleFreqs2;
}

const char *AdmixOptions::getHaplotypeAssociationScoreFilename() const
{
  return HaplotypeAssociationScoreFilename.c_str();
}

bool AdmixOptions::getTestForHaplotypeAssociation() const
{
  return TestForHaplotypeAssociation;
}

void AdmixOptions::setTestForHaplotypeAssociation(bool b){
  TestForHaplotypeAssociation = b;
  if(!b)useroptions.erase("haplotypeassociationscorefile");
}

const char *AdmixOptions::getEtaPriorFilename() const
{
  return EtaPriorFilename.c_str();
}

const char* AdmixOptions::getEYFilename()const{
  return EYFilename.c_str();
}

int AdmixOptions::sizeInitAlpha() const
{
  //unsigned size = 0;
  //if(alpha0.size() > 0)++size;
  //if(alpha1.size() > 0)++size;
  return initalpha.size();
}

std::vector<double> AdmixOptions::getInitAlpha(int gamete) const
{
//   switch(gamete){
//   case 0: return alpha0;
//   case 1: return alpha1;
//   default : 
//     {
//       cerr<<"Error in call to getInitAlpha, g > 1"<<endl;
//       exit(1);
//     }
//   }
  return initalpha[gamete];
}
std::vector<std::vector<double> > AdmixOptions::getInitAlpha()const{
  return initalpha;
}

unsigned int AdmixOptions::getgenotypesSexColumn() const
{
  return genotypesSexColumn;
}

void AdmixOptions::setgenotypesSexColumn(unsigned int i)
{
  genotypesSexColumn = i;
}

bool AdmixOptions::isSymmetric()const{
  return _symmetric;
}
bool AdmixOptions::isAdmixed(unsigned gamete)const{
  return _admixed[gamete];
}

const char*AdmixOptions::getLikRatioFilename() const{
  return LikRatioFilename.c_str();
}
const char* AdmixOptions::getIndAdmixModeFilename()const{
  return IndAdmixModeFilename.c_str();
}
bool AdmixOptions::CheckData()const{
    return checkData;
}
bool AdmixOptions::PopAdmixturePropsAreEqual()const{
  return PopAdmixPropsAreEqual;
}
const vector<float>& AdmixOptions::getrhoSamplerParams()const{
  return rhoSamplerParams;
}
const vector<float>& AdmixOptions::getrhoPriorParamSamplerParams() const{
  return rhoPriorParamSamplerParams;
}

const std::vector<double> & AdmixOptions::getAlleleFreqPriorParams()const{
  return allelefreqprior;
}

int AdmixOptions::ReadArgsFromFile(const char* filename, UserOptions& opt){
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
int AdmixOptions::assign(OptionPair& opt, const string value){
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
  else if(opt.second =="old"){//deprecated option - return signal to erase
    return 2;
  }
  else if(opt.second != "null" && opt.second != "outputfile"){
    //skipping output file names as they are set later, prefixing resultsdir
    return 1;//unrecognised option
  }
  return 0;//success
}

void AdmixOptions::SetOptions()
{
  //set up Option map
  OptionMap Options;
  // Required options
  Options["samples"] = OptionPair(&TotalSamples, "int");
  Options["burnin"] = OptionPair(&burnin, "int");
  Options["every"] = OptionPair(&SampleEvery, "int");
  Options["locusfile"] = OptionPair(&LocusFilename, "string");
  Options["genotypesfile"] = OptionPair(&GenotypesFilename, "string");
  // Must specify one of the following
  Options["populations"] = OptionPair(&Populations, "int");
  Options["priorallelefreqfile"] = OptionPair(&PriorAlleleFreqFilename, "string");
  Options["allelefreqfile"] = OptionPair(&alleleFreqFilename, "string");
  Options["historicallelefreqfile"] = OptionPair(&HistoricalAlleleFreqFilename, "string");
  //regression data files
  Options["outcomevarfile"] = OptionPair(&OutcomeVarFilename, "string");
  Options["coxoutcomevarfile"] = OptionPair(&CoxOutcomeVarFilename, "string");
  Options["covariatesfile"] = OptionPair(&CovariatesFilename, "string");
  Options["outcomes"] = OptionPair(&NumberOfOutcomes, "int");
  Options["targetindicator"] = OptionPair(&TargetIndicator, "int");
  Options["reportedancestry"] = OptionPair(&ReportedAncestryFilename, "string");
  //standard output files (optional)
  Options["paramfile"] = OptionPair(&ParameterFilename, "outputfile");
  Options["regparamfile"] = OptionPair(&RegressionOutputFilename, "outputfile");
  Options["dispparamfile"] = OptionPair(&EtaOutputFilename, "outputfile");
  Options["logfile"] = OptionPair(&LogFilename, "outputfile");
  Options["indadmixturefile"] = OptionPair(&IndAdmixtureFilename, "outputfile");
  Options["allelefreqoutputfile"] = OptionPair(&AlleleFreqOutputFilename, "outputfile");
  Options["allelefreqprioroutputfile"] = OptionPair(&AlleleFreqPriorOutputFilename, "outputfile");
  Options["ergodicaveragefile"] = OptionPair(&ErgodicAverageFilename, "outputfile");
  //optional results directory name option - default is 'results'
  Options["resultsdir"] = OptionPair(&ResultsDir, "string");
  //prior and model specification
  Options["randommatingmodel"] = OptionPair(&RandomMatingModel, "bool");
  Options["globalrho"] = OptionPair(&GlobalRho, "bool");
  Options["indadmixhiermodel"] = OptionPair(&IndAdmixHierIndicator, "bool");
  Options["hapmixmodel"] = OptionPair(&HapMixModelIndicator, "bool");
  Options["etapriorfile"] = OptionPair(&EtaPriorFilename, "string");
  Options["globalsumintensitiesprior"] = OptionPair(&globalrhoPrior, "dvector");
  Options["sumintensitiesprior"] = OptionPair(&rhoPrior, "dvector");
  Options["hapmixrhopriormeans"] = OptionPair(&hapmixrhopriormeans, "dvector");
  Options["hapmixrhopriorvars"] = OptionPair(&hapmixrhopriorvars, "dvector");
  Options["allelefreqprior"] = OptionPair(&allelefreqprior, "dvector");
  Options["etapriormean"] = OptionPair(&etamean, "double");
  Options["etapriorvar"] = OptionPair(&etavar, "double");
  Options["regressionpriorprecision"] = OptionPair(&regressionPriorPrecision, "double");
  Options["admixtureprior"] = OptionPair(&initalpha[0], "dvector");
  Options["admixtureprior1"] = OptionPair(&initalpha[1], "dvector");
  Options["fixedallelefreqs"] = OptionPair(&fixedallelefreqs, "bool");
  Options["correlatedallelefreqs"] = OptionPair(&correlatedallelefreqs, "bool");
  Options["popadmixproportionsequal"] = OptionPair(&PopAdmixPropsAreEqual, "bool");
  //sampler settings
  Options["rhosamplerparams"] = OptionPair(&rhoSamplerParams, "fvector");
  Options["rhopriorsamplerparams"] = OptionPair(&rhoPriorParamSamplerParams, "fvector");
  // test options
  Options["allelicassociationscorefile"] = OptionPair(&AllelicAssociationScoreFilename, "outputfile");
  Options["residualallelicassocscorefile"] = OptionPair(&ResidualAllelicAssocScoreFilename, "outputfile");
  Options["ancestryassociationscorefile"] = OptionPair(&AncestryAssociationScoreFilename, "outputfile");
  Options["affectedsonlyscorefile"] = OptionPair(&AffectedsOnlyScoreFilename, "outputfile");
  Options["admixturescorefile"] = OptionPair(&AssocScoreFilename, "outputfile");
  Options["haplotypeassociationscorefile"] = OptionPair(&HaplotypeAssociationScoreFilename, "outputfile");
  Options["stratificationtestfile"] = OptionPair(&StratTestFilename, "outputfile");
  Options["allelefreqscorefile"] = OptionPair(&AlleleFreqScoreFilename, "outputfile");
  Options["allelefreqscorefile2"] = OptionPair(&AlleleFreqScoreFilename2, "outputfile");
  Options["dispersiontestfile"] = OptionPair(&DispersionTestFilename, "outputfile");
  Options["fstoutputfile"] = OptionPair(&FSTOutputFilename, "outputfile");
  Options["hwscoretestfile"] = OptionPair(&HWTestFilename, "outputfile");
  Options["likratiofile"] = OptionPair(&LikRatioFilename, "outputfile");
  Options["indadmixmodefile"] = OptionPair(&IndAdmixModeFilename, "outputfile");
  Options["testgenotypesfile"] = OptionPair(0, "null");
  Options["locusfortest"] = OptionPair(&LocusForTest, "int");
  // Other options
  Options["numannealedruns"] = OptionPair(&NumAnnealedRuns, "int");// number of coolnesses 
  Options["displaylevel"] = OptionPair(&displayLevel, "int");// output detail, 0 to 3
  Options["seed"] = OptionPair(&Seed, "long");// random number seed
  Options["chib"] = OptionPair(&chibIndicator, "bool");// Marginal likelihood by Chib algo
  Options["thermo"] = OptionPair(&thermoIndicator, "bool");// Marginal likelihood by thermodynamic integration
  Options["testoneindiv"] = OptionPair(&TestOneIndivIndicator, "bool");// ML for one individual in a collection 
  Options["checkdata"] = OptionPair(&checkData, "bool");// set to 0 to skip some data checks
  //old options - do nothing but kept for backward-compatibility with old scripts
  Options["analysistypeindicator"] = OptionPair(0, "old");
  Options["coutindicator"] = OptionPair(0, "old");
  Options["truncationpoint"] = OptionPair(0, "old");

  //parse user options
  bool badOptions = false;
  for(UserOptions::iterator i = useroptions.begin(); i != useroptions.end(); ++i){
    int status = assign(Options[i->first], i->second);
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
    if(Options[i->first].second == "outputfile")
      *((string*)(Options[i->first].first)) = ResultsDir + "/" + i->second;
  }
  EYFilename = ResultsDir + "/" + EYFilename;//TODO: make a user option?

  //set indicators
  OutputAlleleFreq = (AlleleFreqOutputFilename.size()>0);
  TestForAllelicAssociation = (AllelicAssociationScoreFilename.size()>0);
  TestForResidualAllelicAssoc = (ResidualAllelicAssocScoreFilename.size()>0);
  TestForLinkageWithAncestry = (AncestryAssociationScoreFilename.size()>0);
  TestForAffectedsOnly = (AffectedsOnlyScoreFilename.size()>0);
  TestForAdmixtureAssociation = (AssocScoreFilename.size()>0);
  TestForHaplotypeAssociation = (HaplotypeAssociationScoreFilename.size()>0);
  StratificationTestIndicator = (StratTestFilename.size()>0);
  TestForMisspecifiedAlleleFreqs = (AlleleFreqScoreFilename.size()>0);
  TestForMisspecifiedAlleleFreqs2 = (AlleleFreqScoreFilename2.size()>0);
  TestForDispersion = (DispersionTestFilename.size()>0);
  OutputFST = (FSTOutputFilename.size()>0);
  HWTest = (HWTestFilename.size()>0);
  locusForTestIndicator = (LocusForTest>-1);

}

void AdmixOptions::PrintOptions(){
  //set populations value in case it has changed
  //NB do similar for any option that can be changed outside AdmixOptions
  std::ostringstream s;
  if (s << Populations) // conversion worked
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

int AdmixOptions::checkOptions(LogWriter &Log, int NumberOfIndividuals){
  bool badOptions = false;//to indicate invalid options. Prog will exit at end of function if true.
  Log.setDisplayMode(Quiet);

  // **** analysis type  ****
  if(CoxOutcomeVarFilename.length() ){
    Log << "Cox Regression\n";
    if(NumberOfOutcomes>-1)++NumberOfOutcomes;
    else NumberOfOutcomes = 1;
    if(RegType == None)RegType = Cox;
    else RegType = Multiple;
  }

  if (NumberOfIndividuals==1) {
    IndAdmixHierIndicator = false;
    useroptions["indadmixhiermodel"]="0";
    Log << "One individual analysis";
  } else if (RegType == None) { //no regression
    NumberOfOutcomes = 0;
    if(AffectedsOnlyScoreFilename.length()>0) {
      Log << "Affecteds-only analysis";
    } else {
      Log << "Cross sectional analysis, no outcome";
    }
  } else if (RegType == Linear) {
    Log << "Cross sectional analysis, continuous outcome";
  } else if (RegType == Logistic) {
    Log << "Case-control or cross sectional analysis, binary outcome";
  } else{ //if (RegType == Both) {
    Log << "Cross sectional analysis, multiple outcome";
  }
  
  if(chibIndicator){
    Log << " with marginal likelihood calculation ";
    if(NumberOfIndividuals >1 )Log << "for first individual";
  }
  Log << "\n";

  if(TestOneIndivIndicator && RegType != None){
    badOptions = true;
    Log << "ERROR: cannot have a test individual with a regression model\n";
  }

  if(OutcomeVarFilename.length() == 0 && CoxOutcomeVarFilename.length()==0){
    if(NumberOfOutcomes > 0){
      Log.setDisplayMode(On);
      Log << "ERROR: 'outcomes' > 0 and no outcomevarfile specified\n";
      badOptions = true;
    }
    if(CovariatesFilename.length()){
      Log << "ERROR: covariatesfile specified without outcomevarfile\n";
    }
    //should check for specified targetindicator too, simply ignoring for now
    if(RegressionOutputFilename.length() > 0){
      Log << "ERROR: regparamfile option is not valid without a regression model\n"
	  << "\tThis option will be ignored\n";
      RegressionOutputFilename = "";
      useroptions.erase("regparamfile");
    }
  }
  
  // **** Hierarchical model on ind admixture ****
  if (!IndAdmixHierIndicator)
    {
      Log << "No hierarchical model for individual admixture or sum-intensities.\n";
      
      if(ParameterFilename.length() > 0 ){
	Log << "ERROR: paramfile option is not valid with indadmixhierindicator = 0\n"
	    << "\tThis option will be ignored\n";
	ParameterFilename = "";
	useroptions.erase("paramfile");
      }

      if(EtaOutputFilename.length() > 0 ){
	Log << "ERROR: dispparamfile option is not valid with indadmixhierindicator = 0\n"
	    << "\tThis option will be ignored\n";
	EtaOutputFilename = "";
	useroptions.erase("dispparamfile");
      }
      GlobalRho = false;
      useroptions["globalrho"] = "0";
    }
  
  // **** sumintensities ****
  if(HapMixModelIndicator){
    Log << "Haplotype mixture model with " << Populations << " block state";if(Populations>1)Log << "s"; Log << "\n";
    RandomMatingModel = false;useroptions["randommatingmodel"] = "0";
    GlobalRho = false; useroptions["globalrho"] = "0";
  }
  else {
    // **** Mating Model **** 
    if(RandomMatingModel )
      Log << "Model assuming random mating.\n";
    else 
      Log << "Model assuming assortative mating.\n";
    
    if( GlobalRho ) {
      Log << "Model with global sum-intensities\n";
      if(globalrhoPrior.size() != 2) {
	Log.setDisplayMode(On);
	Log << "ERROR: globalsumintensitiesprior must have length 2\n";
	badOptions = true;
	if(globalrhoPrior[0] <= 0.0 || globalrhoPrior[1] <= 0.0) {
	  Log.setDisplayMode(On);
	  Log << "ERROR: all elements of globalsumintensitiesprior must be > 0\n";
	  badOptions = true;
	}  
      }
    } else { // sumintensities at individual or gamete level
      if( RandomMatingModel )
	Log << "Model with gamete specific sum-intensities.\n";
      else
	Log << "Model with individual-specific sum-intensities.\n";
      
      if( (rhoPrior.size() != 3) && IndAdmixHierIndicator ) {
	Log.setDisplayMode(On);
	Log << "ERROR: for hierarchical model, sumintensitiesprior must have length 3\n";
	badOptions = true;
      }
      if(rhoPrior[0] <= 0.0 || rhoPrior[1] <= 0.0 || rhoPrior[2]<=0.0) {
	Log.setDisplayMode(On);
	Log << "ERROR: all elements of sumintensitiesprior must be > 0\n";
	badOptions = true;
      }
    }
  }
  if(HapMixModelIndicator){
    //TODO:check length of rhoprior vectors

    if(useroptions["allelefreqprior"].size() && (allelefreqprior.size() !=3)) {
      Log << "Error: 'allelefreqprior' must have length 3\n";
      badOptions = true;
    }

    useroptions.erase("sumintensitiesprior") ;
    useroptions.erase("globalsumintensitiesprior") ;  
  }
  else{
    useroptions.erase("hapmixrhopriormeans") ;
    useroptions.erase("hapmixrhopriorvars") ;
    if((GlobalRho || !IndAdmixHierIndicator ) ) {
      Log << "Gamma prior on sum-intensities with shape parameter: " << globalrhoPrior[0] << "\n"
	  << "and rate (1 / location) parameter " << globalrhoPrior[1] << "\n";
      Log << "Effective prior mean of sum-intensities is " << globalrhoPrior[0] / globalrhoPrior[1] << "\n";
      Log << "Effective prior variance of sum-intensities is " 
	  << globalrhoPrior[0] / (globalrhoPrior[1]*globalrhoPrior[1]) << "\n";
      useroptions.erase("sumintensitiesprior") ;
    } else {  
      double rhopriormean = rhoPrior[0] * rhoPrior[2] / (rhoPrior[1] - 1.0);
      Log << "Population prior distribution of sum-intensities specified as Gamma with shape parameter "
	  << rhoPrior[0] << "\n"
	  << "and Gamma prior on rate (1 / location) parameter with shape and rate parameters: "
	  << rhoPrior[1] << " & "
	  << rhoPrior[2] << "\n"
	  << "Effective prior mean of sum-intensities is ";
      if(rhoPrior[1]>1)Log << rhopriormean << "\n";else Log << "undefined\n";
      Log << "Effective prior variance of sum-intensities is ";
      if(rhoPrior[1]>2)
	Log << rhopriormean * (rhopriormean + 1.0) / (rhoPrior[1] - 2) << "\n";
      else Log <<"undefined\n";
      useroptions.erase("globalsumintensitiesprior") ;  
    }
  //if(allelefreqprior.size())
  //Log << "Warning: option 'allelefreqprior' is valid only with a hapmixmodel. This option will be ignored\n";

  }

  //Prior on admixture
  setInitAlpha(Log);
  if(Populations > 1 && IndAdmixHierIndicator && !HapMixModelIndicator){
    Log << "Gamma(1, 1) prior on population admixture Dirichlet parameters.\n";
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

  // **** model for allele freqs ****

  //fixed allele freqs
  if( alleleFreqFilename.length() ||
           (PriorAlleleFreqFilename.length() && fixedallelefreqs ) ){
    Log << "Analysis with fixed allele frequencies.\n";
    if(OutputAlleleFreq){
      Log << "ERROR: allelefreqoutputfile option is invalid with fixed allele frequencies\n"
	  << "       this option will be ignored\n";
      useroptions.erase("allelefreqoutputfile");
      OutputAlleleFreq = false;
    }
  }
  //prior allele freqs
  else if( PriorAlleleFreqFilename.length() && !fixedallelefreqs ){
    Log << "Analysis with prior allele frequencies.\n";
    if(correlatedallelefreqs) {
      Log << "Analysis with correlated allele frequencies\n";
    }
  }
  //historic allele freqs
  else if( HistoricalAlleleFreqFilename.length() > 0 ){
    Log << "Analysis with dispersion model for allele frequencies.\n";
  }
  //default priors ('populations' option)
  else if(Populations > 0 )
    {
      Log << "No allelefreqfile, priorallelefreqfile or historicallelefreqfile supplied;\n"
	  << "Default priors will be set for the allele frequencies";
      if(!HapMixModelIndicator)
	Log << " with " << Populations << " population(s)";
      Log << "\n";
      if(correlatedallelefreqs) {
	Log << "Analysis with correlated allele frequencies\n";
      }
    }
  
  if( (FSTOutputFilename.length() > 0) && (HistoricalAlleleFreqFilename.length() == 0) ){
    Log << "ERROR: fstoutputfile option is only valid with historicallelefreqfile option\n"
	<< "       this option will be ignored\n";
    OutputFST = false;
    FSTOutputFilename = "";
    useroptions.erase("fstoutputfile");
  }

  // **** score tests ****
  if( TestForLinkageWithAncestry && Populations == 1 ){
    Log << "Cannot test for linkage with ancestry with 1 population.\n";
    badOptions = true;
  }
  if(TestForAdmixtureAssociation &&
      ( TestForLinkageWithAncestry || TestForAllelicAssociation ) ) {
    Log << "Cannot test for linkage with ancestry or allelic association\n"
	<< "with score test for association. Can only use affecteds only test\n"
	<< "for linkage.\n"
	<< "If admixturescorefile is selected, then please unselect both\n"
	<< "allelicassociationscorefile and ancestryassociationscorefile\n";
    badOptions = true;
  }

  if( TestForMisspecifiedAlleleFreqs &&
      ( alleleFreqFilename.length()==0 && !(fixedallelefreqs) ) ){
    Log << "Cannot test for mis-specified allele frequencies unless allele frequencies are fixed.\n";
    badOptions = true;
  }

  if( TestForAffectedsOnly )
    if( RegType == Linear || RegType == Mlinear){
      Log << "ERROR: affectedsonly score test is not valid with a linear regression only."
	  << " This option will be ignored.\n";
      setTestForAffectedsOnly(false);
    }
    else   useroptions["likratiofile"] = "LikRatioFile.txt";
  if( TestForLinkageWithAncestry ){
    if(NumberOfOutcomes < 1){
      Log << "ERROR: ancestryassociation score test is not valid without a regression model."
	  << " This option will be ignored.\n";
      setTestForLinkageWithAncestry(false);
    }
  }
  if( TestForAllelicAssociation ){
    if( NumberOfOutcomes < 1 ){
      Log << "ERROR: allelic association score test is not valid without a regression model."
	  << " This option will be ignored.\n";
      setTestForAllelicAssociation(false);
    }
  }
  if( TestForHaplotypeAssociation && !TestForAllelicAssociation ){
    Log << "ERROR: Can't test for haplotype associations if allelicassociationscorefile is not specified"
	<< " This option will be ignored.\n";
    setTestForHaplotypeAssociation(false);
    }
  
  ScoreTestIndicator = (TestForAffectedsOnly || TestForLinkageWithAncestry || TestForAllelicAssociation || 
			TestForAdmixtureAssociation || TestForHaplotypeAssociation || TestForResidualAllelicAssoc);
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


  if(thermoIndicator) {
    // for thermo integration, NumAnnealedRuns is set to default value of 100 
    // if not specified as an option
    //if(NumAnnealedRuns==0) NumAnnealedRuns = 100;
    Log << "\nUsing thermodynamic integration to calculate marginal likelihood ";
    if(!TestOneIndivIndicator) Log << "for all individuals\n";
    else Log << "for first individual\n"; 
  }

  if(badOptions) return 1;
  else return 0;
}

//Note: requires Populations option to have already been set
void AdmixOptions::setInitAlpha(LogWriter &Log){
  _admixed.resize(2,(bool)(Populations>1));
  _symmetric = true;
  vector<double> alphatemp(Populations);
  Log.setDisplayMode(Quiet);

  //if no initalpha is specified, alpha for both gametes is initialised to 1.0 for each population  
  if( initalpha[0].size() == 0 && initalpha[1].size() == 0 ){
    fill( alphatemp.begin(), alphatemp.end(), 1.0);//fill alphatemp with 1s
    initalpha[0] = alphatemp; initalpha[1] = alphatemp;//put 2 copies of alphatemp in alpha
    if(!HapMixModelIndicator){
      if(!IndAdmixHierIndicator)
	Log << "Dirichlet parameters of prior on admixture: ";
      else 
	Log << "Initial value for population admixture (Dirichlet) parameter vector: ";
      for(int k = 0;k < Populations; ++k){Log << alphatemp[k] << " " ;}
      Log << "\n";
    }
  }
  //if only initalpha0 specified, sets initial values of alpha parameter vector for both gametes
  // if indadmixhiermodel=0, alpha values stay fixed
  else if( initalpha[0].size() > 0 && initalpha[1].size() == 0 ){
    _admixed[0] = CheckInitAlpha( initalpha[0] );
    initalpha[1] = initalpha[0];//put 2 copies of alpha[0] in alpha
    if(!IndAdmixHierIndicator)
      Log << "Dirichlet parameters of prior on admixture: ";
    else 
      Log << "Initial value for population admixture (Dirichlet) parameter vector: ";
    for(size_t k = 0;k < initalpha[0].size(); ++k){Log << initalpha[0][k] << " ";}
    Log << "\n";
  }
  //if both are specified and there is no hierarchical model on admixture,
  //paternal/gamete1 and maternal/gamete2 alphas are set to initalpha0 and initalpha1
  else if( !IndAdmixHierIndicator ){ 
    _admixed[0] = CheckInitAlpha( initalpha[0] );    //gamete 1
    _admixed[1] = CheckInitAlpha( initalpha[1] );    //gamete 2

    Log << "Dirichlet parameters of prior on maternal gamete admixture: ";
    for(size_t k = 0;k < initalpha[0].size(); ++k){Log << initalpha[0][k] << " ";}
    Log << "\n";
    
    Log << "Dirichlet parameters of prior on paternal gamete admixture: ";
    for(size_t k = 0;k < initalpha[1].size(); ++k){Log << initalpha[1][k] << " " ;}
    Log << "\n";
    
    _symmetric = false;
  }
  else{
    Log.setDisplayMode(On);
    Log << "ERROR: Can specify separate priors on admixture of each gamete only if indadmixhierindicator = 0\n";
    exit(1);
  }
}

bool AdmixOptions::CheckInitAlpha( const vector<double> &alphatemp)const
//returns indicator for admixture as indicated by initalpha   
// also check that Dirichlet parameter vector, if specified by user, has correct length
{
   bool admixed = true;
   int count = 0;
   for( size_t i = 0; i < alphatemp.size(); i++ )
      if( alphatemp[i] > 0.0 )
         count++;
   if( count == 1 )
     admixed = false;//unadmixed if eg 1 0 0
   if( alphatemp.size() != (unsigned)Populations ){
     cerr << "Error in specification of initalpha.\n";
     copy(alphatemp.begin(), alphatemp.end(), ostream_iterator<double>(cerr));//prints alphatemp
     cerr << endl;
     exit(0);
   }
   return admixed;
}

void AdmixOptions::ReadCommandLineArgs(const int argc, char** argv){
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
