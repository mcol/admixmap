/** 
 *   ADMIXMAP
 *   AdmixOptions.cc 
 *   Class to hold program options
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "AdmixOptions.h"
#include "LogWriter.h"
#include "StringSplitter.h"
#include "StringConvertor.h"
#include <getopt.h>    /* for getopt and getopt_long */
#include <string.h>
#include <sstream>
#include <numeric> // for checkInitAlpha

using namespace std;


/**
 *   Convert Cstrings to vector (used on initalpha option)
 */
static std::vector<double> CstrToVec2(const char* str)
{
  // Split string to separate elems.
  StringSplitter splitter;
  vector<string> elems = splitter.split(str);

  // Convert elements to doubles and fill resulting vector.
  vector<double> x(elems.size());
  for (size_t i = 0; i < elems.size(); i++) {
    x[i] = StringConvertor::toFloat(elems[i]);
  }    
  return x;
}

AdmixOptions::AdmixOptions()
{
  burnin = 1000;
  TotalSamples = 1000000;
  SampleEvery = 100;
  Seed = 1;
  AnalysisTypeIndicator = 0;
  TargetIndicator = 0;
  TruncPt = 99;
  Populations = 0;

  use_cout = 1; //should be bool
  TextIndicator = 1;//should be bool
  OutputFST = false;
  XOnlyAnalysis = false;
  isPedFile = false; 
  genotypesSexColumn = 0;
  locusForTestIndicator = false;
  LocusForTest = 0;
  fixedallelefreqs = false;
  RandomMatingModel = false;
  RhoIndicator = false;//corresponds to globalrho = 1;
  IndAdmixHierIndicator = true;//hierarchical model on ind admixture
  MLIndicator = false;//calculate marginal likelihood - valid only for analysistypeindicator < 0
  AnnealIndicator = false;
  ScoreTestIndicator = false; //indicator for any of the score tests in ScoreTests class
  TestForAdmixtureAssociation = false;
  StratificationTestIndicator = false;
  TestForAffectedsOnly = false;
  TestForAllelicAssociation = false;
  TestForSNPsInHaplotype = false;
  TestForDispersion = false;
  TestForLinkageWithAncestry = false;
  TestForMisspecifiedAlleleFreqs = false;
  TestForMisspecifiedAlleleFreqs2 = false;
  HWTest = false;
  OutputAlleleFreq = false;

  Rhoalpha = 3.0; //should be 3.0
  Rhobeta = 0.5;  // and 0.5
  alphamean = 1;  //  gamma(0.25, 0.25)
  alphavar = 16;   // 

  ResultsDir = "results";
  LogFilename = "log.txt";

  OptionValues["burnin"] = "1000";
  OptionValues["samples"] = "1000000";
  OptionValues["targetindicator"] = "0";
  OptionValues["coutindicator"] = "1";
  OptionValues["analysistypeindicator"] = "0";
  OptionValues["every"] = "100";
  OptionValues["fixedallelefreqs"] = "0";
  OptionValues["locusfortest"] = "0";
  OptionValues["logfile"] = "log.txt";
  OptionValues["resultsdir"] = "results";
  OptionValues["randommatingmodel"] = "0";
  OptionValues["globalrho"] = "1";
  OptionValues["indadmixhiermodel"] = "1";
  OptionValues["marglikelihood"] = "0";
  OptionValues["truncationpoint"] = "99";
  OptionValues["populations"] = "0";
  OptionValues["seed"] = "1";
  OptionValues["sumintensitiesalpha"] = "5.0";
  OptionValues["sumintensitiesbeta"] = "1.0";
  OptionValues["popadmixpriormean"] = "1.0";
  OptionValues["popadmixpriorvar"] = "1.0";
  OptionValues["xonlyanalysis"] = "0";
  OptionValues["textindicator"] = "1";
}

AdmixOptions::~AdmixOptions()
{
}

const char *AdmixOptions::getResultsDir() const{
  return ResultsDir.c_str();
}
const char *AdmixOptions::getAlleleFreqFilename() const
{
  return alleleFreqFilename.c_str();
}

long AdmixOptions::getBurnIn() const
{
  return burnin;
}

const char *AdmixOptions::getDICoutputFilename() const
{
  return DICoutputFilename.c_str();
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
  return fixedallelefreqs;
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

int AdmixOptions::useCOUT() const
{
  return use_cout;
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

bool AdmixOptions::getOutputAlleleFreq() const
{
  return OutputAlleleFreq;
}

int AdmixOptions::getAnalysisTypeIndicator() const
{
  return AnalysisTypeIndicator;
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

const char *AdmixOptions::getMLEFilename() const
{
  return MLEFilename.c_str();
}

bool AdmixOptions::isRandomMatingModel() const
{
  return RandomMatingModel;
}

bool AdmixOptions::getIndAdmixHierIndicator() const{
  return IndAdmixHierIndicator;
}
bool AdmixOptions::getMLIndicator()const{
  return MLIndicator;
}
bool AdmixOptions::getAnnealIndicator()const{
  return AnnealIndicator;
}

double AdmixOptions::getTruncPt() const
{
  return TruncPt;
}

bool AdmixOptions::getRhoIndicator() const
{
  return RhoIndicator;
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

long AdmixOptions::getSeed() const
{
  return Seed;
}
double AdmixOptions::getRho() const
{
  return Rhoalpha;
}
double AdmixOptions::getRhoalpha() const
{
  return Rhoalpha;
}
double AdmixOptions::getRhobeta() const
{
  return Rhobeta;
}
double AdmixOptions::getAlphamean() const{
  return alphamean;
}
double AdmixOptions::getAlphavar() const{
  return alphavar;
}
bool AdmixOptions::RhoFlatPrior() const{
  if( Rhoalpha==99 || ((Rhoalpha==1) && (Rhobeta==0)) ) return true;
  else return false;
}
bool AdmixOptions::logRhoFlatPrior() const{
  if( Rhoalpha==98 || ((Rhoalpha==0) && (Rhobeta==0)) ) return true;
  else return false;
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
    OptionValues.erase("affectedsonlyscorefile");
}

bool AdmixOptions::getTestForAllelicAssociation() const
{
  return TestForAllelicAssociation;
}

void AdmixOptions::setTestForAllelicAssociation(bool b){
  TestForAllelicAssociation = b;
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
    OptionValues.erase("ancestryassociationscorefile");
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs() const
{
  return TestForMisspecifiedAlleleFreqs;
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs2() const
{
  return TestForMisspecifiedAlleleFreqs2;
}

const char *AdmixOptions::getTestsForSNPsInHaplotypeOutputFilename() const
{
  return TestsForSNPsInHaplotypeOutputFilename.c_str();
}

bool AdmixOptions::getTestForSNPsInHaplotype() const
{
  return TestForSNPsInHaplotype;
}

void AdmixOptions::setTestForSNPsInHaplotype(bool b){
  TestForSNPsInHaplotype = b;
}

const char *AdmixOptions::getEtaPriorFilename() const
{
  return EtaPriorFilename.c_str();
}

int AdmixOptions::getTextIndicator() const
{
  return TextIndicator;
}

int AdmixOptions::sizeInitAlpha() const
{
  return alpha.size();
}

std::vector<double> AdmixOptions::getInitAlpha(int gamete) const
{
  return alpha[gamete];
}

bool AdmixOptions::IsPedFile() const
{
  return isPedFile;
}

void AdmixOptions::IsPedFile(bool i)
{
  isPedFile = i;
}

unsigned int AdmixOptions::getgenotypesSexColumn() const
{
  return genotypesSexColumn;
}

void AdmixOptions::setgenotypesSexColumn(unsigned int i)
{
  genotypesSexColumn = i;
}

bool AdmixOptions::getXOnlyAnalysis() const
{
  return XOnlyAnalysis;
}
bool AdmixOptions::isSymmetric()const{
  return _symmetric;
}
bool AdmixOptions::isAdmixed(unsigned gamete)const{
  return _admixed[gamete];
}

void AdmixOptions::SetOptions(int nargs,char** args)
{

  // This is the command-line parsing
  int c;

  /**
   * long_options is a pointer to the first element of an array of struct option
   * declared in <getopt.h> as
   * 
   *   struct option {
   *     const char *name;
   *     int has_arg;
   *     int *flag;
   *     int val;
   *   };
   * 
   * The meanings of the different fields are:
   * 
   *   name
   *     the name of the long option.
   * 
   *   has_arg
   *     no_argument (or 0) if the option does not take an argument,
   *     required_argument (or 1) if the option requires an argument,  or
   *     optional_argument  (or  2) if the option takes an optional argu-
   *     ment.
   * 
   *   flag
   *     specifies how results are returned for a long option.   If  flag
   *     is  NULL,  then  getopt_long()  returns  val.  (For example, the
   *     calling program may set val to the equivalent short option char-
   *     acter.)   Otherwise, getopt_long() returns 0, and flag points to
   *     a variable which is set to val if the option is found, but  left
   *     unchanged if the option is not found.
   * 
   *   val
   *     the value to return, or to load into the variable pointed to
   *     by flag.
   *
   * The last element of the array has to be filled with zeroes.
   */

  static struct option long_options[] = {
    // Required options
    {"analysistypeindicator",                 1, 0,  0 }, // int 0: 4
    {"locusfile",                             1, 0, 'l'}, // string
    {"genotypesfile",                         1, 0, 'g'}, // string
    {"burnin",                                1, 0, 'b'}, // long
    {"samples",                               1, 0, 's'}, // long
    {"every",                                 1, 0,  0 }, // long

    // Must specify one of the following
    {"populations",                           1, 0,  0 }, // int 1:populations
    {"priorallelefreqfile",                   1, 0,  0 }, // string
    {"allelefreqfile",                        1, 0, 'a'}, // string
    {"historicallelefreqfile",                1, 0,  0 }, // string

    // Required with certain analysistypeindicator values
    {"outcomevarfile",                        1, 0,  0 }, // string
    {"covariatesfile",                        1, 0,  0 }, // string
    {"mlefile",                               1, 0,  0 }, // string

    // Optional if specify outcomevarfile
    {"targetindicator",                       1, 0, 't'}, // int
    // 0: no. of cols in outcomvarfile

    //standard output files (optional)
    {"paramfile",                             1, 0, 'p'}, // string
    {"regparamfile",                          1, 0,  0 }, // string
    {"dispparamfile",                         1, 0,  0 }, // string
    {"logfile",                               1, 0,  0 }, // string
    {"indadmixturefile",                      1, 0, 'i'}, // string
    {"allelefreqoutputfile",                  1, 0,  0 }, // string
    {"ergodicaveragefile",                    1, 0, 'e'}, // string

    //optional results directory name option - default is 'results'
    {"resultsdir",                            1, 0,  0 }, //string
   
    // Extra output options
    {"allelicassociationscorefile",           1, 0,  0 }, // string
    {"ancestryassociationscorefile",          1, 0,  0 }, // string
    {"affectedsonlyscorefile",                1, 0,  0 }, // string
    {"stratificationtestfile",                1, 0, 'd'}, // string
    {"admixturescorefile",                    1, 0,  0 }, // string
    {"haplotypeassociationscorefile",         1, 0,  0 }, // string
    {"locusfortest",                          1, 0,  0 }, // int 0: no. of composite loci - 1
    {"allelefreqscorefile",                   1, 0,  0 }, // string
    {"allelefreqscorefile2",                  1, 0,  0 }, // string
    {"dispersiontestfile",                    1, 0,  0 }, // string
    {"fstoutputfile",                         1, 0,  0 }, // string
    {"hwscoretestfile",                       1, 0,  0 }, // string

    // Other options
    {"coutindicator",                         1, 0, 'c'}, // int 0: 1
    {"help",                                  0, 0, 'h'}, // NONE
    {"randommatingmodel",                     1, 0,  0 }, // int 0: 1
    {"globalrho",                             1, 0,  0 }, // int 0: 1
    {"indadmixhiermodel",                     1, 0,  0 }, // int 0: 1
    {"marglikelihood",                        1, 0,  0 }, // int 0: 1
    {"anneal",                                1, 0,  0 }, // int 0: 1
    {"reportedancestry",                      1, 0, 'r'}, // string 
    {"seed",                                  1, 0,  0 }, // long
    {"etapriorfile",                          1, 0,  0 }, // string      
    {"textindicator",                         1, 0,  0 }, // int
    {"sumintensitiesalpha",                   1, 0,  0 }, // double
    {"sumintensitiesbeta",                    1, 0,  0 }, // double
    {"popadmixpriormean",                     1, 0,  0 }, //double
    {"popadmixpriorvar",                      1, 0,  0 }, //double
    {"truncationpoint",                       1, 0,  0 }, // double
    {"initalpha0",                            1, 0,  0 }, // double
    {"initalpha1",                            1, 0,  0 }, // double
    {"fixedallelefreqs",                      1, 0,  0 }, // long
    {"xonlyanalysis",                         1, 0,  0 }, // long
    {0, 0, 0, 0}    // marks end of array
  };

  //options specified on command line
  while (1) {
    int option_index = 0;

    c = getopt_long (nargs, args, "a:b:c:d:e:g:h:i:l:o:p:r:s:t:",
		     long_options, &option_index);

    string long_option_name = long_options[option_index].name;

    if (c == -1)
      break;

    switch (c) {
    case 'a': // allelefreqfile
      { alleleFreqFilename = optarg;OptionValues["allelefreqfile"]=optarg;}
      break;

    case 'b': // burnin
      { burnin = strtol(optarg, NULL, 10);OptionValues["burnin"]=optarg;}
      break;

    case 'c': // coutindicator
      { use_cout = (int)strtol(optarg, NULL, 10);OptionValues["coutindicator"]=optarg;}
      break;

    case 'd': // stratificationtestfile
      DICoutputFilename = optarg;OptionValues["stratificationtestfile"]=optarg;
      StratificationTestIndicator = true;//OptionValues["StratificationTestIndicator"]="1";
      break;

    case 'e': // ergodicaveragefile
      { ErgodicAverageFilename = optarg;OptionValues["ergodicaveragefile"]=optarg;}
      break;

    case 'g': // genotypesfile
      { GenotypesFilename = optarg;OptionValues["genotypesfile"]=optarg;}
      break;

    case 'h': // help
      cout << "Usage: " << args[0] << " [options]" << endl;
      cout << " some help text goes here" << endl;
//outputlist of user options to console
      exit(0);

    case 'i': // indadmixturefile
      { IndAdmixtureFilename = optarg;OptionValues["indadmixturefile"]=optarg;}
      break;

    case 'l': // locusfile
      { LocusFilename = optarg;OptionValues["locusfile"]=optarg;}
      break;

    case 'p': // paramfile
      { ParameterFilename = optarg;OptionValues["paramfile"]=optarg;}
      break;

    case 'r': // reportedancestry
      { ReportedAncestryFilename = optarg;OptionValues["reportedancestry"]=optarg;}
      break;

    case 's': // samples
      { TotalSamples = strtol(optarg, NULL, 10);OptionValues["samples"]=optarg;}
      break;

    case 't': // targetindicator
      { TargetIndicator = (int)strtol(optarg, NULL, 10);OptionValues["targetindicator"]=optarg;}
      break;

    case '?':
      exit(1);

    case 0:
      if(long_option_name == "resultsdir"){
	 ResultsDir = optarg;OptionValues["resultsdir"]=optarg;
      }else if (long_option_name == "regparamfile") {
	 RegressionOutputFilename = optarg;OptionValues["regparamfile"]=optarg;
      }else if (long_option_name == "dispparamfile") {
	 EtaOutputFilename = optarg; OptionValues["dispparamfile"]=optarg;
      }else  if (long_option_name == "admixturescorefile") {
	 AssocScoreFilename = optarg;OptionValues["admixturescorefile"]=optarg;
	 TestForAdmixtureAssociation = true; ScoreTestIndicator = true;
      } else if (long_option_name == "affectedsonlyscorefile") {
	 AffectedsOnlyScoreFilename = optarg;OptionValues["affectedsonlyscorefile"]=optarg;
	 TestForAffectedsOnly = true; ScoreTestIndicator = true;
      } else if (long_option_name == "allelicassociationscorefile") {
	 AllelicAssociationScoreFilename = optarg;OptionValues["allelicassociationscorefile"]=optarg;
	 TestForAllelicAssociation = true; ScoreTestIndicator = true;
      } else if (long_option_name == "allelefreqoutputfile") {
	 AlleleFreqOutputFilename = optarg;OptionValues["allelefreqoutputfile"]=optarg;
	 OutputAlleleFreq = true;
      } else if (long_option_name == "allelefreqscorefile") {
	 AlleleFreqScoreFilename = optarg;OptionValues["allelefreqscorefile"]=optarg;
	 TestForMisspecifiedAlleleFreqs = true;
      } else if (long_option_name == "allelefreqscorefile2"){
	 AlleleFreqScoreFilename2 = optarg;OptionValues["allelefreqscorefile2"]=optarg;
	 TestForMisspecifiedAlleleFreqs2 = true;
      } else if (long_option_name == "analysistypeindicator") {
	 AnalysisTypeIndicator = (int)strtol(optarg, NULL, 10);OptionValues["analysistypeindicator"]=optarg;
      } else if (long_option_name == "covariatesfile") {
	 CovariatesFilename = optarg;OptionValues["covariatesfile"]=optarg;
      } else if (long_option_name == "mlefile") {
	 MLEFilename = optarg;OptionValues["mlefile"]=optarg;
      } else if (long_option_name == "dispersiontestfile") {
	 DispersionTestFilename = optarg;OptionValues["dispersiontestfile"]=optarg;
	 TestForDispersion = true;
      } else if (long_option_name == "every") {
	 SampleEvery = strtol(optarg, NULL, 10);OptionValues["every"]=optarg;
      } else if (long_option_name == "fstoutputfile") {
	 FSTOutputFilename = optarg;OptionValues["fstoutputfile"]=optarg;
	 OutputFST = true;
      } else if (long_option_name == "hwscoretestfile") {
	 HWTestFilename = optarg;OptionValues["hwscoretestfile"]=optarg;
	 HWTest = true;
      } else if (long_option_name == "fixedallelefreqs") {
	if (strtol(optarg, NULL, 10) == 1) {
	  fixedallelefreqs = true;OptionValues["fixedallelefreqs"]="1";
	}
      } else if (long_option_name == "xonlyanalysis") {
	if (strtol(optarg, NULL, 10) == 1) {
	  XOnlyAnalysis = true;OptionValues["xonlyanalysis"]="1";
	}
      } else if (long_option_name == "historicallelefreqfile") {
	 HistoricalAlleleFreqFilename = optarg;OptionValues["historicallelefreqfile"]=optarg;
      } else if (long_option_name == "locusfortest") {
	 LocusForTest = (int)strtol(optarg, NULL, 10);OptionValues["locusfortest"]=optarg;
	 locusForTestIndicator = true;
      } else if (long_option_name == "ancestryassociationscorefile") {
	 AncestryAssociationScoreFilename = optarg;OptionValues["ancestryassociationscorefile"]=optarg;
	 TestForLinkageWithAncestry = true; ScoreTestIndicator = true;
      } else if (long_option_name == "logfile") {
	 LogFilename = optarg;OptionValues["logfile"]=optarg;
      } else if (long_option_name == "randommatingmodel") {
	if (strtol(optarg, NULL, 10) == 1) {
	  RandomMatingModel = true;OptionValues["randommatingmodel"]="1";
	}
      } else if (long_option_name == "indadmixhiermodel") {
	if (strtol(optarg, NULL, 10) == 0) {
	  IndAdmixHierIndicator = false;OptionValues["indadmixhiermodel"]="0";
	}
      }else if (long_option_name == "marglikelihood") {
	if (strtol(optarg, NULL, 10) == 1) {
	  MLIndicator = true;OptionValues["marglikelihood"]="1";
	}
      }else if (long_option_name == "anneal") {
	if (strtol(optarg, NULL, 10) == 1) {
	  AnnealIndicator = true;OptionValues["anneal"]="1";
	}
      }else if (long_option_name == "globalrho") {
	if (strtol(optarg, NULL, 10) == 1) {
	  RhoIndicator = false;OptionValues["globalrho"]="1";
	} else if (strtol(optarg, NULL, 10) == 0) {
	  RhoIndicator = true;OptionValues["globalrho"]="0";
	} else {
	  cerr << "Set global rho to 0 or 1.\n";
	  exit(1);
	}
      } else if (long_option_name == "outcomevarfile") {
	 OutcomeVarFilename = optarg;OptionValues["outcomevarfile"]=optarg;
      } else if (long_option_name == "populations") {
	setPopulations((int)strtol(optarg, NULL, 10));OptionValues["populations"]=optarg;
      } else if (long_option_name == "priorallelefreqfile") {
	 PriorAlleleFreqFilename = optarg;OptionValues["priorallelefreqfile"]=optarg;
      } else if (long_option_name == "seed") {
	 Seed = strtol(optarg, NULL, 10);OptionValues["seed"]=optarg;
      } else if (long_option_name == "sumintensitiesalpha") {
	 Rhoalpha = strtod(optarg, NULL);OptionValues["sumintensitiesalpha"]=optarg;
      } else if (long_option_name == "sumintensitiesbeta") {
	 Rhobeta = strtod(optarg, NULL);OptionValues["sumintensitiesbeta"]=optarg;
      } else if (long_option_name == "popadmixpriormean") {
	 alphamean = strtod(optarg, NULL);OptionValues["popadmixpriormean"]=optarg;
      } else if (long_option_name == "popadmixpriorvar") {
	 alphavar = strtod(optarg, NULL);OptionValues["popadmixpriorvar"]=optarg;
      } else if (long_option_name == "truncationpoint") {
	 TruncPt = strtod(optarg, NULL);OptionValues["truncationpoint"]=optarg;
      } else if (long_option_name == "haplotypeassociationscorefile") {
	 TestsForSNPsInHaplotypeOutputFilename = optarg;OptionValues["haplotypeassociationscorefile"]=optarg;
	 TestForSNPsInHaplotype = true; ScoreTestIndicator = true;
      } else if (long_option_name == "etapriorfile") {
	 EtaPriorFilename = optarg;OptionValues["etapriorfile"]=optarg;
      } else if (long_option_name == "initalpha0" ) {
	 alpha.push_back(CstrToVec2(optarg));OptionValues["initalpha0"]=optarg;
      } else if (long_option_name == "initalpha1") {
	 alpha.push_back(CstrToVec2(optarg));OptionValues["initalpha1"]=optarg;
      } else if (long_option_name == "textindicator") {
	 TextIndicator = (int)strtol(optarg, NULL, 10);OptionValues["textindicator"]=optarg;
      } else {
	cerr << "Unknown option: " << long_option_name;
	if (optarg) {
	  cerr << " with arg: " << optarg;
	}
	cerr << endl;
	exit(1);
      }
      break;

    default:
      fprintf(stderr, "?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (optind < nargs) {
    printf ("non-option ARGV-elements: ");
        
    while (optind < nargs) {
      printf ("%s ", args[optind++]);
    }

    printf ("\n");
  }
  // command-line parsing ends here
  SetOutputNames();
}

void AdmixOptions::SetOutputNames(){

  //prefix output files with ResultsDir
  if ( LogFilename != "") LogFilename = ResultsDir + "/" + LogFilename;
  if ( ParameterFilename != "") ParameterFilename = ResultsDir + "/" + ParameterFilename;
  if ( RegressionOutputFilename != "") RegressionOutputFilename = ResultsDir + "/" + RegressionOutputFilename;
  if ( EtaOutputFilename != "") EtaOutputFilename = ResultsDir + "/" + EtaOutputFilename;
  if ( AllelicAssociationScoreFilename != "") AllelicAssociationScoreFilename = ResultsDir + "/" + AllelicAssociationScoreFilename;
  if ( AncestryAssociationScoreFilename != "") AncestryAssociationScoreFilename = ResultsDir + "/" + AncestryAssociationScoreFilename;
  if ( AffectedsOnlyScoreFilename != "") AffectedsOnlyScoreFilename = ResultsDir + "/" + AffectedsOnlyScoreFilename;
  if ( IndAdmixtureFilename != "") IndAdmixtureFilename = ResultsDir + "/" + IndAdmixtureFilename;
  if ( ErgodicAverageFilename != "") ErgodicAverageFilename = ResultsDir + "/" + ErgodicAverageFilename;
  if ( DICoutputFilename != "") DICoutputFilename = ResultsDir + "/" + DICoutputFilename;
  if ( AssocScoreFilename != "") AssocScoreFilename = ResultsDir + "/" + AssocScoreFilename;
  if ( TestsForSNPsInHaplotypeOutputFilename != "") TestsForSNPsInHaplotypeOutputFilename = ResultsDir + "/" + TestsForSNPsInHaplotypeOutputFilename;
  if ( AlleleFreqOutputFilename != "") AlleleFreqOutputFilename = ResultsDir + "/" + AlleleFreqOutputFilename;
  if ( AlleleFreqScoreFilename2 != "") AlleleFreqScoreFilename2 = ResultsDir + "/" + AlleleFreqScoreFilename2;
  if ( AlleleFreqScoreFilename != "") AlleleFreqScoreFilename = ResultsDir + "/" + AlleleFreqScoreFilename;
  if ( DispersionTestFilename != "") DispersionTestFilename = ResultsDir + "/" + DispersionTestFilename;
  if ( FSTOutputFilename != "") FSTOutputFilename = ResultsDir + "/" + FSTOutputFilename;
  if ( HWTestFilename != "") HWTestFilename = ResultsDir + "/" + HWTestFilename;
  
}

void AdmixOptions::PrintOptions(){
  //set populations value in case it has changed
  //NB do similar for any option that can be changed outside AdmixOptions
  std::ostringstream s;
  if (s << Populations) // conversion worked
    {
    OptionValues["populations"] = (char *)s.str().c_str();
    }

  //Now output Options table to args.txt
  string ss;
  ss = ResultsDir + "/args.txt";
  ofstream argstream(ss.c_str());

  for( OptionMap::iterator p= OptionValues.begin(); p!=OptionValues.end(); p++) {
    argstream << (*p).first << "=" << (*p).second <<endl;
  }
  argstream.close();
}

int AdmixOptions::checkOptions(LogWriter *Log){

  // **** analysis type  ****
  if (AnalysisTypeIndicator == 0)
    {
      Log->logmsg(true,"Affecteds only analysis.\n");
      if( !TestForAffectedsOnly ){
        Log->logmsg(true,"Must specify affectedsonlyscorefile.\n");
      }
    }
  else if (AnalysisTypeIndicator == 1)
    {
      Log->logmsg(true,"Cross sectional analysis, no outcome.\n");
    }
  else if (AnalysisTypeIndicator == 2)
    {
      Log->logmsg(true,"Cross sectional analysis, continuous outcome.\n");
      if( OutcomeVarFilename.length() == 0 )
	{
	   Log->logmsg(true,"Must specify outcomevar file.\n");
	   exit(0);
	}
    }
  else if (AnalysisTypeIndicator == 3)
    {
      Log->logmsg(true,"Cross sectional analysis, binary outcome.\n");
      if( OutcomeVarFilename.length() == 0 )
	{
	  Log->logmsg(true,"Must specify outcomevar file.\n");
	  exit(0);
	}
    }
  else if (AnalysisTypeIndicator == 4)
    {
      Log->logmsg(true,"Case control analysis.\n");
      if( OutcomeVarFilename.length() == 0  )
	{
	   Log->logmsg(true,"Must specify outcomevar file.\n");
	   exit(0);
	}
    }
  else if (AnalysisTypeIndicator == 5)
    {
      Log->logmsg(true,"Cross sectional analysis, multiple outcome.\n");
      if( OutcomeVarFilename.length() == 0  )
	{
	  Log->logmsg(true,"Must specify outcomevar file.\n");
	  exit(0);
	}
    }
  else if (AnalysisTypeIndicator == -1 || AnalysisTypeIndicator == -2)
    {
      Log->logmsg(true,"One individual analysis");
      if(MLIndicator)Log->logmsg(true, " with marginal likelihood calculation");
      Log->logmsg(true, "\n");
    }

  else
    {
      Log->logmsg(true,"Unknown analysis type: ");
      Log->logmsg(true, AnalysisTypeIndicator);
      Log->logmsg(true, "\n");
      exit(0);
    }
  if(AnalysisTypeIndicator < 2 && RegressionOutputFilename.length() > 0){
    Log->logmsg(true, "ERROR: regparamfile option is not valid without a regression model\n");
    Log->logmsg(true, "\tThis option will be ignored\n");
    RegressionOutputFilename = "";
    OptionValues.erase("regparamfile");
  }


  // **** Hierarchical model on ind admixture ****
  if (!IndAdmixHierIndicator)
    {
      Log->logmsg(true,"No hierarchical model for individuals.\n");

      if(ParameterFilename.length() > 0 ){
	Log->logmsg(true, "ERROR: paramfile option is not valid with indadmixhierindicator = 0\n");
	Log->logmsg(true, "\tThis option will be ignored\n");
	 ParameterFilename = "";
	OptionValues.erase("paramfile");
      }
      if(RegressionOutputFilename.length() > 0){
	Log->logmsg(true, "ERROR: regparamfile option is not valid with indadmixhierindicator = 0\n");
	Log->logmsg(true, "\tThis option will be ignored\n");
	RegressionOutputFilename = "";
	OptionValues.erase("regparamfile");
	 }
      if(EtaOutputFilename.length() > 0 ){
	Log->logmsg(true, "ERROR: dispparamfile option is not valid with indadmixhierindicator = 0\n");
	Log->logmsg(true, "\tThis option will be ignored\n");
	EtaOutputFilename = "";
	OptionValues.erase("dispparamfile");
      }
    }

  // **** Random Mating Model **** 
  if(RandomMatingModel )
    Log->logmsg(true,"Model assuming random mating.\n");
  else 
    Log->logmsg(true,"Model assuming assortative mating.\n");

  // **** global rho ****
  if( !RhoIndicator )
    Log->logmsg(true,"Model with global sumintensities.\n");
  else if( RandomMatingModel )
    Log->logmsg(true,"Model with gamete specific sumintensities.\n");
  else
    Log->logmsg(true,"Model with individual specific sumintensities.\n");

  // **** Marginal Likelihood ****
  if(MLIndicator){
    if(AnalysisTypeIndicator >= 0){
      Log->logmsg(true, "Error: Cannot calculate marginal likelihood with analysis type ");
      Log->logmsg(true, AnalysisTypeIndicator);
      Log->logmsg(true, "\n");
      exit(0);
    }
    //change this when marginal likelihood can be calculated for other type of model
    //else Log->logmsg(true,"Analysis with marginal likelihood calculation\n");
  }


  // **** Check whether genotypes file has been specified ****
  if ( GenotypesFilename.length() == 0 )
    {
      Log->logmsg(true,"Must specify genotypesfile.\n");
      exit( 1 );
    }
  // **** Check whether locus file has been specified ****
  if ( LocusFilename.length() == 0 )
    {
      Log->logmsg(true,"Must specify locusfile.\n");
      exit( 1 );
    }

  // **** model for allele freqs ****
  if(Populations > 0 )
    {
      Log->logmsg(true,"No allelefreq filename or priorallelefreq filename given.\n");
      Log->logmsg(true,"Default priors will be set for the allele frequencies with ");
      Log->logmsg(true, Populations);
      Log->logmsg(true," populations.\n");
    }
  else if( alleleFreqFilename.length() ||
           (PriorAlleleFreqFilename.length() && fixedallelefreqs ) ){
    Log->logmsg(true,"Analysis with fixed allele frequencies.\n");
    if(OutputAlleleFreq){
      Log->logmsg(true, "ERROR: allelefreqoutputfile option is invalid with fixed allele frequencies\n");
      Log->logmsg(true, "       this option will be ignored\n");
      OptionValues.erase("allelefreqoutputfile");
      OutputAlleleFreq = false;
    }
  }
  else if( PriorAlleleFreqFilename.length() && !fixedallelefreqs ){
    Log->logmsg(true,"Analysis with prior allele frequencies.\n");
  }
  else if( HistoricalAlleleFreqFilename.length() > 0 ){
    Log->logmsg(true,"Analysis with dispersion model for allele frequencies.\n");
  }
  if( (FSTOutputFilename.length() > 0) && (HistoricalAlleleFreqFilename.length() == 0) ){
    Log->logmsg(true, "ERROR: fstoutputfile option is only valid with historicallelefreqfile option\n");
    Log->logmsg(true, "       this option will be ignored\n");
    OutputFST = false;
    FSTOutputFilename = "";
    OptionValues.erase("fstoutputfile");
  }



  // **** score tests ****
  if( TestForLinkageWithAncestry && Populations == 1 ){
    Log->logmsg(true,"Cannot test for linkage with ancestry with 1 population.\n");
    exit(0);
  }
  if(TestForAdmixtureAssociation &&
      ( TestForLinkageWithAncestry || TestForAllelicAssociation ) ){
    Log->logmsg(true,"Cannot test for linkage with ancestry or allelic association\n");
    Log->logmsg(true,"with score test for association. Can only use affecteds only test\n");
    Log->logmsg(true,"for linkage.\n");
    Log->logmsg(true,"If admixturescorefile is selected, then please unselect both\n");
    Log->logmsg(true,"allelicassociationscorefile and ancestryassociationscorefile\n");
    exit(1);
  }

  if( TestForMisspecifiedAlleleFreqs &&
      ( alleleFreqFilename.length()==0 && !(fixedallelefreqs) ) ){
    Log->logmsg(true,"Cannot test for mis-specified allele frequencies unless allele frequencies are fixed.\n");
    exit(1);
  }

  if( TestForAffectedsOnly )
    if(!(AnalysisTypeIndicator==0 || AnalysisTypeIndicator==3 || AnalysisTypeIndicator==4)){
      Log->logmsg(true,"ERROR: affectedsonly score test is only valid with analysistypeindicator 0, 3, or 4.");
      Log->logmsg(true," This option will be ignored.\n");
      setTestForAffectedsOnly(false);
    }
  if( TestForLinkageWithAncestry ){
    if(AnalysisTypeIndicator < 2){
      Log->logmsg(true,"ERROR: ancestryassociation score test is not valid with analysistypeindicator < 2 .");
      Log->logmsg(true," This option will be ignored.\n");
      setTestForLinkageWithAncestry(false);
    }
  }
  if( TestForAllelicAssociation ){
    if( AnalysisTypeIndicator < 2 ){
      Log->logmsg(true,"ERROR: allelic association score test is not valid with analysistypeindicator < 2.");
      Log->logmsg(true," This option will be ignored.\n");
      setTestForAllelicAssociation(false);
    }
  }
      
  return 1;
}
//Note: requires Populations option to have already been set
vector<vector<double> > AdmixOptions::getAndCheckInitAlpha(LogWriter *Log){
  _admixed.resize(2,true);
  _symmetric = true;
  vector<double> alphatemp(Populations);
  vector<vector<double> > initalpha;

  //if no initalpha is specified, alpha for both gametes is initialised to 1.0 for each population  
  if( alpha.size() == 0 ){
    fill( alphatemp.begin(), alphatemp.end(), 1.0);//fill alphatemp with 1s
    initalpha.push_back(alphatemp);initalpha.push_back(alphatemp);//put 2 copies of alphatemp in alpha
    Log->logmsg(false,  "Initial value for population admixture (Dirichlet) parameter vector: ");
    for(int k = 0;k < Populations; ++k){Log->logmsg(false,alphatemp[k]);Log->logmsg(false," ");}
    Log->logmsg(false,"\n");
  }
  //if exactly one of initalpha0 or initalpha1 is specified, sets initial values of alpha parameter vector for both gametes
  // if indadmixhiermodel=0, alpha values stay fixed
  else if( alpha.size() == 1 ){
    _admixed[0] = CheckInitAlpha( alpha[0] );
    initalpha.push_back(alpha[0]);initalpha.push_back(alpha[0]);//put 2 copies of alpha[0] in alpha
    Log->logmsg(false, "Initial value for population admixture (Dirichlet) parameter vector: ");
    for(size_t k = 0;k < alpha[0].size(); ++k){Log->logmsg(false,alpha[0][k]);Log->logmsg(false," ");}
    Log->logmsg(false,"\n");
  }
  //if both are specified and analysis is for a single individual,
  //paternal/gamete1 and maternal/gamete2 alphas are set to initalpha0 and initalpha1
  else if( AnalysisTypeIndicator < 0 ){ // should be if indadmixhiermodel=0
    _admixed[0] = CheckInitAlpha( alpha[0] );    //gamete 1
    _admixed[1] = CheckInitAlpha( alpha[1] );    //gamete 2

    initalpha = alpha;
    Log->logmsg(false, "Dirichlet prior for paternal gamete admixture: ");
    for(size_t k = 0;k < alpha[0].size(); ++k){Log->logmsg(false,alpha[0][k]);Log->logmsg(false," ");}
    Log->logmsg(false,"\n");
    
    Log->logmsg(false, "Dirichlet prior for maternal gamete admixture: ");
    for(size_t k = 0;k < alpha[1].size(); ++k){Log->logmsg(false,alpha[1][k]);Log->logmsg(false," ");}
    Log->logmsg(false,"\n");
    
    _symmetric = false;
  }
  else{
    Log->logmsg(true,"ERROR: Can specify separate priors on admixture of each gamete only if analysing single individual\n");
    exit(1);
  }
  return initalpha;
}

bool AdmixOptions::CheckInitAlpha( vector<double> &alphatemp)
// check that Dirichlet parameter vector, if specified by user, has correct length
//returns indicator for admixture as indicated by initalpha   
{
   bool admixed = true;
   int count = 0;
   for( size_t i = 0; i < alphatemp.size(); i++ )
      if( alphatemp[i] != 0.0 )
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
