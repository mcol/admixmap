#include "AdmixOptions.h"
#include "LogWriter.h"
#include "StringSplitter.h"
#include "StringConvertor.h"
#include <getopt.h>    /* for getopt and getopt_long */
#include <string.h>
#include <sstream>

using namespace std;


/**
 *   Convert Cstrings to vector (used on initalpha option)
 */
static Vector_d CstrToVec2(const char* str)
{
  // Split string to separate elems.
  StringSplitter splitter;
  vector<string> elems = splitter.split(str);

  // Convert elements to doubles and fill resulting vector.
  Vector_d x;
  x.SetNumberOfElements(elems.size());
  for (size_t i = 0; i < elems.size(); i++) {
    x(i) = StringConvertor::toFloat(elems[i]);
  }    
  return x;
}

struct AdmixOptions::Options
{
  string alleleFreqFilename;
  long burnin;
  string DICoutputFilename;
  string ErgodicAverageFilename;
  string GeneInfoFilename;
  string GeneticDataFilename;
  string IndAdmixtureFilename;
  string ParameterFilename;
  string RegressionOutputFilename;
  string EtaOutputFilename;
  string ReportedAncestryFilename;
  long TotalSamples;
  int TargetIndicator;
  int use_cout;

  string ResultsDir;
  string AffectedsOnlyScoreFilename;
  string AlleleFreqOutputFilename;
  string AlleleFreqScoreFilename;
  string AlleleFreqScoreFilename2;
  int AnalysisTypeIndicator;
  string AssocScoreFilename;
  long SampleEvery;
  string DispersionTestFilename;
  string HistoricalAlleleFreqFilename;
  string InputFilename;
  string MLEFilename;
  bool locusForTestIndicator;
  bool fixedallelefreqs;
  int LocusForTest;
  string LocusScoreFilename;
  string AncestryAssociationScoreFilename;
  string LogFilename;
  bool modelIndicator;//random mating model
  bool RhoIndicator;// global rho
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  bool MLIndicator;//calculate marginal likelihood - valid only for analysistypeindicator < 0
  double TruncPt;
  int Populations;
  string PriorAlleleFreqFilename;
  bool ScoreTestIndicator;
  long Seed;
  double Rho;
  std::vector<Vector_d> alpha;
  bool StratificationTestIndicator;
  string TargetFilename;
  bool TestForAffectedsOnly;
  bool TestForAllelicAssociation;
  bool TestForSNPsInHaplotype;
  bool TestForDispersion;
  bool TestForLinkageWithAncestry;
  bool TestForMisspecifiedAlleleFreqs;
  bool TestForMisspecifiedAlleleFreqs2;
  bool OutputAlleleFreq;
  string TestsForSNPsInHaplotypeOutputFilename;
  string EtaPriorFilename;
  string TestTargetFilename;
  int TextIndicator;
  string FSTOutputFilename;
  bool OutputFST;
  bool runRScript;
  bool XOnlyAnalysis;
  string argsfile;
  unsigned int isPedFile;
  unsigned int genotypesSexColumn;
};

AdmixOptions::AdmixOptions()
{
  imp = new Options;
  imp->burnin = 1000;
  imp->TotalSamples = 1000000;
  imp->TargetIndicator = 0;
  imp->use_cout = 1;

  imp->AnalysisTypeIndicator = 0;
  imp->SampleEvery = 100;
  imp->locusForTestIndicator = false;
  imp->fixedallelefreqs = false;
  imp->LocusForTest = 0;
  imp->LogFilename = "log.txt";
  imp->ResultsDir = "results";

  imp->modelIndicator = false;
  imp->RhoIndicator = false;//corresponds to globalrho = 1;
  imp->IndAdmixHierIndicator = true;
  imp->MLIndicator = false;
  imp->TruncPt = 99;
  imp->Populations = 0;

  imp->ScoreTestIndicator = false;
  imp->Seed = 1;
  imp->Rho = 5.0;
  imp->StratificationTestIndicator = false;

  imp->TestForAffectedsOnly = false;
  imp->TestForAllelicAssociation = false;
  imp->TestForSNPsInHaplotype = false;
  imp->TestForDispersion = false;
  imp->TestForLinkageWithAncestry = false;
  imp->TestForMisspecifiedAlleleFreqs = false;
  imp->TestForMisspecifiedAlleleFreqs2 = false;
  imp->OutputAlleleFreq = false;
  imp->OutputFST=false;    
  imp->XOnlyAnalysis = false;
  imp->TextIndicator = 1;
  imp->isPedFile = 0;
  imp->genotypesSexColumn = 0;
  OptionValues["burnin"] = "1000";
  OptionValues["samples"] = "1000000";
  OptionValues["targetindicator"] = "0";
  OptionValues["coutindicator"] = "1";
  OptionValues["analysistypeindicator"] = "0";
  OptionValues["every"] = "100";
  OptionValues["locusForTestIndicator"] = "0";
  OptionValues["fixedallelefreqs"] = "0";
  OptionValues["locusfortest"] = "0";
  OptionValues["logfile"] = "log.txt";
  OptionValues["resultsdir"] = "results";
  OptionValues["randommatingmodel"] = "0";
  OptionValues["globalrho"] = "1";
  OptionValues["IndAdmixHierModel"] = "1";
  OptionValues["MargLikelihood"] = "0";
  OptionValues["truncationpoint"] = "99";
  OptionValues["populations"] = "0";
  OptionValues["ScoreTestIndicator"] = "0";
  OptionValues["seed"] = "1";
  OptionValues["sumintensities"] = "5.0";
  OptionValues["StratificationTestIndicator"] = "0";
  OptionValues["TestForAffectedsOnly"] = "0";
  OptionValues["TestForAllelicAssociation"] = "0";
  OptionValues["TestForSNPsInHaplotype"] = "0";
  OptionValues["TestForDispersion"] = "0";
  OptionValues["TestForLinkageWithAncestry"] = "0";
  OptionValues["TestForMisspecifiedAlleleFreqs"] = "0";
  OptionValues["TestForMisspecifiedAlleleFreqs2"] = "0";
  OptionValues["OutputFST"]="0";
  OptionValues["xonlyanalysis"] = "0";
  OptionValues["textindicator"] = "1";
  OptionValues["isPedFile"] = "0";
  OptionValues["genotypesSexColumn"]= "0";
}

AdmixOptions::~AdmixOptions()
{
  delete imp;
}

const char *AdmixOptions::getResultsDir() const{
  return imp->ResultsDir.c_str();
}
const char *AdmixOptions::getAlleleFreqFilename() const
{
  return imp->alleleFreqFilename.c_str();
}

long AdmixOptions::getBurnIn() const
{
  return imp->burnin;
}

const char *AdmixOptions::getDICoutputFilename() const
{
  return imp->DICoutputFilename.c_str();
}

const char *AdmixOptions::getErgodicAverageFilename() const
{
  return imp->ErgodicAverageFilename.c_str();
}

const char *AdmixOptions::getFSTOutputFilename() const
{
  return imp->FSTOutputFilename.c_str();
}

bool AdmixOptions::getOutputFST() const
{
  return imp->OutputFST;
}

bool AdmixOptions::getFixedAlleleFreqs() const
{
  return imp->fixedallelefreqs;
}

const char *AdmixOptions::getGeneInfoFilename() const
{
  return imp->GeneInfoFilename.c_str();
}

const char *AdmixOptions::getGeneticDataFilename() const
{
  return imp->GeneticDataFilename.c_str();
}

const char *AdmixOptions::getIndAdmixtureFilename() const
{
  return imp->IndAdmixtureFilename.c_str();
}

const char *AdmixOptions::getParameterFilename() const
{
  return imp->ParameterFilename.c_str();
}

const char *AdmixOptions::getRegressionOutputFilename() const
{
  return imp->RegressionOutputFilename.c_str();
}

const char *AdmixOptions::getEtaOutputFilename() const
{
  return imp->EtaOutputFilename.c_str();
}

const char *AdmixOptions::getReportedAncestryFilename() const
{
  return imp->ReportedAncestryFilename.c_str();
}

long AdmixOptions::getTotalSamples() const
{
  return imp->TotalSamples;
}

int AdmixOptions::getTargetIndicator() const
{
  return imp->TargetIndicator;
}

int AdmixOptions::useCOUT() const
{
  return imp->use_cout;
}

const char *AdmixOptions::getAffectedsOnlyScoreFilename() const
{
  return imp->AffectedsOnlyScoreFilename.c_str();
}

const char *AdmixOptions::getAlleleFreqScoreFilename() const
{
  return imp->AlleleFreqScoreFilename.c_str();
}

const char *AdmixOptions::getAlleleFreqScoreFilename2() const
{
  return imp->AlleleFreqScoreFilename2.c_str();
}

const char *AdmixOptions::getAlleleFreqOutputFilename() const
{
  return imp->AlleleFreqOutputFilename.c_str();
}

bool AdmixOptions::getOutputAlleleFreq() const
{
  return imp->OutputAlleleFreq;
}

int AdmixOptions::getAnalysisTypeIndicator() const
{
  return imp->AnalysisTypeIndicator;
}

const char *AdmixOptions::getAssocScoreFilename() const
{
  return imp->AssocScoreFilename.c_str();
}

long AdmixOptions::getSampleEvery() const
{
  return imp->SampleEvery;
}

const char *AdmixOptions::getDispersionTestFilename() const
{
  return imp->DispersionTestFilename.c_str();
}

const char *AdmixOptions::getHistoricalAlleleFreqFilename() const
{
  return imp->HistoricalAlleleFreqFilename.c_str();
}

const char *AdmixOptions::getInputFilename() const
{
  return imp->InputFilename.c_str();
}

const char *AdmixOptions::getMLEFilename() const
{
  return imp->MLEFilename.c_str();
}

bool AdmixOptions::getModelIndicator() const
{
  return imp->modelIndicator;
}

bool AdmixOptions::getIndAdmixHierIndicator() const{
  return imp->IndAdmixHierIndicator;
}
bool AdmixOptions::getMLIndicator()const{
  return imp->MLIndicator;
}

double AdmixOptions::getTruncPt() const
{
  return imp->TruncPt;
}

bool AdmixOptions::getRhoIndicator() const
{
  return imp->RhoIndicator;
}

bool AdmixOptions::getLocusForTestIndicator() const
{
  return imp->locusForTestIndicator;
}

int AdmixOptions::getLocusForTest() const
{
  return imp->LocusForTest;
}

const char *AdmixOptions::getLocusScoreFilename() const
{
  return imp->LocusScoreFilename.c_str();
}

const char *AdmixOptions::getAncestryAssociationScoreFilename() const
{
  return imp->AncestryAssociationScoreFilename.c_str();
}

const char *AdmixOptions::getLogFilename() const
{
  return imp->LogFilename.c_str();
}

int AdmixOptions::getPopulations() const
{
  return imp->Populations;
}

void AdmixOptions::setPopulations(int num)
{
  imp->Populations = num;
}
const char *AdmixOptions::getPriorAlleleFreqFilename() const
{
  return imp->PriorAlleleFreqFilename.c_str();
}

bool AdmixOptions::getScoreTestIndicator() const
{
  return imp->ScoreTestIndicator;
}

long AdmixOptions::getSeed() const
{
  return imp->Seed;
}

double AdmixOptions::getRho() const
{
  return imp->Rho;
}

bool AdmixOptions::getStratificationTest() const
{
  return imp->StratificationTestIndicator;
}

const char *AdmixOptions::getTargetFilename() const
{
  return imp->TargetFilename.c_str();
}

bool AdmixOptions::getTestForAffectedsOnly() const
{
  return imp->TestForAffectedsOnly;
}

bool AdmixOptions::getTestForAllelicAssociation() const
{
  return imp->TestForAllelicAssociation;
}

bool AdmixOptions::getTestForDispersion() const
{
  return imp->TestForDispersion;
}

bool AdmixOptions::getTestForLinkageWithAncestry() const
{
  return imp->TestForLinkageWithAncestry;
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs() const
{
  return imp->TestForMisspecifiedAlleleFreqs;
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs2() const
{
  return imp->TestForMisspecifiedAlleleFreqs2;
}

const char *AdmixOptions::getTestsForSNPsInHaplotypeOutputFilename() const
{
  return imp->TestsForSNPsInHaplotypeOutputFilename.c_str();
}

bool AdmixOptions::getTestForSNPsInHaplotype() const
{
  return imp->TestForSNPsInHaplotype;
}

const char *AdmixOptions::getEtaPriorFilename() const
{
  return imp->EtaPriorFilename.c_str();
}

const char *AdmixOptions::getTestTargetFilename() const
{
  return imp->TestTargetFilename.c_str();
}

int AdmixOptions::getTextIndicator() const
{
  return imp->TextIndicator;
}

int AdmixOptions::sizeInitAlpha() const
{
  return imp->alpha.size();
}

Vector_d AdmixOptions::getInitAlpha(int gamete) const
{
  return imp->alpha[gamete];
}

unsigned int AdmixOptions::IsPedFile() const
{
  return imp->isPedFile;
}

void AdmixOptions::IsPedFile(unsigned int i)
{
  imp->isPedFile = i;
}

unsigned int AdmixOptions::genotypesSexColumn() const
{
  return imp->genotypesSexColumn;
}

void AdmixOptions::genotypesSexColumn(unsigned int i)
{
  imp->genotypesSexColumn = i;
}

bool AdmixOptions::getXOnlyAnalysis() const
{
  return imp->XOnlyAnalysis;
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

    {"paramfile",                             1, 0, 'p'}, // string
    {"regparamfile",                         1, 0, 0}, // string
    {"dispparamfile",                         1, 0, 0}, // string
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

    //optional results directory name option - default is 'results'
    {"resultsdir",                            1, 0,  0 }, //string
   
    // Extra output options
    {"allelicassociationscorefile",           1, 0,  0 }, // string
    {"ancestryassociationscorefile",          1, 0,  0 }, // string
    {"affectedsonlyscorefile",                1, 0,  0 }, // string
    {"indadmixturefile",                      1, 0, 'i'}, // string
    {"ergodicaveragefile",                    1, 0, 'e'}, // string
    {"stratificationtestfile",                1, 0, 'd'}, // string
    {"admixturescorefile",                    1, 0,  0 }, // string
    {"haplotypeassociationscorefile",         1, 0,  0 }, // string
    {"locusfortest",                          1, 0,  0 }, // int 0: no. of composite loci - 1

    // Extra subpopulations options
    {"allelefreqoutputfile",                  1, 0,  0 }, // string
    {"allelefreqscorefile",                   1, 0,  0 }, // string
    {"allelefreqscorefile2",                  1, 0,  0 }, // string
    {"dispersiontestfile",                    1, 0,  0 }, // string
    {"FSToutputfilename",                     1, 0,  0 }, // string

    // Other options
    {"coutindicator",                         1, 0, 'c'}, // int 0: 1
    {"help",                                  0, 0, 'h'}, // NONE
    {"logfile",                               1, 0,  0 }, // string
    {"randommatingmodel",                     1, 0,  0 }, // int 0: 1
    {"globalrho",                             1, 0,  0 }, // int 0: 1
    {"indadmixhiermodel",                     1, 0,  0 }, // int 0: 1
    {"marglikelihood",                        1, 0,  0 }, // int 0: 1
    {"reportedancestry",                      1, 0, 'r'}, // string 
    {"seed",                                  1, 0,  0 }, // long
    {"etapriorfile",                          1, 0,  0 }, // string      
    {"textindicator",                         1, 0,  0 }, // int
    {"sumintensities",                        1, 0,  0 }, // double
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
      {imp->alleleFreqFilename = optarg;OptionValues["allelefreqfile"]=optarg;}
      break;

    case 'b': // burnin
      {imp->burnin = strtol(optarg, NULL, 10);OptionValues["burnin"]=optarg;}
      break;

    case 'c': // coutindicator
      {imp->use_cout = (int)strtol(optarg, NULL, 10);OptionValues["coutindicator"]=optarg;}
      break;

    case 'd': // stratificationtestfile
      imp->DICoutputFilename = optarg;OptionValues["stratificationtestfile"]=optarg;
      imp->StratificationTestIndicator = true;OptionValues["StratificationTestIndicator"]="1";
      break;

    case 'e': // ergodicaveragefile
      {imp->ErgodicAverageFilename = optarg;OptionValues["ergodicaveragefile"]=optarg;}
      break;

    case 'g': // genotypesfile
      {imp->GeneticDataFilename = optarg;OptionValues["genotypesfile"]=optarg;}
      break;

    case 'h': // help
      cout << "Usage: " << args[0] << " [options]" << endl;
      cout << " some help text goes here" << endl;
//outputlist of user options to console
      exit(0);

    case 'i': // indadmixturefile
      {imp->IndAdmixtureFilename = optarg;OptionValues["indadmixturefile"]=optarg;}
      break;

    case 'l': // locusfile
      {imp->GeneInfoFilename = optarg;OptionValues["locusfile"]=optarg;}
      break;

    case 'p': // paramfile
      {imp->ParameterFilename = optarg;OptionValues["paramfile"]=optarg;}
      break;

    case 'r': // reportedancestry
      {imp->ReportedAncestryFilename = optarg;OptionValues["reportedancestry"]=optarg;}
      break;

    case 's': // samples
      {imp->TotalSamples = strtol(optarg, NULL, 10);OptionValues["samples"]=optarg;}
      break;

    case 't': // targetindicator
      {imp->TargetIndicator = (int)strtol(optarg, NULL, 10);OptionValues["targetindicator"]=optarg;}
      break;

    case '?':
      exit(1);

    case 0:
      if(long_option_name == "resultsdir"){
	imp->ResultsDir = optarg;OptionValues["resultsdir"]=optarg;
      }else if (long_option_name == "regparamfile") {
	imp->RegressionOutputFilename = optarg;OptionValues["regparamfile"]=optarg;
      }else if (long_option_name == "dispparamfile") {
	imp->EtaOutputFilename = optarg; OptionValues["dispparamfile"]=optarg;
      }else  if (long_option_name == "admixturescorefile") {
	imp->AssocScoreFilename = optarg;OptionValues["admixturescorefile"]=optarg;
	imp->ScoreTestIndicator = true;OptionValues["ScoreTestIndicator"]="1";
      } else if (long_option_name == "affectedsonlyscorefile") {
	imp->AffectedsOnlyScoreFilename = optarg;OptionValues["affectedsonlyscorefile"]=optarg;
	imp->TestForAffectedsOnly = true;OptionValues["TestForAffectedsOnly"]="1";
      } else if (long_option_name == "allelicassociationscorefile") {
	imp->LocusScoreFilename = optarg;OptionValues["allelicassociationscorefile"]=optarg;
	imp->TestForAllelicAssociation = true;OptionValues["TestForAllelicAssociation"]="1";
      } else if (long_option_name == "allelefreqoutputfile") {
	imp->AlleleFreqOutputFilename = optarg;OptionValues["allelefreqoutputfile"]=optarg;
	imp->OutputAlleleFreq = true;OptionValues["OutputAlleleFreq"]="1";
      } else if (long_option_name == "allelefreqscorefile2") {
	imp->AlleleFreqScoreFilename = optarg;OptionValues["allelefreqscorefile2"]=optarg;
	imp->TestForMisspecifiedAlleleFreqs = true;OptionValues["TestForMisspecifiedAlleleFreqs2"]="1";
      } else if (long_option_name == "allelefreqscorefile"){
	imp->AlleleFreqScoreFilename2 = optarg;OptionValues["allelefreqscorefile"]=optarg;
	imp->TestForMisspecifiedAlleleFreqs2 = true;OptionValues["TestForMisspecifiedAlleleFreqs"]="1";
      } else if (long_option_name == "analysistypeindicator") {
	imp->AnalysisTypeIndicator = (int)strtol(optarg, NULL, 10);OptionValues["analysistypeindicator"]=optarg;
      } else if (long_option_name == "covariatesfile") {
	imp->InputFilename = optarg;OptionValues["covariatesfile"]=optarg;
      } else if (long_option_name == "mlefile") {
	imp->MLEFilename = optarg;OptionValues["mlefile"]=optarg;
      } else if (long_option_name == "dispersiontestfile") {
	imp->DispersionTestFilename = optarg;OptionValues["dispersiontestfile"]=optarg;
	imp->TestForDispersion = true;OptionValues["TestForDispersion"]="1";
      } else if (long_option_name == "every") {
	imp->SampleEvery = strtol(optarg, NULL, 10);OptionValues["every"]=optarg;
      } else if (long_option_name == "FSToutputfilename") {
	imp->FSTOutputFilename = optarg;OptionValues["FSToutputfilename"]=optarg;
	imp->OutputFST = true;OptionValues["OutputFST"]="1";
      } else if (long_option_name == "fixedallelefreqs") {
	if (strtol(optarg, NULL, 10) == 1) {
	  imp->fixedallelefreqs = true;OptionValues["fixedallelefreqs"]="1";
	}
      } else if (long_option_name == "xonlyanalysis") {
	if (strtol(optarg, NULL, 10) == 1) {
	  imp->XOnlyAnalysis = true;OptionValues["xonlyanalysis"]="1";
	}
      } else if (long_option_name == "historicallelefreqfile") {
	imp->HistoricalAlleleFreqFilename = optarg;OptionValues["historicallelefreqfile"]=optarg;
      } else if (long_option_name == "locusfortest") {
	imp->LocusForTest = (int)strtol(optarg, NULL, 10);OptionValues["locusfortest"]=optarg;
	imp->locusForTestIndicator = true;OptionValues["locusForTestIndicator"]="1";
      } else if (long_option_name == "ancestryassociationscorefile") {
	imp->AncestryAssociationScoreFilename = optarg;OptionValues["ancestryassociationscorefile"]=optarg;
	imp->TestForLinkageWithAncestry = true;OptionValues["TestForLinkageWithAncestry"]="1";
      } else if (long_option_name == "logfile") {
	imp->LogFilename = optarg;OptionValues["logfile"]=optarg;
      } else if (long_option_name == "randommatingmodel") {
	if (strtol(optarg, NULL, 10) == 1) {
	  imp->modelIndicator = true;OptionValues["randommatingmodel"]="1";
	}
      } else if (long_option_name == "indadmixhiermodel") {
	if (strtol(optarg, NULL, 10) == 0) {
	  imp->IndAdmixHierIndicator = false;OptionValues["indadmixhiermodel"]="0";
	}
      }else if (long_option_name == "marglikelihood") {
	if (strtol(optarg, NULL, 10) == 1) {
	  imp->MLIndicator = true;OptionValues["marglikelihood"]="1";
	}
      }else if (long_option_name == "globalrho") {
	if (strtol(optarg, NULL, 10) == 1) {
	  imp->RhoIndicator = false;OptionValues["globalrho"]="0";
	} else if (strtol(optarg, NULL, 10) == 0) {
	  imp->RhoIndicator = true;OptionValues["globalrho"]="1";
	} else {
	  cerr << "Set global rho to 0 or 1.\n";
	  exit(1);
	}
      } else if (long_option_name == "outcomevarfile") {
	imp->TargetFilename = optarg;OptionValues["outcomevarfile"]=optarg;
      } else if (long_option_name == "populations") {
	setPopulations((int)strtol(optarg, NULL, 10));OptionValues["populations"]=optarg;
      } else if (long_option_name == "priorallelefreqfile") {
	imp->PriorAlleleFreqFilename = optarg;OptionValues["priorallelefreqfile"]=optarg;
      } else if (long_option_name == "scoretestindicator") {
	// TODO:
    OptionValues["ScoreTestIndicator"]=optarg;
      } else if (long_option_name == "seed") {
	imp->Seed = strtol(optarg, NULL, 10);OptionValues["seed"]=optarg;
      } else if (long_option_name == "sumintensities") {
	imp->Rho = strtod(optarg, NULL);OptionValues["sumintensities"]=optarg;
      } else if (long_option_name == "truncationpoint") {
	imp->TruncPt = strtod(optarg, NULL);OptionValues["truncationpoint"]=optarg;
      } else if (long_option_name == "haplotypeassociationscorefile") {
	imp->TestsForSNPsInHaplotypeOutputFilename = optarg;OptionValues["haplotypeassociationscorefile"]=optarg;
	imp->TestForSNPsInHaplotype = true;OptionValues["TestForSNPsInHaplotype"]=optarg;
      } else if (long_option_name == "etapriorfile") {
	imp->EtaPriorFilename = optarg;OptionValues["etapriorfile"]=optarg;
      } else if (long_option_name == "initalpha0" ) {
	imp->alpha.push_back(CstrToVec2(optarg));OptionValues["initalpha0"]=optarg;
      } else if (long_option_name == "initalpha1") {
	imp->alpha.push_back(CstrToVec2(optarg));OptionValues["initalpha1"]=optarg;
      } else if (long_option_name == "textindicator") {
	imp->TextIndicator = (int)strtol(optarg, NULL, 10);OptionValues["textindicator"]=optarg;
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
  if (imp->LogFilename != "")imp->LogFilename = imp->ResultsDir + "/" + imp->LogFilename; 
  if (imp->ParameterFilename != "")imp->ParameterFilename = imp->ResultsDir + "/" + imp->ParameterFilename; 
  if (imp->RegressionOutputFilename != "")imp->RegressionOutputFilename = imp->ResultsDir + "/" + imp->RegressionOutputFilename; 
  if (imp->EtaOutputFilename != "")imp->EtaOutputFilename = imp->ResultsDir + "/" + imp->EtaOutputFilename; 
  if (imp->LocusScoreFilename != "")imp->LocusScoreFilename = imp->ResultsDir + "/" + imp->LocusScoreFilename; 
  if (imp->AncestryAssociationScoreFilename != "")imp->AncestryAssociationScoreFilename = imp->ResultsDir + "/" + imp->AncestryAssociationScoreFilename; 
  if (imp->AffectedsOnlyScoreFilename != "")imp->AffectedsOnlyScoreFilename = imp->ResultsDir + "/" + imp->AffectedsOnlyScoreFilename; 
  if (imp->IndAdmixtureFilename != "")imp->IndAdmixtureFilename = imp->ResultsDir + "/" + imp->IndAdmixtureFilename; 
  if (imp->ErgodicAverageFilename != "")imp->ErgodicAverageFilename = imp->ResultsDir + "/" + imp->ErgodicAverageFilename; 
  if (imp->DICoutputFilename != "")imp->DICoutputFilename = imp->ResultsDir + "/" + imp->DICoutputFilename; 
  if (imp->AssocScoreFilename != "")imp->AssocScoreFilename = imp->ResultsDir + "/" + imp->AssocScoreFilename; 
  if (imp->TestsForSNPsInHaplotypeOutputFilename != "")imp->TestsForSNPsInHaplotypeOutputFilename = imp->ResultsDir + "/" + imp->TestsForSNPsInHaplotypeOutputFilename; 
  if (imp->AlleleFreqOutputFilename != "")imp->AlleleFreqOutputFilename = imp->ResultsDir + "/" + imp->AlleleFreqOutputFilename; 
  if (imp->AlleleFreqScoreFilename2 != "")imp->AlleleFreqScoreFilename2 = imp->ResultsDir + "/" + imp->AlleleFreqScoreFilename2; 
  if (imp->AlleleFreqScoreFilename != "")imp->AlleleFreqScoreFilename = imp->ResultsDir + "/" + imp->AlleleFreqScoreFilename; 
  if (imp->DispersionTestFilename != "")imp->DispersionTestFilename = imp->ResultsDir + "/" + imp->DispersionTestFilename; 
  if (imp->FSTOutputFilename != "")imp->FSTOutputFilename = imp->ResultsDir + "/" + imp->FSTOutputFilename; 
  
}

void AdmixOptions::PrintOptions(){
  //set populations value in case it has changed
  //NB do similar for any option that can be changed outside AdmixOptions
  std::ostringstream s;
  if (s << imp->Populations) // conversion worked
    {
    OptionValues["populations"] = (char *)s.str().c_str();
    }

  //Now output Options table to args.txt
  string ss;
  ss = imp->ResultsDir + "/args.txt";
  ofstream argstream(ss.c_str());

  for( OptionMap::iterator p= OptionValues.begin(); p!=OptionValues.end(); p++) {
    argstream << (*p).first << "=" << (*p).second <<endl;
  }
  argstream.close();
}

int AdmixOptions::checkOptions(LogWriter *Log){
  if(getScoreTestIndicator() &&
      ( getTestForLinkageWithAncestry() || getTestForAllelicAssociation() ) ){
    Log->logmsg(true,"Cannot test for linkage with ancestry or allelic association\n");
    Log->logmsg(true,"with score test for association. Can only use affecteds only test\n");
    Log->logmsg(true,"for linkage.\n");
    Log->logmsg(true,"If admixturescorefile is selected, then please unselect both\n");
    Log->logmsg(true,"allelicassociationscorefile and ancestryassociationscorefile\n");
    exit(1);
  }
  //the rest taken from Latent constructor
  if( getTestForMisspecifiedAlleleFreqs() &&
      ( !strlen( getAlleleFreqFilename() ) && !(getFixedAlleleFreqs()) ) ){
    Log->logmsg(true,"Cannot test for mis-specified allele frequencies with unknown allele frequncies.\n");
    exit(0);
  }
  
  if (getAnalysisTypeIndicator() == 0)
    {
      Log->logmsg(true,"Affecteds only analysis.\n");
      if( !getTestForAffectedsOnly() ){
        Log->logmsg(true,"Must specify affectedsonlyscorefile.\n");
      }
    }
  else if (getAnalysisTypeIndicator() == 1)
    {
      Log->logmsg(true,"Cross sectional analysis, no outcome.\n");
    }
  else if (getAnalysisTypeIndicator() == 2)
    {
      Log->logmsg(true,"Cross sectional analysis, continuous outcome.\n");
      if( strlen(getTargetFilename() ) == 0 )
	{
	   Log->logmsg(true,"Must specify target filename.\n");
	   exit(0);
	}
    }
  else if (getAnalysisTypeIndicator() == 3)
    {
      Log->logmsg(true,"Cross sectional analysis, binary outcome.\n");
      if( strlen(getTargetFilename() ) == 0 )
	{
	  Log->logmsg(true,"Must specify target filename.\n");
	  exit(0);
	}
    }
  else if (getAnalysisTypeIndicator() == 4)
    {
      Log->logmsg(true,"Case control analysis.\n");
      if( strlen(getTargetFilename() ) == 0  )
	{
	   Log->logmsg(true,"Must specify target filename.\n");
	   exit(0);
	}
    }
  else if (getAnalysisTypeIndicator() == 5)
    {
      Log->logmsg(true,"Cross sectional analysis, multiple outcome.\n");
      if( strlen(getTargetFilename() ) == 0  )
	{
	  Log->logmsg(true,"Must specify target filename.\n");
	  exit(0);
	}
    }
  else if (getAnalysisTypeIndicator() == -1 || getAnalysisTypeIndicator() == -2)
    {
      Log->logmsg(true,"One individual analysis");
      if(getMLIndicator())Log->logmsg(true, " with marginal likelihood calculation");
      Log->logmsg(true, "\n");
    }

  else
    {
      Log->logmsg(true,"Unknown analysis type: ");
      Log->logmsg(true, getAnalysisTypeIndicator());
      Log->logmsg(true, "\n");
      exit(0);
    }


  if (!getIndAdmixHierIndicator())
    {
      Log->logmsg(true,"No hierarchical model for individuals.\n");
    }
 
  if(getModelIndicator() )
    Log->logmsg(true,"Model assuming random mating.\n");
  else 
    Log->logmsg(true,"Model assuming assortative mating.\n");

  if(getMLIndicator()){
    if(getAnalysisTypeIndicator() >= 0){
      Log->logmsg(true, "Error: Cannot calculate marginal likelihood with analysis type ");
      Log->logmsg(true, getAnalysisTypeIndicator());
      Log->logmsg(true, "\n");
      exit(0);
    }
    //change this when marginal likelihood can be calculated for other type of model
    //else Log->logmsg(true,"Analysis with marginal likelihood calculation\n");
  }
  // Check whether genotypes file has been specified
  if ( strlen(getGeneticDataFilename() ) == 0 )
    {
      Log->logmsg(true,"Must specify geneticdata filename.\n");
      exit( 1 );
    }
  if(getPopulations() > 0 )
    {
      Log->logmsg(true,"No allelefreq filename or priorallelefreq filename given.\n");
      Log->logmsg(true,"Default priors will be set for the allele frequencies with ");
      Log->logmsg(true,getPopulations());
      Log->logmsg(true," populations.\n");
    }
  else if( strlen( getAlleleFreqFilename() ) ||
           (strlen( getPriorAlleleFreqFilename() ) && getFixedAlleleFreqs() ) ){
    Log->logmsg(true,"Analysis with fixed allele frequencies.\n");
  }
  else if( strlen(getPriorAlleleFreqFilename() ) && !getFixedAlleleFreqs() ){
    Log->logmsg(true,"Analysis with prior allele frequencies.\n");
  }
  else if( strlen( getHistoricalAlleleFreqFilename() ) ){
    Log->logmsg(true,"Analysis with dispersion model for allele frequencies.\n");
  }
  if ( strlen( getGeneInfoFilename() ) == 0 )
    {
      Log->logmsg(true,"Must specify locusinfo filename.\n");
      exit( 1 );
    }

  if( getTestForLinkageWithAncestry() && getPopulations() == 1 ){
    Log->logmsg(true,"Cannot test for linkage with ancestry with 1 population.\n");
    exit(0);
  }
  return 1;
}
