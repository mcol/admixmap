/** 
 *   ADMIXMAP
 *   AdmixOptions.cc 
 *   Class to hold program options
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AdmixOptions.h"
//#include "utils/StringConvertor.h"
#include <string.h>
#include <sstream>
#include <numeric> // for checkInitAlpha
#include "defs.h"

using namespace std;

// AdmixOptions::AdmixOptions()
// {
//   Initialise();
// }
AdmixOptions::AdmixOptions(int argc,  char** argv){
  //base class sets default base options
  //this call to SetDefaultValues calls this class' options
  SetDefaultValues();
  ReadUserOptions(argc, argv);
  OptionMap ProgOptions;
  SetOptions(ProgOptions);
}

void AdmixOptions::SetDefaultValues(){
  Populations = 1;
  OutputFST = false;
  locusForTestIndicator = false;
  LocusForTest = -1;
  FreqDispersionHierModel = false;
  correlatedallelefreqs = false;
  RandomMatingModel = false;
  GlobalRho = true;//corresponds to globalrho = 1;
  PopAdmixPropsAreEqual = false;
  IndAdmixHierIndicator = true; //hierarchical model on ind admixture
  HapMixModelIndicator = false; //model haplotypes with mixture model
  chibIndicator = false;//calculate marginal likelihood by Chib method
  TestOneIndivIndicator = false; // evaluate marginal likelihood for single individual
  TestForAdmixtureAssociation = false;
  StratificationTestIndicator = false;
  TestForAffectedsOnly = false;
  TestForHaplotypeAssociation = false;
  TestForDispersion = false;
  TestForLinkageWithAncestry = false;
  TestForMisspecifiedAlleleFreqs = false;
  TestForMisspecifiedAlleleFreqs2 = false;
  HWTest = false;

//   globalrhoPrior.push_back(3.0);//rhoalpha 
//   globalrhoPrior.push_back(0.5);//rhobeta

//   rhoPrior.push_back(6.0);//rhoalpha 
//   rhoPrior.push_back(5.0);//rhobeta shape
//   rhoPrior.push_back(4.0);//rhobeta rate

  initalpha.resize(2);//TODO: check if this is necessary
  //gamma(3, 0.01) prior on dispersion parameter
  etamean = 100.0; 
  etavar = 2500.0; 


  //prior on frequency Dirichlet prior params in hapmixmodel
  //allelefreqprior.push_back(0.01);//shape
  //allelefreqprior.push_back(0.1);//prior shape of rate
  //allelefreqprior.push_back(1.0);//prior rate of rate

  LikRatioFilename = "LikRatioFile.txt";//hardcoding for now, can change later

  // option names and default option values are stored as strings in a map container 
  // these are default values
  // other specified options will be appended to this array 
#ifdef __HAPMIXMAP__
  //final values are always output to file, with these default filenames
  FinalFreqPriorFilename = "finalfreqpriors.txt";
  FinalLambdaFilename = "finallambdas.txt";
  AlleleFreqOutputFilename = "finalallelefreqs.txt";

  useroptions["hapmixmodel"] = "1";
  useroptions["hapmixlambdaprior"] = "30, 0.1, 10, 1";
  useroptions["finalfreqpriorfile"] = "finalfreqpriors.txt";
  useroptions["finallambdafile"] = "finallambdas.txt";
#else
  useroptions["hapmixmodel"] = "0";
  useroptions["correlatedallelefreqs"] = "0";
  useroptions["randommatingmodel"] = "0";
  useroptions["globalrho"] = "1";
  useroptions["indadmixhiermodel"] = "1";
  useroptions["chib"] = "0";
#endif
  //global rho: default gamma (3, 0.5) prior has mean 6, variance 12 
  useroptions["globalsumintensitiesprior"] = "3.0,0.5";
  // non-global rho: default gamma-gamma prior with parameters n=6, alpha=5, beta=4
  // effective prior mean is 6*4/(5-1) = 6 and effective prior variance is 6*7 / (5-2) = 14
  useroptions["sumintensitiesprior"] = "6.0,5.0,4.0";
}

AdmixOptions::~AdmixOptions()
{
}

// each option has a function to return its value
const char *AdmixOptions::getStratTestFilename() const
{
  return StratTestFilename.c_str();
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
  return (fixedallelefreqs || (alleleFreqFilename.length() && !HapMixModelIndicator));
}
bool AdmixOptions::getCorrelatedAlleleFreqs() const
{
  return correlatedallelefreqs;
}

const char *AdmixOptions::getIndAdmixtureFilename() const
{
  return IndAdmixtureFilename.c_str();
}

const char *AdmixOptions::getEtaOutputFilename() const
{
  return EtaOutputFilename.c_str();
}

const char *AdmixOptions::getReportedAncestryFilename() const
{
  return ReportedAncestryFilename.c_str();
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

const char *AdmixOptions::getAlleleFreqPriorOutputFilename() const
{
  return AlleleFreqPriorOutputFilename.c_str();
}
bool AdmixOptions::OutputAlleleFreqPrior()const{
  return (bool)(AlleleFreqPriorOutputFilename.size()>0);
}

const char *AdmixOptions::getAssocScoreFilename() const
{
  return AssocScoreFilename.c_str();
}

const char *AdmixOptions::getDispersionTestFilename() const
{
  return DispersionTestFilename.c_str();
}

bool AdmixOptions::getMHTest()const{
  return (bool)(MHTestFilename.size());
}
const char* AdmixOptions::getMHTestFilename()const{
  return MHTestFilename.c_str();
}

const char *AdmixOptions::getHistoricalAlleleFreqFilename() const
{
  return HistoricalAlleleFreqFilename.c_str();
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

bool AdmixOptions::getTestOneIndivIndicator()const{
  return TestOneIndivIndicator;
}

bool AdmixOptions::isGlobalRho() const
{
  return GlobalRho;
}

bool AdmixOptions::isFreqDispersionHierModel()const{
  return FreqDispersionHierModel;
}
bool AdmixOptions::getLocusForTestIndicator() const
{
  return locusForTestIndicator;
}

int AdmixOptions::getLocusForTest() const
{
  return LocusForTest;
}

const char *AdmixOptions::getAncestryAssociationScoreFilename() const
{
  return AncestryAssociationScoreFilename.c_str();
}

int AdmixOptions::getPopulations() const
{
  return Populations;
}

void AdmixOptions::setPopulations(int num)
{
  Populations = num;
}

bool AdmixOptions::getTestForAdmixtureAssociation() const
{
  return TestForAdmixtureAssociation;
}

double AdmixOptions::getRhoPriorMean()const{
  if( !HapMixModelIndicator && (GlobalRho || !IndAdmixHierIndicator ) )
    return globalrhoPrior[0] / globalrhoPrior[1];
  else 
    return rhoPrior[0] * rhoPrior[2] / (rhoPrior[1] - 1.0);
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
const std::vector<double> &AdmixOptions::getHapMixLambdaPrior()const{
  return hapmixlambdaprior;
}
double AdmixOptions::getEtaMean() const{
  return etamean;
}
double AdmixOptions::getEtaVar() const{
  return etavar;
}

bool AdmixOptions::getStratificationTest() const
{
  return StratificationTestIndicator;
}

void AdmixOptions::setStratificationTest(bool b){
  StratificationTestIndicator = b;
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

const char* AdmixOptions::getInitialHapMixLambdaFilename()const{
    return InitialHapMixLambdaFilename.c_str();
}
const char* AdmixOptions::getInitialFreqPriorFilename()const{
  return InitialFreqPriorFile.c_str();
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
bool AdmixOptions::PopAdmixturePropsAreEqual()const{
  return PopAdmixPropsAreEqual;
}
const vector<float>& AdmixOptions::getrhoSamplerParams()const{
  return rhoSamplerParams;
}
const vector<float>& AdmixOptions::getPopAdmixSamplerParams()const{
  return popAdmixSamplerParams;
}

const std::vector<double> & AdmixOptions::getAlleleFreqPriorParams()const{
  return allelefreqprior;
}

const char* AdmixOptions::getHapMixLambdaOutputFilename()const{
    return HapMixLambdaOutputFilename.c_str();
}
const char* AdmixOptions::getCCGenotypesFilename()const{
  return CCGenotypesFilename.c_str();
}

const vector<unsigned>& AdmixOptions::getMaskedIndividuals()const{
  return MaskedIndividuals;
}
const vector<unsigned>& AdmixOptions::getMaskedLoci()const{
  return MaskedLoci;
}
bool AdmixOptions::OutputCGProbs()const{
  return (MaskedLoci.size() && MaskedIndividuals.size()); 
}

const char* AdmixOptions::getFinalFreqPriorFilename()const{
  return FinalFreqPriorFilename.c_str();
}
const char* AdmixOptions::getFinalLambdaFilename()const{
  return FinalLambdaFilename.c_str();
}

unsigned AdmixOptions::GetNumMaskedIndividuals()const{
  return MaskedIndividuals.size();
}
unsigned AdmixOptions::GetNumMaskedLoci()const{
  return MaskedLoci.size();
}

void AdmixOptions::SetOptions(OptionMap& ProgOptions)
{
  //set up Option map
  //first set options for this class
  //A = AdmixOption, H = HapMixOption, C = Common (Base class)

  ProgOptions["populations"] = OptionPair(&Populations, "int");//A
  ProgOptions["states"] = OptionPair(&Populations, "int");//H
  ProgOptions["historicallelefreqfile"] = OptionPair(&HistoricalAlleleFreqFilename, "string");//A
  ProgOptions["reportedancestry"] = OptionPair(&ReportedAncestryFilename, "string");//A
  ProgOptions["ccgenotypesfile"] = OptionPair(&CCGenotypesFilename, "string");//H
  //standard output files (optional)

  //file to write sampled values of dispersion parameter (ADMIXMAP) or mean and variance of dispersion params (HAPMIXMAP)
  ProgOptions["dispparamfile"] = OptionPair(&EtaOutputFilename, "outputfile");// C
  ProgOptions["indadmixturefile"] = OptionPair(&IndAdmixtureFilename, "outputfile");//A

  ProgOptions["finalfreqpriorfile"] = OptionPair(&FinalFreqPriorFilename, "outputfile");//H
  ProgOptions["finallambdafile"] = OptionPair(&FinalLambdaFilename, "outputfile");//H
  ProgOptions["allelefreqprioroutputfile"] = OptionPair(&AlleleFreqPriorOutputFilename, "outputfile");//H
  ProgOptions["hapmixlambdaoutputfile"] = OptionPair(&HapMixLambdaOutputFilename, "outputfile");//H, remove hapmix from name
  //prior and model specification
  ProgOptions["randommatingmodel"] = OptionPair(&RandomMatingModel, "bool");//A
  ProgOptions["globalrho"] = OptionPair(&GlobalRho, "bool");//A
  ProgOptions["indadmixhiermodel"] = OptionPair(&IndAdmixHierIndicator, "bool");//A
  ProgOptions["hapmixmodel"] = OptionPair(&HapMixModelIndicator, "bool");
  ProgOptions["etapriorfile"] = OptionPair(&EtaPriorFilename, "string");//C
  ProgOptions["globalsumintensitiesprior"] = OptionPair(&globalrhoPrior, "dvector");//A
  ProgOptions["sumintensitiesprior"] = OptionPair(&rhoPrior, "dvector");//A
  ProgOptions["hapmixlambdaprior"] = OptionPair(&hapmixlambdaprior, "dvector");//H
  ProgOptions["allelefreqprior"] = OptionPair(&allelefreqprior, "dvector");//H
  ProgOptions["etapriormean"] = OptionPair(&etamean, "double");//A
  ProgOptions["etapriorvar"] = OptionPair(&etavar, "double");//A
  ProgOptions["admixtureprior"] = OptionPair(&initalpha[0], "dvector");//A
  ProgOptions["admixtureprior1"] = OptionPair(&initalpha[1], "dvector");//A
  ProgOptions["correlatedallelefreqs"] = OptionPair(&correlatedallelefreqs, "bool");//A
  ProgOptions["popadmixproportionsequal"] = OptionPair(&PopAdmixPropsAreEqual, "bool");//A
  ProgOptions["initialhapmixlambdafile"] = OptionPair(&InitialHapMixLambdaFilename, "string");//H, remove hapmix from name
  ProgOptions["initialfreqpriorfile"] = OptionPair(&InitialFreqPriorFile, "string");//H
  ProgOptions["freqdispersionhiermodel"] = OptionPair(&FreqDispersionHierModel, "bool");//H
  //sampler settings
  ProgOptions["rhosamplerparams"] = OptionPair(&rhoSamplerParams, "fvector");//A
  ProgOptions["popadmixsamplerparams"] = OptionPair(&popAdmixSamplerParams, "fvector");//A
  // test options
  ProgOptions["ancestryassociationscorefile"] = OptionPair(&AncestryAssociationScoreFilename, "outputfile");//C
  ProgOptions["affectedsonlyscorefile"] = OptionPair(&AffectedsOnlyScoreFilename, "outputfile");//A
  ProgOptions["admixturescorefile"] = OptionPair(&AssocScoreFilename, "outputfile");//A
  ProgOptions["haplotypeassociationscorefile"] = OptionPair(&HaplotypeAssociationScoreFilename, "outputfile");//A
  ProgOptions["stratificationtestfile"] = OptionPair(&StratTestFilename, "outputfile");//A
  ProgOptions["allelefreqscorefile"] = OptionPair(&AlleleFreqScoreFilename, "outputfile");//A
  ProgOptions["allelefreqscorefile2"] = OptionPair(&AlleleFreqScoreFilename2, "outputfile");//A
  ProgOptions["dispersiontestfile"] = OptionPair(&DispersionTestFilename, "outputfile");//A
  ProgOptions["fstoutputfile"] = OptionPair(&FSTOutputFilename, "outputfile");//A
  ProgOptions["mhscoretestfile"] = OptionPair(&MHTestFilename, "outputfile");//H
  ProgOptions["likratiofile"] = OptionPair(&LikRatioFilename, "outputfile");//A, requires affectedsonlytestfile
  ProgOptions["indadmixmodefile"] = OptionPair(&IndAdmixModeFilename, "outputfile");//A
  ProgOptions["testgenotypesfile"] = OptionPair(0, "null");//A
  ProgOptions["locusfortest"] = OptionPair(&LocusForTest, "int");//A
  ProgOptions["maskedindivs"] = OptionPair(&MaskedIndividuals, "range");//H
  ProgOptions["maskedloci"] = OptionPair(&MaskedLoci, "range");//H
  // Other options
  ProgOptions["chib"] = OptionPair(&chibIndicator, "bool");//A,  Marginal likelihood by Chib algo
  ProgOptions["testoneindiv"] = OptionPair(&TestOneIndivIndicator, "bool");//A,  ML for one individual in a collection 
  //old options - do nothing but kept for backward-compatibility with old scripts
  ProgOptions["analysistypeindicator"] = OptionPair(0, "old");//A
  ProgOptions["coutindicator"] = OptionPair(0, "old");//A
  ProgOptions["truncationpoint"] = OptionPair(0, "old");//A

  //now set base options and finish parsing
  Options::SetOptions(ProgOptions);

  //set indicators
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

int AdmixOptions::checkOptions(LogWriter &Log, int NumberOfIndividuals){
  //check base options
  Options::checkOptions(Log, 1);

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

    if(useroptions["allelefreqprior"].size()){

      if(FreqDispersionHierModel && (allelefreqprior.size() !=3)) {
        Log << "Error: 'allelefreqprior' must have length 3\n";
        badOptions = true;
      }
      else if(allelefreqprior.size()< 2){
        Log << "Error: 'allelefreqprior' must have length 2\n";
        badOptions = true;
      }
    }

    useroptions.erase("sumintensitiesprior") ;
    useroptions.erase("globalsumintensitiesprior") ;  
  }
  else{
    useroptions.erase("hapmixlambdaprior") ;
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
  if( (alleleFreqFilename.length() && !HapMixModelIndicator) ||
           (PriorAlleleFreqFilename.length() && fixedallelefreqs ) ){
    Log << "Analysis with fixed allele frequencies.\n";
    if(OutputAlleleFreq && !HapMixModelIndicator){
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

