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
#include "AdmixFilenames.h"
#include <string.h>
#include <sstream>
#include <iterator>
#include <numeric> // for checkInitAlpha

using namespace std;
using namespace bclib;

AdmixOptions::AdmixOptions(int argc,  char** argv){
  SetDefaultValues();
  DefineOptions();
  if(!ReadUserOptions(argc, argv, "-f"))
    exit(1);
}

void AdmixOptions::SetDefaultValues(){
  Populations = 1;
  OutputFST = false;
  locusForTestIndicator = false;
  LocusForTest = -1;
  correlatedallelefreqs = false;
  RandomMatingModel = false;
  GlobalRho = true;//corresponds to globalrho = 1;
  PopAdmixPropsAreEqual = false;
  IndAdmixHierIndicator = true; //hierarchical model on ind admixture
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
  ScoreTestIndicator = false; //indicator for any of the score tests in ScoreTests class
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

  // option names and default option values are stored as strings in a map container 
  // these are default values
  // other specified options will be appended to this array 

  useroptions["correlatedallelefreqs"] = "0";
  useroptions["randommatingmodel"] = "0";
  useroptions["globalrho"] = "1";
  useroptions["indadmixhiermodel"] = "1";
  useroptions["chib"] = "0";

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
bool AdmixOptions::outputParams()const{
  return ( ParameterFilename.length() || RegressionOutputFilename.length() || EtaOutputFilename.length());
}

bool AdmixOptions::getOutputFST() const{
  return OutputFST;
}

bool AdmixOptions::getScoreTestIndicator() const{
  return ScoreTestIndicator;
}

bool AdmixOptions::getTestForAdmixtureAssociation() const{
  return TestForAdmixtureAssociation;
}
bool AdmixOptions::getStratificationTest() const{
  return StratificationTestIndicator;
}

void AdmixOptions::setStratificationTest(bool b){
  StratificationTestIndicator = b;
}

bool AdmixOptions::getTestForAffectedsOnly() const{
  return TestForAffectedsOnly;
}

void AdmixOptions::setTestForAffectedsOnly(bool b){
  TestForAffectedsOnly = b;
  if(!b)
    useroptions.erase("affectedsonlytest");
}

bool AdmixOptions::getTestForDispersion() const{
  return TestForDispersion;
}

bool AdmixOptions::getTestForLinkageWithAncestry() const{
  return TestForLinkageWithAncestry;
}

void AdmixOptions::setTestForLinkageWithAncestry(bool b){
  TestForLinkageWithAncestry = b;
  if(!b)
    useroptions.erase("ancestryassociationtest");
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs() const{
  return TestForMisspecifiedAlleleFreqs;
}

bool AdmixOptions::getTestForMisspecifiedAlleleFreqs2() const{
  return TestForMisspecifiedAlleleFreqs2;
}

bool AdmixOptions::getTestForHaplotypeAssociation() const{
  return TestForHaplotypeAssociation;
}

void AdmixOptions::setTestForHaplotypeAssociation(bool b){
  TestForHaplotypeAssociation = b;
  if(!b)useroptions.erase("haplotypeassociationtest");
}

bool AdmixOptions::getFixedAlleleFreqs() const
{
  return (fixedallelefreqs || alleleFreqFilename.length() );
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

const char *AdmixOptions::getAlleleFreqFilename() const
{
  return alleleFreqFilename.c_str();
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

bool AdmixOptions::getLocusForTestIndicator() const
{
  return locusForTestIndicator;
}

int AdmixOptions::getLocusForTest() const
{
  return LocusForTest;
}

int AdmixOptions::getPopulations() const
{
  return Populations;
}

void AdmixOptions::setPopulations(int num)
{
  Populations = num;
}

const string& AdmixOptions::getPopLabelString()const{
  return PopLabels;
}

double AdmixOptions::getRhoPriorMean()const{
  if( GlobalRho || !IndAdmixHierIndicator  )
    return globalrhoPrior[0] / globalrhoPrior[1];
  else 
    return rhoPrior[0] * rhoPrior[2] / (rhoPrior[1] - 1.0);
}

double AdmixOptions::getRhoalpha() const {
  if(GlobalRho || !IndAdmixHierIndicator) {
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

double AdmixOptions::getEtaMean() const{
  return etamean;
}
double AdmixOptions::getEtaVar() const{
  return etavar;
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

bool AdmixOptions::isSymmetric()const{
  return _symmetric;
}
bool AdmixOptions::isAdmixed(unsigned gamete)const{
  return _admixed[gamete];
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

void AdmixOptions::DefineOptions(){
  //set up Option map

  addOption("populations", intOption, &Populations);
  addOption("allelefreqfile", stringOption, &alleleFreqFilename);
  addOption("historicallelefreqfile", stringOption, &HistoricalAlleleFreqFilename);
  addOption("poplabels", stringOption, &PopLabels);
  addOption("reportedancestry", stringOption, &ReportedAncestryFilename);
  //standard output files (optional)

  //file to write sampled values of dispersion parameter
  addOption("dispparamfile", outputfileOption, &EtaOutputFilename);
  addOption("indadmixturefile", outputfileOption, &IndAdmixtureFilename);

  //prior and model specification
  addOption("hapmixmodel", nullOption, 0);//this will cause program to abort if 1
  addOption("randommatingmodel", boolOption, &RandomMatingModel);
  addOption("globalrho", boolOption, &GlobalRho);
  addOption("indadmixhiermodel", boolOption, &IndAdmixHierIndicator);
  addOption("etapriorfile", stringOption, &EtaPriorFilename);
  addOption("globalsumintensitiesprior", dvectorOption, &globalrhoPrior);
  addOption("sumintensitiesprior", dvectorOption, &rhoPrior);
  addOption("etapriormean", doubleOption, &etamean);
  addOption("etapriorvar", doubleOption, &etavar);
  addOption("admixtureprior", dvectorOption, &initalpha[0]);
  addOption("admixtureprior1", dvectorOption, &initalpha[1]);
  addOption("correlatedallelefreqs", boolOption, &correlatedallelefreqs);
  addOption("popadmixproportionsequal", boolOption, &PopAdmixPropsAreEqual);

  //sampler settings
  addOption("rhosamplerparams", fvectorOption, &rhoSamplerParams);
  addOption("popadmixsamplerparams", fvectorOption, &popAdmixSamplerParams);
  // test options
  addOption("ancestryassociationtest", boolOption, &TestForLinkageWithAncestry);
  addOption("affectedsonlytest", boolOption, &TestForAffectedsOnly);
  addOption("admixtureassoctest", boolOption, &TestForAdmixtureAssociation);
  addOption("haplotypeassociationtest", boolOption, &TestForHaplotypeAssociation);
  addOption("stratificationtest", boolOption, &StratificationTestIndicator);
  addOption("allelefreqtest", boolOption, &TestForMisspecifiedAlleleFreqs);
  addOption("allelefreqtest2", boolOption, &TestForMisspecifiedAlleleFreqs2);
  addOption("dispersiontest", boolOption, &TestForDispersion);
  addOption("fstoutput", boolOption, &OutputFST);

  addOption("indadmixmodefile", outputfileOption, &IndAdmixModeFilename);
  addOption("testgenotypesfile", nullOption, 0);
  addOption("locusfortest", intOption, &LocusForTest);
  // Other options
  addOption("chib", boolOption, &chibIndicator);//  Marginal likelihood by Chib algo
  addOption("testoneindiv", boolOption, &TestOneIndivIndicator);//  ML for one individual in a collection 

  /*
    old names for test options, kept for compatibility with R script
  */
  addOption("allelicassociationscorefile", oldOption, 0);
  addOption("affectedsonlyscorefile", oldOption, 0);
  addOption("admixturescorefile", oldOption, 0);
  addOption("ancestryassociationscorefile", oldOption, 0);
  addOption("haplotypeassociationscorefile", oldOption, 0);
  addOption("allelefreqscorefile", oldOption, 0);
  addOption("allelefreqscorefile2", oldOption, 0);
  addOption("stratificationtestfile", oldOption, 0);
  addOption("fstoutputfile", oldOption, 0);
  addOption("hwtestfile", oldOption, 0);
  addOption("residualallelicassocscorefile", oldOption, 0);
}
bool AdmixOptions::SetOptions(){
  if(!Options::SetOptions())
    return false;

  //set indicators
  locusForTestIndicator = (LocusForTest>-1);

  return true;
}

int AdmixOptions::checkOptions(bclib::LogWriter &Log, int NumberOfIndividuals){
  //check base options
  Options::checkOptions(Log, 1);

  bool badOptions = false;//to indicate invalid options. Prog will exit at end of function if true.
  Log.setDisplayMode(Quiet);

  if(useroptions["hapmixmodel"].size() && useroptions["hapmixmodel"] != "0"){
    Log << On << "ERROR: 'hapmixmodel' option is not supported. Consider trying HAPMIXMAP instead\n" << Quiet;
    badOptions = true;
  }

  // **** analysis type  ****
  if(CoxOutcomeVarFilename.length() ){
    Log << "Cox Regression\n";
    ++NumberOfOutcomes;
    if(RegType == None)RegType = Cox;
    else RegType = Multiple;
  }

  if (NumberOfIndividuals==1) {
    IndAdmixHierIndicator = false;
    useroptions["indadmixhiermodel"]="0";
    Log << "One individual analysis";
  } else if (RegType == None) { //no regression
    NumberOfOutcomes = 0;
    if(TestForAffectedsOnly) {
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

  if(OutcomeVarFilename.length() == 0 && OutcomeVarColumns.size()){
    Log << "ERROR: outcomevarcols option specified but no outcomevarfile\n";
    useroptions.erase("outcomevarcols");
    OutcomeVarColumns.clear();
  }
 if(OutcomeVarFilename.length() == 0 && CoxOutcomeVarFilename.length()==0){
    if(CovariatesFilename.length()){
      Log << "ERROR: covariatesfile specified without outcomevarfile\n";
      useroptions.erase("covariatesfile");
      CovariatesFilename.clear();
    }
     if(RegressionOutputFilename.length() > 0){
       Log << "ERROR: regparamfile option is not valid without a regression model\n";
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

      if(OutputAlleleFreq){
      Log << "ERROR: allelefreqoutputfile options is not valid with indadmixhierindicator = 0\n"
	      << "\tThis option will be ignored\n";
         AlleleFreqOutputFilename = "";
         OutputAlleleFreq = false;
      }
      GlobalRho = false;
      useroptions["globalrho"] = "0";
    }
  
  // **** sumintensities ****

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

  //Prior on admixture
  setInitAlpha(Log);
  if(Populations > 1 && IndAdmixHierIndicator){
    Log << "Gamma(1, 1) prior on population admixture Dirichlet parameters.\n";
  }

  // **** model for allele freqs ****

  //fixed allele freqs
  if( (alleleFreqFilename.length()) ||
           (PriorAlleleFreqFilename.length() && fixedallelefreqs ) ){
    Log << "Analysis with fixed allele frequencies.\n";
    if(OutputAlleleFreq ){
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
	  << "Default priors will be set for the allele frequencies"
          << " with " << Populations << " population(s)"
          << "\n";
      if(correlatedallelefreqs) {
	Log << "Analysis with correlated allele frequencies\n";
      }
    }
  
  if( OutputFST && !TestForHaplotypeAssociation ){
    Log << "ERROR: fstoutputfile option is only valid with historicallelefreqfile option\n"
	<< "       this option will be ignored\n";
    OutputFST = false;
    useroptions.erase("fstoutput");
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
			TestForAdmixtureAssociation || TestForHaplotypeAssociation );

  AddFilenamesToUserOptions();

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
void AdmixOptions::setInitAlpha(bclib::LogWriter &Log){
  _admixed.resize(2,(bool)(Populations>1));
  _symmetric = true;
  vector<double> alphatemp(Populations);
  Log.setDisplayMode(Quiet);

  //if no initalpha is specified, alpha for both gametes is initialised to 1.0 for each population  
  if( initalpha[0].size() == 0 && initalpha[1].size() == 0 ){
    fill( alphatemp.begin(), alphatemp.end(), 1.0);//fill alphatemp with 1s
    initalpha[0] = alphatemp; initalpha[1] = alphatemp;//put 2 copies of alphatemp in alpha
    if(!IndAdmixHierIndicator)
      Log << "Dirichlet parameters of prior on admixture: ";
    else 
      Log << "Initial value for population admixture (Dirichlet) parameter vector: ";
    for(int k = 0;k < Populations; ++k){Log << alphatemp[k] << " " ;}
    Log << "\n";

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

void AdmixOptions::PrintUserOptions(const char* filename){
  //set populations value in case it has changed or not specified
  //NB do similar for any option that can be changed outside Options
  std::ostringstream s;
  if (s << getPopulations()) // conversion worked
    {
    useroptions["populations"] = (char *)s.str().c_str();
    }
  useroptions["hapmixmodel"] = "0";
  //Now output Options table to file
  Options::PrintUserOptions(filename);
}

///add names of output files to useroptions so they can be read by R script in args.txt
//TODO: ?? add final table filenames too
void AdmixOptions::AddFilenamesToUserOptions(){
  if(TestForAllelicAssociation )
    useroptions["allelicassociationscorefile"] = ALLELICASSOCTEST_PVALUES;

  if(TestForAffectedsOnly)
    useroptions["affectedsonlyscorefile"] = AFFECTEDSONLYTEST_PVALUES;

  if(TestForAdmixtureAssociation)
    useroptions["admixturescorefile"] = ADMIXTUREASSOCTESTFILE;

  if(TestForLinkageWithAncestry)
    useroptions["ancestryassociationscorefile"] = ANCESTRYASSOCTEST_PVALUES;

  if(TestForHaplotypeAssociation)
    useroptions["haplotypeassociationscorefile"] = HAPLOTYPEASSOCTEST_PVALUES;

  if(TestForMisspecifiedAlleleFreqs)
    useroptions["allelefreqscorefile"] = MISSPECALLELEFREQTEST_1;

  if(TestForMisspecifiedAlleleFreqs2)
    useroptions["allelefreqscorefile2"] = MISSPECALLELEFREQTEST_2;

  if(StratificationTestIndicator)
    useroptions["stratificationtestfile"] = STRAT_TEST_FILE;

  if(TestForDispersion)
    useroptions["dispersiontestfile"] = DISPERSION_TEST_FILE;

  if(OutputFST)
    useroptions["fstoutputfile"] = FST_OUTPUT_FILE;

  if(HWTest)
    useroptions["hwtestfile"] = HARDY_WEINBERG_TEST;

  if(TestForResidualAllelicAssoc)
    useroptions["residualallelicassocscorefile"] = RESIDUAL_LD_TEST_PVALUES;
}
