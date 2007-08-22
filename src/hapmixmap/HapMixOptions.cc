/** 
 *   HAPMIXMAP
 *   HapMixOptions.cc 
 *   Class to hold HAPMIXMAP options
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "HapMixOptions.h"
#include "HapMixFilenames.h"
#include <string.h>
#include <sstream>

using namespace std;
using namespace bclib;

HapMixOptions::HapMixOptions(int argc,  char** argv){
  SetDefaultValues();
  DefineOptions();
  if(!ReadUserOptions(argc, argv, "-f"))
    exit(1);
}

void HapMixOptions::SetDefaultValues(){
  NumBlockStates = 8;
  FreqPrecisionHierModel = false;
  FixedMixtureProps = true;
  FixedMixturePropsPrecision = true;
  MixturePropsPrecision = 0.0;
  MHTestIndicator = false;
  NumStarts = 0;

  // option names and default option values are stored as strings in a map container 
  // these are default values
  // other specified options will be appended to this array 

  useroptions["states"] = "8";
  useroptions["arrivalrateprior"] = "12, 0.5, 2";
  useroptions["residualadhiermodel"] = "0";
  useroptions["fixedmixtureprops"] = "1";
  useroptions["fixedmixturepropsprecision"] = "1";
}

HapMixOptions::~HapMixOptions(){
}

// each option has a function to return its value

bool HapMixOptions::outputParams()const{
  return ( ParameterFilename.length() || RegressionOutputFilename.length() || FreqPrecisionOutputFilename.length());
}

bool HapMixOptions::getFixedAlleleFreqs() const{
  return fixedallelefreqs;
}

bool HapMixOptions::getFixedMixtureProps()const{
  return FixedMixtureProps;
}

bool HapMixOptions::getFixedMixturePropsPrecision()const{
  return FixedMixturePropsPrecision;
}

const char *HapMixOptions::getFreqPrecisionOutputFilename() const
{
  return FreqPrecisionOutputFilename.c_str();
}


const char *HapMixOptions::getAlleleFreqPriorOutputFilename() const{
  return AlleleFreqPriorOutputFilename.c_str();
}
bool HapMixOptions::OutputAlleleFreqPrior()const{
  return (bool)(AlleleFreqPriorOutputFilename.size()>0);
}


bool HapMixOptions::getMHTest()const{
  return MHTestIndicator;
}
bool HapMixOptions::getHapMixModelIndicator() const{
  return true;
}
bool HapMixOptions::isFreqPrecisionHierModel()const{
  return FreqPrecisionHierModel;
}
int HapMixOptions::getNumberOfBlockStates() const{
  return NumBlockStates;
}
int HapMixOptions::getPopulations() const{
  return NumBlockStates;
}

void HapMixOptions::setPopulations(int num){
  NumBlockStates = num;
}

unsigned HapMixOptions::GetNumStarts()const{
  return NumStarts;
}
const vector<float>& HapMixOptions::getLambdaSamplerParams()const{
  return rhoSamplerParams;
}

const std::vector<double> &HapMixOptions::getLambdaPrior()const{
  return lambdaprior;
}

const std::vector<double>& HapMixOptions::getMixturePropsPrecisionPrior()const{
  return MixturePropsPrecisionPrior;
}
float HapMixOptions::getMixturePropsPrecision()const{
  return MixturePropsPrecision;
}

string HapMixOptions::getInitialValuePath(unsigned startindex, const string& filename)const{
  stringstream ss;
  ss << InitialValueDir << "/start" << startindex+1  << "/" << filename;
  return ss.str();
}
std::string HapMixOptions::getInitialAlleleFreqFilename(unsigned startindex)const{
  if(NumStarts < 2)//no initial values (filename is empty) or a single start
    return InitialAlleleFreqFilename;
  else {
    return getInitialValuePath(startindex, InitialAlleleFreqFilename);
  }
}
std::string HapMixOptions::getInitialArrivalRateFilename(unsigned startindex)const{
  if(NumStarts < 2)//no initial values (filename is empty) or a single start
    return InitialArrivalRateFilename;
  else {
    return getInitialValuePath(startindex, InitialArrivalRateFilename);
  }
}
std::string HapMixOptions::getInitialMixturePropsFilename(unsigned startindex)const{
  if(NumStarts < 2)//no initial values (filename is empty) or a single start
    return InitialMixturePropsFilename;
  else {
    return getInitialValuePath(startindex, InitialMixturePropsFilename);
  }
}
std::string HapMixOptions::getInitialFreqPriorFilename(unsigned startindex)const{
  if(NumStarts < 2)//no initial values (filename is empty) or a single start
    return InitialFreqPriorFilename;
  else {
    return getInitialValuePath(startindex, InitialFreqPriorFilename);
  }
}

const std::vector<double> & HapMixOptions::getAlleleFreqPriorParams()const{
  return allelefreqprecisionprior;
}

const char* HapMixOptions::getArrivalRateOutputFilename()const{
    return ArrivalRateOutputFilename.c_str();
}
const char* HapMixOptions::getTestGenotypesFilename()const{
  return TestGenotypesFilename.c_str();
}

const vector<unsigned>& HapMixOptions::getMaskedIndividuals()const{
  return MaskedIndividuals;
}
const vector<unsigned>& HapMixOptions::getMaskedLoci()const{
  return MaskedLoci;
}
bool HapMixOptions::OutputCGProbs()const{
  return (MaskedLoci.size() && MaskedIndividuals.size()); 
}
string HapMixOptions::getFinalAlleleFreqFilename()const{
  return (FinalValueDir + "/" + ALLELEFREQSTATEFILE);
}
string HapMixOptions::getFinalFreqPriorFilename()const{
  return (FinalValueDir + "/" + FREQPRIORSTATEFILE);
}
string HapMixOptions::getFinalLambdaFilename()const{
  return (FinalValueDir + "/" + ARRIVALRATESSTATEFILE);
}
string HapMixOptions::getFinalMixturePropsFilename()const{
  return (FinalValueDir + "/" + MIXTUREPROPSSTATEFILE);
}
const string& HapMixOptions::getFinalValueDir()const{
  return FinalValueDir;
}
unsigned HapMixOptions::GetNumMaskedIndividuals()const{
  return MaskedIndividuals.size();
}
unsigned HapMixOptions::GetNumMaskedLoci()const{
  return MaskedLoci.size();
}

void HapMixOptions::DefineOptions()
{
  //set up Option map

  /*
    data and initial values
  */
  addOption("ccgenotypesfile", stringOption, &TestGenotypesFilename);//obsolescent
  addOption("testgenotypesfile", stringOption, &TestGenotypesFilename);
  addOption("initialvaluedir", stringOption, &InitialValueDir);
  addOption("numberofstarts", intOption, &NumStarts);

  /*
     model specification
  */
  addOption("hapmixmodel", nullOption, 0);//does nothing
  //number of block states
  addOption("states", intOption, &NumBlockStates);
  //toggle hierarchical model on allele freq precision
  addOption("freqprecisionhiermodel", boolOption, &FreqPrecisionHierModel);
  addOption("residualadhiermodel", boolOption, &FreqPrecisionHierModel);
  //toggle sampling of mixture props
  addOption("fixedmixtureprops", boolOption, &FixedMixtureProps);
  //toggle sampling of mixture props precision, valid only if fixedmixtureprops=0
  addOption("fixedmixturepropsprecision", boolOption, &FixedMixturePropsPrecision);
  //set mixture props precision, initial value only if fixedmixturepropsprecision=0
  addOption("mixturepropsprecision", floatOption, &MixturePropsPrecision);

  /*
    prior specification
  */
  //vector of length 3 or 4, 2 Gamma priors on Gamma parameters
  addOption("arrivalrateprior", dvectorOption, &lambdaprior);
  //vector of length 2, Gamma prior on Dirichlet precision
  addOption("mixturepropsprecisionprior", dvectorOption, &MixturePropsPrecisionPrior);
  //vector of length 2 or 3, Gamma/Gamma-Gamma prior on allele freq prior precision
  addOption("allelefreqprecisionprior", dvectorOption, &allelefreqprecisionprior);
  addOption("residualadprior", dvectorOption, &allelefreqprecisionprior);

  //sampler settings
  addOption("arrivalratesamplerparams", fvectorOption, &rhoSamplerParams);

  /*
    Output files
  */
  //file to write mean and variance of frequency precision params
  addOption("freqprecisionfile", outputfileOption, &FreqPrecisionOutputFilename);
  addOption("residualadfile", outputfileOption, &FreqPrecisionOutputFilename);

  //final values
  addOption("finalvaluedir", stringOption, &FinalValueDir);

  //posterior means
  addOption("allelefreqprecisionposteriormeanfile", outputfileOption, &AlleleFreqPriorOutputFilename);
  addOption("residualadposteriormeanfile", outputfileOption, &AlleleFreqPriorOutputFilename);
  addOption("arrivalrateposteriormeanfile", outputfileOption, &ArrivalRateOutputFilename);


  /*
    test options
  */
  //Mantel-Haentszel test
  addOption("mhtest", boolOption, &MHTestIndicator);
  //for output of posterior predictive genotype probs
  addOption("maskedindivs", rangeOption, &MaskedIndividuals);
  addOption("maskedloci", rangeOption, &MaskedLoci);

  /*
    old names for test options, kept for compatibility with R script
  */
  addOption("mhscoretestfile", oldOption, 0);
  addOption("allelicassociationscorefile", oldOption, 0);
  addOption("residualallelicassocscorefile", oldOption, 0);
  addOption("hwtestfile", oldOption, 0);

}

int HapMixOptions::checkOptions(bclib::LogWriter &Log, int ){
  //check base options
  Options::checkOptions(Log, 1);

  bool badOptions = false;//to indicate invalid options. Prog will exit at end of function if true.
  Log.setDisplayMode(Quiet);

  // **** analysis type  ****
  if(CoxOutcomeVarFilename.length() ){
    Log << "Cox Regression\n";
    ++NumberOfOutcomes;
    if(RegType == None)RegType = Cox;
    else RegType = Multiple;
  }

  if (RegType == None) { //no regression
    NumberOfOutcomes = 0;
    Log << "Cross sectional analysis, no outcome";
  } else if (RegType == Linear) {
    Log << "Cross sectional analysis, continuous outcome";
  } else if (RegType == Logistic) {
    Log << "Case-control or cross sectional analysis, binary outcome";
  } else{ //if (RegType == Both) {
    Log << "Cross sectional analysis, multiple outcome";
  }
  Log << "\n";

  if(TestGenotypesFilename.size() && !OutcomeVarFilename.size()){
    Log << "WARNING: observed genotypes supplied but no outcome variable!\n";
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
  
  Log << "Haplotype mixture model with " << NumBlockStates << " block states\n";
  if(NumBlockStates < 2){
    Log << "ERROR: states must be at least 2\n";
    badOptions = true;
  }

  if(lambdaprior.size()<3 || lambdaprior.size()>4){
    Log << "ERROR: 'arrivalrateprior' should have 3 or 4 elements\n";
    badOptions = true;
  }

  if(FixedMixtureProps){
    if(!FixedMixturePropsPrecision){
      useroptions.erase("mixturepropsprecisionprior");
      useroptions["fixedmixturepropsprecision"] = "1"; 
      FixedMixturePropsPrecision = true;
      Log << "WARNING: option 'fixedmixtureprops=1' requires 'fixedmixturepropsprecision=1'\n";
    }
    //    useroptions.erase["finalmixturepropsfile"];
  }else{//we are sampling mixture props
    if(!FixedMixturePropsPrecision){//we are also sampling mixture props precision
      //check mixture props prior has length 2 and both elements are positive
      if( MixturePropsPrecisionPrior.size() && //if user has specified a prior
	  (MixturePropsPrecisionPrior.size()!=2 || MixturePropsPrecisionPrior[0] <= 0.0 || MixturePropsPrecisionPrior[1] <= 0.0) ){
	Log << "ERROR: 'mixturepropsprecisionprior' should be 2 positive numbers\n";
	badOptions = true;
      }

      if(!MixturePropsPrecisionPrior.size()){//user has not specified a prior
	/*//redundant since set in PopHapMix
	  MixtureProps.push_back(NumBlockStates);
	  MixtureProps.push_back(1);
	*/
	std::stringstream ss;
	ss << NumBlockStates << ", 1";
	useroptions["mixturepropsprecisionprior"] = ss.str();
      }
    }
  }// end random mixture props block
  if(FixedMixturePropsPrecision  && MixturePropsPrecisionPrior.size()){
      Log << "WARNING: option 'mixturepropsprecisionprior' is not valid with 'fixedmixturepropsprecision'\n";
      useroptions.erase("mixturepropsprecisionprior");
  }

  // **** model for allele freqs ****

  if(useroptions.find("allelefreqprecisionprior") != useroptions.end() ||
     useroptions.find("residualadprior") != useroptions.end()){
    
    if(FreqPrecisionHierModel && (allelefreqprecisionprior.size() !=3)) {
      Log << "ERROR: 'residualadprior' must have length 3\n";
      badOptions = true;
    }
    else if(allelefreqprecisionprior.size()< 2){
      Log << "ERROR: 'residualadprior' must have length 2\n";
      badOptions = true;
    }
  }

  //freqprecisionfile not allowed if frequencies are fixed
  if(fixedallelefreqs){
      //TODO: ?? warnings
    if( FreqPrecisionOutputFilename.size()){
      FreqPrecisionOutputFilename.clear();
      useroptions.erase("freqprecisionfile");
      useroptions.erase("residualadfile");
    }
    if(AlleleFreqPriorOutputFilename.size()){
      AlleleFreqPriorOutputFilename.clear();
      useroptions.erase("allelefreqprecisionposteriormeanfile");
      useroptions.erase("residualadposteriormeanfile");
    }
  }


  //fixed allele freqs
  if( PriorAlleleFreqFilename.length() && fixedallelefreqs  ){
    Log << "Analysis with fixed allele frequencies.\n";
  }
  //prior allele freqs
  else if( PriorAlleleFreqFilename.length() && !fixedallelefreqs ){
    Log << "Analysis with prior allele frequencies.\n";
  }
  //default priors ('states' option)
  else if(NumBlockStates > 0 )  {
    Log << "Default priors will be set for the allele frequencies\n";
  }
  
  if(thermoIndicator) {
    // for thermo integration, NumAnnealedRuns is set to default value of 100 
    // if not specified as an option
    //if(NumAnnealedRuns==0) NumAnnealedRuns = 100;
    Log << "\nUsing thermodynamic integration to calculate marginal likelihood ";
  }

  //check initial value dir is specified with multiple starts
  if(NumStarts > 0){
    if(!InitialValueDir.size()){
      Log << "ERROR: multiple starts specified but initial value directory not set\n";
      badOptions = true;
    }
  }

  if(NumStarts == 0)
    NumStarts = 1;
  //set initial value filenames
  if(InitialValueDir.size()){
    if(NumStarts < 2){
      InitialArrivalRateFilename = InitialValueDir + "/" + ARRIVALRATESSTATEFILE;
      InitialMixturePropsFilename = InitialValueDir + "/" + MIXTUREPROPSSTATEFILE;
      InitialAlleleFreqFilename = InitialValueDir + "/" + ALLELEFREQSTATEFILE;
      InitialFreqPriorFilename = InitialValueDir + "/" + FREQPRIORSTATEFILE;
    }
    else{
      InitialArrivalRateFilename = ARRIVALRATESSTATEFILE;
      InitialMixturePropsFilename = MIXTUREPROPSSTATEFILE;
      InitialAlleleFreqFilename = ALLELEFREQSTATEFILE;
      InitialFreqPriorFilename = FREQPRIORSTATEFILE;
    }
  }

  //default final value dir to results dir if not specified
  if(!FinalValueDir.size()){
    FinalValueDir = ResultsDir;
    useroptions["finalvaluedir"] = ResultsDir;
  }
//check maskedloci and maskedindivs options
  if(MaskedLoci.size()){
      if(!MaskedIndividuals.size()){
          Log << On << "Error: maskedloci option specified without maskedindivs\n";
          MaskedLoci.clear();
      }
      else if(MaskedLoci[0] <= 0){
          Log << On << "Warning: " << MaskedLoci[0] 
              << " is not a valid value for a masked locus and will be ignored\n";
          MaskedLoci.erase(MaskedLoci.begin());
      }
  }
  if(MaskedIndividuals.size()){
      if(!MaskedLoci.size()){
          Log << On << "Error: maskedindivs option specified without maskedloci\n";
          MaskedIndividuals.clear();
      }
      else if(MaskedIndividuals[0] <= 0){
          Log << On << "Warning: " << MaskedIndividuals[0] 
              << " is not a valid value for a masked individual and will be ignored\n";
          MaskedIndividuals.erase(MaskedIndividuals.begin());
      }
  }

  AddFilenamesToUserOptions();

  if(badOptions) return 1;
  else return 0;
}

void HapMixOptions::PrintUserOptions(const char* filename){
  //set states value in case it has changed or not specified
  std::ostringstream s;
  if (s << getNumberOfBlockStates()) // conversion worked
    {
    useroptions["states"] = (char *)s.str().c_str();
    }
  //Now output Options table to file
  Options::PrintUserOptions(filename);
}

///add names of output files to useroptions so they can be read by R script in args.txt
//TODO: ?? add final table filenames too
void HapMixOptions::AddFilenamesToUserOptions(){
  if(TestForAllelicAssociation )
    useroptions["allelicassociationscorefile"] = ALLELICASSOCTEST_PVALUES;

  if(HWTest)
    useroptions["hwtestfile"] = HARDY_WEINBERG_TEST;

  if(TestForResidualAllelicAssoc)
    useroptions["residualallelicassocscorefile"] = RESIDUAL_LD_TEST_PVALUES;

  if(MHTestIndicator)
    useroptions["mhscoretestfile"] = MH_TEST_PVALUES;
}
