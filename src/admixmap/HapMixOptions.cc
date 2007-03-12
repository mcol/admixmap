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
//#include "utils/StringConvertor.h"
#include <string.h>
#include <sstream>

using namespace std;

// HapMixOptions::HapMixOptions()
// {
//   Initialise();
// }
HapMixOptions::HapMixOptions(int argc,  char** argv){
  //base class sets default base options
  //this call to SetDefaultValues calls this class' options
  SetDefaultValues();
  ReadUserOptions(argc, argv);
  OptionMap ProgOptions;
  SetOptions(ProgOptions);
}

void HapMixOptions::SetDefaultValues(){
  NumBlockStates = 8;
  FreqDispersionHierModel = true;

  // option names and default option values are stored as strings in a map container 
  // these are default values
  // other specified options will be appended to this array 

  //final values are always output to file, with these default filenames
  FinalFreqPriorFilename = "state-freqpriors.txt";
  FinalLambdaFilename = "state-lambdas.txt";
  AlleleFreqOutputFilename = "state-allelefreqs.txt";

  useroptions["states"] = "4";
  useroptions["hapmixlambdaprior"] = "30, 0.1, 10, 1";
  //useroptions["mixturepropsprior"] = "1,1,1,1,1,1,1,1";
  useroptions["finalfreqpriorfile"] = FinalFreqPriorFilename;
  useroptions["finallambdafile"] = FinalLambdaFilename;
  useroptions["finalallelefreqfile"] = AlleleFreqOutputFilename;
  useroptions["freqdispersionhiermodel"] = "1";
}

HapMixOptions::~HapMixOptions()
{
}

// each option has a function to return its value
bool HapMixOptions::getFixedAlleleFreqs() const
{
  return (fixedallelefreqs);
}

const char *HapMixOptions::getEtaOutputFilename() const
{
  return EtaOutputFilename.c_str();
}


const char *HapMixOptions::getAlleleFreqPriorOutputFilename() const
{
  return AlleleFreqPriorOutputFilename.c_str();
}
bool HapMixOptions::OutputAlleleFreqPrior()const{
  return (bool)(AlleleFreqPriorOutputFilename.size()>0);
}


bool HapMixOptions::getMHTest()const{
  return (bool)(MHTestFilename.size());
}
const char* HapMixOptions::getMHTestFilename()const{
  return MHTestFilename.c_str();
}

bool HapMixOptions::getHapMixModelIndicator() const{
  return true;
}
bool HapMixOptions::isFreqDispersionHierModel()const{
  return FreqDispersionHierModel;
}
int HapMixOptions::getNumberOfBlockStates() const
{
  return NumBlockStates;
}
int HapMixOptions::getPopulations() const
{
  return NumBlockStates;
}

void HapMixOptions::setPopulations(int num)
{
  NumBlockStates = num;
}
const vector<float>& HapMixOptions::getLambdaSamplerParams()const{
  return rhoSamplerParams;
}

const std::vector<double> &HapMixOptions::getHapMixLambdaPrior()const{
  return hapmixlambdaprior;
}

const std::vector<double>& HapMixOptions::getMixturePropsPrior()const{
  return MixturePropsPrior;
}
const char* HapMixOptions::getInitialAlleleFreqFilename()const{
  return InitialAlleleFreqFilename.c_str();
}
const char* HapMixOptions::getInitialHapMixLambdaFilename()const{
    return InitialHapMixLambdaFilename.c_str();
}
const char* HapMixOptions::getInitialFreqPriorFilename()const{
  return InitialFreqPriorFile.c_str();
}

const std::vector<double> & HapMixOptions::getAlleleFreqPriorParams()const{
  return allelefreqprior;
}

const char* HapMixOptions::getHapMixLambdaOutputFilename()const{
    return HapMixLambdaOutputFilename.c_str();
}
const char* HapMixOptions::getCCGenotypesFilename()const{
  return CCGenotypesFilename.c_str();
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

const char* HapMixOptions::getFinalFreqPriorFilename()const{
  return FinalFreqPriorFilename.c_str();
}
const char* HapMixOptions::getFinalLambdaFilename()const{
  return FinalLambdaFilename.c_str();
}

unsigned HapMixOptions::GetNumMaskedIndividuals()const{
  return MaskedIndividuals.size();
}
unsigned HapMixOptions::GetNumMaskedLoci()const{
  return MaskedLoci.size();
}

void HapMixOptions::SetOptions(OptionMap& ProgOptions)
{
  //set up Option map
  //first set options for this class

  ProgOptions["states"] = OptionPair(&NumBlockStates, "int");
  ProgOptions["ccgenotypesfile"] = OptionPair(&CCGenotypesFilename, "string");
  //standard output files (optional)

  //file to write sampled values of dispersion parameter (ADMIXMAP) or mean and variance of dispersion params (HAPMIXMAP)
  ProgOptions["dispparamfile"] = OptionPair(&EtaOutputFilename, "outputfile");// C

  //final values
  ProgOptions["finalfreqpriorfile"] = OptionPair(&FinalFreqPriorFilename, "outputfile");
  ProgOptions["finallambdafile"] = OptionPair(&FinalLambdaFilename, "outputfile");
  ProgOptions["finalallelefreqfile"] = OptionPair(&AlleleFreqOutputFilename, "outputfile");// synonym for allelefreqoutputfile

  //posterior means
  ProgOptions["allelefreqprioroutputfile"] = OptionPair(&AlleleFreqPriorOutputFilename, "outputfile");
  ProgOptions["hapmixlambdaoutputfile"] = OptionPair(&HapMixLambdaOutputFilename, "outputfile");// remove hapmix from name

  //prior and model specification
  ProgOptions["hapmixmodel"] = OptionPair(0, "null");//does nothing
  ProgOptions["hapmixlambdaprior"] = OptionPair(&hapmixlambdaprior, "dvector");
  ProgOptions["mixturepropsprior"] = OptionPair(&MixturePropsPrior, "dvector");
  ProgOptions["allelefreqprior"] = OptionPair(&allelefreqprior, "dvector");
  ProgOptions["freqdispersionhiermodel"] = OptionPair(&FreqDispersionHierModel, "bool");

  //initial values
  ProgOptions["initiallambdafile"] = OptionPair(&InitialHapMixLambdaFilename, "string");
  ProgOptions["initialfreqpriorfile"] = OptionPair(&InitialFreqPriorFile, "string");
  ProgOptions["initialallelefreqfile"] = OptionPair(&InitialAlleleFreqFilename, "string");// synonym for allelefreqfile

  //sampler settings
  ProgOptions["lambdasamplerparams"] = OptionPair(&rhoSamplerParams, "fvector");//A

  // test options
  ProgOptions["mhscoretestfile"] = OptionPair(&MHTestFilename, "outputfile");
  ProgOptions["maskedindivs"] = OptionPair(&MaskedIndividuals, "range");
  ProgOptions["maskedloci"] = OptionPair(&MaskedLoci, "range");

  //now set base options and finish parsing
  Options::SetOptions(ProgOptions);

}

int HapMixOptions::checkOptions(LogWriter &Log, int ){
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

  if (RegType == None) { //no regression
    NumberOfOutcomes = 0;
    {
      Log << "Cross sectional analysis, no outcome";
    }
  } else if (RegType == Linear) {
    Log << "Cross sectional analysis, continuous outcome";
  } else if (RegType == Logistic) {
    Log << "Case-control or cross sectional analysis, binary outcome";
  } else{ //if (RegType == Both) {
    Log << "Cross sectional analysis, multiple outcome";
  }
  Log << "\n";

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
  
  // **** sumintensities ****
  Log << "Haplotype mixture model with " << NumBlockStates << " block states\n";
  
  //check mixture props prior has length K
  if(MixturePropsPrior.size()>0 & MixturePropsPrior.size()!= (unsigned)NumBlockStates){
    Log << "Error: 'mixturepropsprior' has incorrect length\n";
    badOptions = true;
  }
  if(!MixturePropsPrior.size()){//user has not specified a prior
    //MixtureProps.assign(1, NumBlockStates);//redundant since set in PopHapMix
    useroptions["mixturepropsprior"] = "1";
    for(int k = 1; k < NumBlockStates; ++k)
      useroptions["mixturepropsprior"].append(",1");
  }

  // **** model for allele freqs ****

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

  //fixed allele freqs
  if( PriorAlleleFreqFilename.length() && fixedallelefreqs  ){
    Log << "Analysis with fixed allele frequencies.\n";
  }
  //prior allele freqs
  else if( PriorAlleleFreqFilename.length() && !fixedallelefreqs ){
    Log << "Analysis with prior allele frequencies.\n";
  }
  //default priors ('states' option)
  else if(NumBlockStates > 0 )
    {
      Log << "Default priors will be set for the allele frequencies\n";
      }
  
  ScoreTestIndicator = (TestForAllelicAssociation || TestForResidualAllelicAssoc);

  if(thermoIndicator) {
    // for thermo integration, NumAnnealedRuns is set to default value of 100 
    // if not specified as an option
    //if(NumAnnealedRuns==0) NumAnnealedRuns = 100;
    Log << "\nUsing thermodynamic integration to calculate marginal likelihood ";
  }

  if(badOptions) return 1;
  else return 0;
}

void HapMixOptions::PrintOptions(){
  //set states value in case it has changed or not specified
  std::ostringstream s;
  if (s << getNumberOfBlockStates()) // conversion worked
    {
    useroptions["states"] = (char *)s.str().c_str();
    }
  //Now output Options table to args.txt
  Options::PrintOptions();
}
