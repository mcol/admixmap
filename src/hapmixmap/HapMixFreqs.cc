/** 
 *   HAPMIXMAP
 *   HapMixFreqs.cc 
 *   Class to sample prior parameters of frequencies in a hapmixmodel
 *   Copyright (c) 2006, 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HapMixFreqs.h"
#include "Genome.h"
#include "HapMixOptions.h"
#include "bclib/misc.h"
#include "bclib/Exceptions.h"
#include <iomanip>

HapMixFreqs::HapMixFreqs(){
  DirichletParams = 0;
  Eta = 0;
  EtaSampler = 0;
  SumLambda = 0;
  SumEta = 0;
  NumEtaUpdates = 0;
  accumulateEta = false;
  etaHierModel = true;
  DiploidGenotypeProbs.array = 0;
  FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
}

HapMixFreqs::~HapMixFreqs(){
  delete[] DirichletParams;
  delete[] Eta;
  delete[] EtaSampler;
  delete[] SumEta;
  if( allelefreqprioroutput.is_open()) allelefreqprioroutput.close();
  //TODO: delete DiploidGenotypeProbs
}

void HapMixFreqs::Initialise(HapMixOptions* const options, InputData* const data, Genome *pLoci, bclib::LogWriter &Log ){
  Loci = pLoci;
  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  //set model indicators
  RandomAlleleFreqs = !options->getFixedAlleleFreqs();
  
  LoadAlleleFreqs(options, data, Log);
  
  InitialisePrior(Populations, NumberOfCompositeLoci, options, Log );  
  
  //open freqpriorfile
  if(IsRandom() ){
    OpenOutputFile(options->getFreqPrecisionOutputFilename());
  }
  
  AllocateAlleleCountArrays(options->getPopulations());
}

void HapMixFreqs::setSampler(bool thermo, bool AllHaploid, bool /*DefaultPriors*/){
//   // set which sampler will be used for allele freqs
//   // current version uses conjugate sampler if annealing without thermo integration
//   if( (options->getThermoIndicator() ) ||
//       //using default allele freqs or CAF model
//       (  !strlen(options->getPriorAlleleFreqFilename()) && !strlen(options->getInitialAlleleFreqFilename()) ) ) {
//     FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
//   } else {
//     FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
//   }

  if(!thermo && (AllHaploid /*|| !DefaultPriors*/))
    FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
  else//thermo, some diploid data, default priors
    FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
  
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    if(RandomAlleleFreqs){
      if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER){
	//set up samplers for allelefreqs
	FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
						    &(DirichletParams[i]), true));
      }
    }
    //set AlleleProbs pointers in CompositeLocus objects to point to Freqs
    //allocate HapPairProbs and calculate them using AlleleProbs
    (*Loci)(i)->InitialiseHapPairProbs(Freqs[i], AllHaploid);
    
  }//end comp locus loop

}

void HapMixFreqs::AllocateDiploidGenotypeProbs(){
  DiploidGenotypeProbs.array = new double*[NumberOfCompositeLoci];
  for(int i = 0; i < NumberOfCompositeLoci; ++i)
     DiploidGenotypeProbs.array[i] = new double[3*Populations*Populations];
}

void HapMixFreqs::SetDiploidGenotypeProbs(){
if(DiploidGenotypeProbs.array)
  for(int i = 0; i < NumberOfCompositeLoci; ++i){
    for(int k0 = 0; k0 < Populations; ++k0)
      for(int k1 = 0; k1 < Populations; ++k1){

      //1,1   =  phi_k0,1 * phi_k1,1
       DiploidGenotypeProbs[i][k0*Populations*3 + k1*3   ] = Freqs[i][k0*2] * Freqs[i][k1*2];
       //2,2  =  phi_k0,2 * phi_k1,2
       DiploidGenotypeProbs[i][k0*Populations*3 + k1*3 +1] = Freqs[i][k0*2 + 1] * Freqs[i][k1*2+1];
       //1,2  =  phi_k0,1 * phi_k1,2 + phi_k0,2 *  phi_k1,1
       DiploidGenotypeProbs[i][k0*Populations*3 + k1*3 +2] = Freqs[i][k0*2] * Freqs[i][k1*2 + 1] + Freqs[i][k0*2 + 1] * Freqs[i][k1*2];
     }
  }
}

void HapMixFreqs::InitialisePrior(unsigned Populations, unsigned L, const HapMixOptions* const  options, bclib::LogWriter& Log ){
  NumberOfCompositeLoci = L;
  //allocate prior arrays
  DirichletParams = new double[NumberOfCompositeLoci];//1D array of prior params
  Eta = new double[NumberOfCompositeLoci];
  if(options->OutputAlleleFreqPrior()){
    accumulateEta = true;
    SumEta = new double[NumberOfCompositeLoci];
  }
  
  EtaSampler = new bclib::StepSizeTuner[NumberOfCompositeLoci];

  //set parameters of prior on frequency Dirichlet prior params
  const std::vector<double> &params = options->getAlleleFreqPriorParams();
  etaHierModel = options->isFreqPrecisionHierModel();

  if(etaHierModel){//hierarchical model on precision
    if(params.size()==3) {//user-specified prior
      EtaShape = params[0];
      EtaRatePriorShape = params[1];
      EtaRatePriorRate = params[2];
    }
    else{//use defaults
      EtaShape = 1.0;//0.25;
      EtaRatePriorShape = 40.0;//10.0;
      EtaRatePriorRate = 10.0;
    }
    EtaRate = EtaRatePriorShape / EtaRatePriorRate;
  }
  else{//no hierarchical model on precision
    //the values of these 2 should be irrelevant
    EtaRatePriorRate = 10.0;    
    EtaRatePriorShape = 10.0;

    if(params.size()>=2) {//user-specified prior
      EtaShape = params[0];
      EtaRate = params[1];
    }
    else{//use defaults
      EtaShape = 1.0;
      EtaRate = 4.0;
    }
  }

  HapMixMuArgs.K = Populations;
  //initialise sampler for freq precision Dirichlet prior proportions, bounded by 0 and 1
  MuSampler.Initialise(true, true, 1.0, 0.0, fmu_hapmix, dfmu_hapmix );

  const string initialvaluefilename = options->getInitialFreqPriorFilename();
  if(initialvaluefilename.size()){
    //TODO??: check prior is consistent with initial values
    ReadInitialPriorParamsFromFile(initialvaluefilename.c_str(), Log);
  }
  else{
    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      //set precision to prior mean
      Eta[i] = (EtaShape / EtaRate);
      //set proportions to 0.5
      DirichletParams[i] =  Eta[i] * 0.5;
    }
  }
  //initialise sampler and set cumulative sums to 0
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    EtaSampler[i].SetParameters(0.1, 0.00001, 100.0, 0.26);
    if(accumulateEta)SumEta[i] = 0.0;
  }

}

void HapMixFreqs::ReadInitialPriorParamsFromFile(const char* filename, bclib::LogWriter& Log){
  std::ifstream initialvaluefile(filename);
  if(initialvaluefile.is_open()){
    Log << bclib::Quiet << "Reading initial values of allele freq prior from " << filename << "\n";

    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      //read eta (precision) and alpha (mean)
      if(!( initialvaluefile >> Eta[i] >> DirichletParams[i]))
	{//reached end-of-file too early
	  throw string("ERROR: Too few entries in initialfreqpriorfile\n");
	}
      //check values read are ok
      if(Eta[i] < 0.0)
	throw bclib::DataOutOfRangeException("residual allelic diversity parameter", ">0", "initialfreqpriorfile");
      if( DirichletParams[i] < 0.0 || DirichletParams[i] > 1.0)
	throw bclib::DataOutOfRangeException("allele frequency mean", "between 0 and 1", "initialfreqpriorfile");
      DirichletParams[i] *= Eta[i]; //multiply mean by precision to get Dir params
    }

   //see if anything more than whitespace left in file
    string test;
    initialvaluefile >> test;
    if(test.find_first_not_of(" \t\n\r") != string::npos){
      throw string("ERROR: too many entries in initialfreqpriorfile\n");
    }
    initialvaluefile.close();
    
  }else{
    string err("ERROR: cannot open initialfreqpriorfile: ");
    err.append(filename);
    throw err;
  }

}
void HapMixFreqs::PrintPrior(bclib::LogWriter& Log)const{
  Log << "Dirichlet prior on allele frequencies. \n"
      << "Uniform prior on Dirichlet proportions.\n"
      << "Gamma prior on Dirichlet precision (Residual Allelic Diversity) with shape " << EtaShape;
  if(etaHierModel)
    Log << " and Gamma("<< EtaRatePriorShape << ", " << EtaRatePriorRate << ") prior on rate.\n\n";
  else
    Log << " and rate " << EtaRate << ".\n\n";

}

void HapMixFreqs::LoadAlleleFreqs(HapMixOptions* const options, InputData* const data_, bclib::LogWriter &Log)
{
  int newrow;
  int row = 0;
  const Matrix_s* temporary = 0;

  //allocate frequency arrays
  Freqs.array = new double*[NumberOfCompositeLoci];
  AlleleFreqsMAP.array = Freqs.array;

  //set static members of CompositeLocus
  CompositeLocus::SetRandomAlleleFreqs(RandomAlleleFreqs);
  CompositeLocus::SetNumberOfPopulations(Populations);

  //read initial values from file
  const string initfilename = options->getInitialAlleleFreqFilename();
  if(initfilename.size() ){
    LoadInitialAlleleFreqs(initfilename.c_str(), Log);
  }
  else{
    int offset = 0;
    bool file = false;//flag to indicate if priors have been supplied in a file
    
    //priorallelefreqfile option
    if( strlen( options->getPriorAlleleFreqFilename() ) ){
      offset = 0;
      file = true;
      temporary = &(data_->getPriorAlleleFreqData());
    }
    
    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      Freqs.array[i] = new double[Loci->GetNumberOfStates(i)* Populations];
      
      if(file){//read allele freqs from file
        newrow = row + (*Loci)(i)->GetNumberOfStates() - offset;
        AlleleFreqs::LoadAlleleFreqs( *temporary, i, row+1, false);//row+1 is the first row for this locus (+1 for the header)
        row = newrow;
      }
      else {  //set default Allele Freqs
        SetDefaultAlleleFreqs(i);
      }
    }
  }//end else
}

///opens file for output of precision mean and var and prior rate
void HapMixFreqs::OpenOutputFile(const char* filename){
  if(strlen(filename)){
    allelefreqprioroutput.open(filename);
    allelefreqprioroutput << "RAD.Mean\tRAD.Var";
    if(etaHierModel)allelefreqprioroutput << "\tRAD.Prior.Rate";
    allelefreqprioroutput << std::endl;
  }
}

/// samples allele frequencies and prior parameters.
void HapMixFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness, bool AllHaploid){
  
  // Sample allele frequencies conditional on Dirichlet priors 
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized

  if(etaHierModel){//sample rate parameter of gamma prior on eta
    double scalarSumEta = 0.0;
    for( int i = 0; i < NumberOfCompositeLoci; i++ )scalarSumEta += Eta[i];
    SampleEtaRate(afterBurnIn, scalarSumEta);
  }


  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    const unsigned NumberOfStates = Loci->GetNumberOfStates(i);
    //accumulate summary stats for update of priors
    double sumlogfreqs1 = 0.0,sumlogfreqs2 = 0.0; 
    for(int k = 0; k < Populations; ++k){
      sumlogfreqs1 += bclib::eh_log(Freqs[i][k*NumberOfStates]);//allele 1
      sumlogfreqs2 += bclib::eh_log(Freqs[i][k*NumberOfStates + 1]);//allele 2
    }
    SamplePriorPrecision(i, Populations, sumlogfreqs1, sumlogfreqs2);
    if(accumulateEta){
      SumEta[i] += Eta[i];
      ++NumEtaUpdates;
    }
    SamplePriorProportions(i, sumlogfreqs1, sumlogfreqs2 );
    
    //Sample prior parameters
    if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER && (NumberOfStates > 2 ||
						  accumulate(hetCounts[i], hetCounts[i]+Populations*Populations, 0, std::plus<int>()) > 0 )) {
      if(NumberOfStates==2) //shortcut for SNPs
	FreqSampler[i]->SampleSNPFreqs(Freqs[i], AlleleCounts[i], hetCounts[i], i, Populations, 
				       coolness);
      else FreqSampler[i]->SampleAlleleFreqs(Freqs[i], IC, i, NumberOfStates, Populations, 
					     coolness);
    }
    else //if (FREQSAMPLER==FREQ_CONJUGATE_SAMPLER)
      SampleAlleleFreqs(i, coolness);
    
    if(afterBurnIn)
      (*Loci)(i)->AccumulateAlleleProbs();

    //no need to update alleleprobs, they are the same as Freqs
    //set HapPair probs using updated alleleprobs
    if(!AllHaploid)//skip if all haploid data
      (*Loci)(i)->SetHapPairProbs();

  }
  
}

/** samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
 and stores result in array Freqs 
*/
void HapMixFreqs::SampleAlleleFreqs(int i, double coolness)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  //if there is, the Dirichlet params are common across populations
  for( int j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    for(unsigned s = 0; s < NumStates; ++s)
      temp[s] = DirichletParams[i] + coolness*AlleleCounts[i][s*Populations +j];

    bclib::Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
    for(unsigned s = 0; s < NumStates; ++s){
	if(Freqs[i][j*NumStates+s]==0.0) Freqs[i][j*NumStates+s] = 0.000001;
	if(Freqs[i][j*NumStates+s]==1.0) Freqs[i][j*NumStates+s] = 0.999999;
    }
  }
  delete[] freqs;  
  delete[] temp;
}

///sample prior params using random walk
//NOTE: assuming 2 (allelic) states
void HapMixFreqs::SamplePriorPrecision(unsigned locus, unsigned Populations, double sumlogfreqs1, double sumlogfreqs2){
  try{
    using bclib::lngamma;
    //for( unsigned locus = 0; locus < NumberOfCompositeLoci; ++locus ){
      double LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0;
      const double step = EtaSampler[locus].getStepSize();
      const double current = Eta[locus];
      const double logcurrent = bclib::eh_log(current);
      const double proposal = exp(bclib::Rand::gennor(logcurrent, step));
      //const int NumberOfStates = Loci->GetNumberOfStates(locus);
      //    const double S = (double)NumberOfStates;
      
      //if(proposal > 0.0){
      LogPriorRatio = (EtaShape - 1.0) *( log(proposal) - log(current) )
	- EtaRate *(proposal - current);
      //        LogLikelihoodRatio += Populations *  ( gsl_sf_lngamma(proposal) - gsl_sf_lngamma(current) );
      
      //        LogLikelihoodRatio -= Populations * S * ( gsl_sf_lngamma(proposal / S ) 
      //  							     - gsl_sf_lngamma(current / S ) ) ;
      //        double sumlogfreqs = 0.0;
      //        for(int k = 0; k < Populations; ++k)
      //  	for(int s = 0; s < NumberOfStates; ++s)
      //  	  sumlogfreqs += mylog(Freqs[locus][k*NumberOfStates + s]);
      
      //        LogLikelihoodRatio += (sumlogfreqs / S) * (proposal - current );

      const double mu = DirichletParams[locus]/current;
      LogLikelihoodRatio = (proposal - current)*( mu*sumlogfreqs1 + (1.0-mu)*sumlogfreqs2 )
	+ Populations*(lngamma(proposal) - lngamma(current) -lngamma(mu*proposal) 
	     - lngamma(proposal*(1.0-mu)) + lngamma(mu*current) + lngamma(current*(1.0-mu)));
      
      const double LogAcceptratio = LogLikelihoodRatio + LogPriorRatio;
      if( log(bclib::Rand::myrand()) < LogAcceptratio ){
	Eta[locus] = proposal;
	DirichletParams[locus] *= proposal / current;
	EtaSampler[locus].UpdateStepSize( exp(LogAcceptratio)  );//update stepsize
      }
      //}
      else EtaSampler[locus].UpdateStepSize( 0.0  );//update stepsize
      //    }//end locus loop

  }
  catch(string s){
    throw string ("Error encountered while sampling residual allelic diversity: " + s);
  }
}

void HapMixFreqs::SampleEtaRate(bool afterBurnIn, double sum){
    //sample rate parameter of prior on precision params
    //cout << "Gamma (" << EtaRatePriorShape + (double)NumberOfCompositeLoci * EtaShape << ", " << EtaRatePriorRate + sum << ")" << std::endl;
  EtaRate = bclib::Rand::gengam( EtaRatePriorShape + (double)NumberOfCompositeLoci * EtaShape, 
    			  EtaRatePriorRate + sum);
  if(afterBurnIn) SumLambda += EtaRate;

}

void HapMixFreqs::SamplePriorProportions(unsigned locus, double sumlogfreqs1, double sumlogfreqs2){
  try{
    //for( unsigned locus = 0; locus < NumberOfCompositeLoci; ++locus ){
    //const int NumberOfStates = Loci->GetNumberOfStates(locus);
      HapMixMuArgs.eta = Eta[locus];
      HapMixMuArgs.sumlogfreqs1 = sumlogfreqs1, HapMixMuArgs.sumlogfreqs2 = sumlogfreqs2;
      double mu = MuSampler.Sample(&HapMixMuArgs, d2fmu_hapmix);
      DirichletParams[locus] = mu*Eta[locus];
      //}
  }
  catch(string s){
    throw string ("Error encountered while sampling frequency prior proportions:\n " + s);
  }
}

double HapMixFreqs::fmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (hapmixmuargs*)vargs;
    double f = 0.0;
    const double eta = args->eta;
//loglikelihood
    f += eta*(mu*args->sumlogfreqs1 + (1.0-mu)*args->sumlogfreqs2);
    f -= args->K * (bclib::lngamma(eta*mu) + bclib::lngamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f += lngamma(4) -2*lngamma(2) + log(mu) + log(1.0-mu);

    return f;
}
double HapMixFreqs::dfmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (hapmixmuargs*)vargs;
    double f = 0.0;
    const double eta = args->eta;
//derivative of loglikelihood
    f += eta * (args->sumlogfreqs1 - args->sumlogfreqs2);
    f -= args->K * eta * (bclib::digamma(eta*mu) - bclib::digamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f += 1/(mu) - 1/(1.0-mu);

    return f;
}
double HapMixFreqs::d2fmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (hapmixmuargs*)vargs;
    double f = 0.0;
    const double eta = args->eta;
//2nd derivative of loglikelihood
    f -= args->K * eta * eta * (bclib::trigamma(eta*mu) + bclib::trigamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f -= 1/(mu*mu) - 1/((1.0-mu)*(1.0-mu));

    return f;
}

void HapMixFreqs::OutputErgodicAvg( int samples, std::ofstream *avgstream)const{
  if (SumLambda >0){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumLambda / samples << "\t";
  }
}

// void HapMixFreqs::OutputPriorParams(){
//   if(allelefreqprioroutput.is_open()){
//     OutputPriorParams(allelefreqprioroutput, false);
//   }
// }

void HapMixFreqs::OutputPriorParams(bool tofile, bool toscreen){
  if(DirichletParams){
    const unsigned L = NumberOfCompositeLoci;
    double sumeta = 0.0, sumetasq = 0.0;//, summu = 0.0, summusq = 0.0;
//    double sumobs = 0.0, sumexp = 0.0;
    for(unsigned j = 0; j < L; ++j){
	//const unsigned NumStates = Loci->GetNumberOfStates(j);
	//double sum1 = 0.0;
	//double sumsq1 = 0.0;
	//for(int k = 0; k < Populations; ++k){
	//sum1 += Freqs[j][k*NumStates];//freq allele 1
	// sumsq1 += Freqs[j][k*NumStates]*Freqs[j][k*NumStates];
	//}
	//double obsvar = sumsq1 / (double)(Populations) - (sum1*sum1) / (double)(Populations*Populations);
	//double alpha = DirichletParams[j] ;
	//double beta = Eta[j] - alpha;
	//double expvar = alpha*beta / ((alpha +  beta)*(alpha+beta)*(alpha+beta+ 1.0));

	sumeta += Eta[j];
	sumetasq += Eta[j] *Eta[j];
        //summu += DirichletParams[j]/Eta[j];
        //summusq += (DirichletParams[j]*DirichletParams[j])/(Eta[j]*Eta[j]);
	//sumobs += obsvar;
	//sumexp += expvar;
    }
   
    double meaneta = sumeta / (double) L;
    double vareta = sumetasq / (double)L - meaneta*meaneta;
    //double meanmu = summu / (double) L;
    //double varmu = summusq/ (double) L - meanmu*meanmu;

    if(toscreen){
      std::cout << meaneta << "\t" << vareta; //<< "\t" << meanmu << "\t" << varmu 
      if(etaHierModel)
        std::cout << "\t" << EtaRate;
      std::cout << std::endl;
    }

    if(tofile && allelefreqprioroutput.is_open()){
      allelefreqprioroutput << meaneta << "\t" << vareta;
      // << "\t" << meanmu << "\t" << varmu 
      if(etaHierModel)
        allelefreqprioroutput << "\t" << EtaRate;
      allelefreqprioroutput << std::endl;

    }
    //std::cout << sumobs / (double)L << "\t" << sumexp / (double)L << std::endl;
    //if(tofile && allelefreqprioroutput.is_open())
    //	allelefreqprioroutput <<sumobs / (double)L << "\t" << sumexp / (double)L << std::endl;
  }
}
float HapMixFreqs::getAcceptanceRate()const{
  float sum = 0.0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    sum += EtaSampler[j].getExpectedAcceptanceRate();
  return sum / (float)NumberOfCompositeLoci;
}
float HapMixFreqs::getStepSize()const{
  float sum = 0.0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    sum += EtaSampler[j].getStepSize();
  return sum / (float)NumberOfCompositeLoci;

}

double HapMixFreqs::getParams(unsigned locus)const{
  return DirichletParams[locus];
}

///outputs posterior means of Eta to file
void HapMixFreqs::OutputPosteriorMeans(const char* filename, bclib::LogWriter& Log)const{
  if(IsRandom() && strlen(filename)){
    std::ofstream outfile(filename);
    if(outfile.is_open()){
      Log << bclib::Quiet << "Writing posterior means of residual allelic diversity to " << filename << "\n"; 
      
      for(int j = 0; j < NumberOfCompositeLoci; ++j){
	outfile << SumEta[j] / (double)(NumEtaUpdates) << " ";
      }
      outfile.close();
    }
    else{
      //throw string("Error: cannot open " + filename);
      Log << bclib::On << "Error: cannot open " <<  filename
	  << ". Not writing residual allelic diversity posterior means.\n";
    }
  }
}
///Outputs final values of prior params to file
void HapMixFreqs::OutputFinalValues(const char* filename, bclib::LogWriter& Log)const{
  if(IsRandom() && strlen(filename)){
    std::ofstream outfile(filename);
    if(outfile.is_open()){
      Log << bclib::Quiet << "Writing final values of allele freq prior to " << filename << "\n"; 

      for(int j = 0; j < NumberOfCompositeLoci; ++j){
        outfile << Eta[j] << " " << DirichletParams[j] / Eta[j] << " ";
      }
      outfile.close();
    }
    else{
      //throw string("Error: cannot open " + filename);
      Log << bclib::On << "Error: cannot open " <<  filename
	  << ". Not writing freq prior final values.\n";

    }
  }
}
const FreqArray& HapMixFreqs::getHaploidGenotypeProbs()const{
  return Freqs;
}

const FreqArray& HapMixFreqs::getDiploidGenotypeProbs()const{
  return DiploidGenotypeProbs;
}
