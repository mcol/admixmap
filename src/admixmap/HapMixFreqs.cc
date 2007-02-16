#include "HapMixFreqs.h"
#include "Genome.h"
#include "HapMixOptions.h"
#include "utils/misc.h"
#include "Comms.h"

HapMixFreqs::HapMixFreqs(){
  DirichletParams = 0;
  Eta = 0;
  EtaSampler = 0;
  SumLambda = 0;
  SumEta = 0;
  NumEtaUpdates = 0;
  accumulateEta = false;
  DiploidGenotypeProbs.array = 0;
}
HapMixFreqs::~HapMixFreqs(){
  delete[] DirichletParams;
  delete[] Eta;
  delete[] EtaSampler;
  delete[] SumEta;
  if( allelefreqprioroutput.is_open()) allelefreqprioroutput.close();
  //TODO: delete DiploidGenotypeProbs
}

void HapMixFreqs::Initialise(HapMixOptions* const options, InputData* const data, Genome *pLoci, LogWriter &Log ){
  Loci = pLoci;
  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  //set model indicators
  hapmixmodel = options->getHapMixModelIndicator();
  RandomAlleleFreqs = !options->getFixedAlleleFreqs();
  
  if(Comms::isFreqSampler()){
    LoadAlleleFreqs(options, data, Log);
    Log.setDisplayMode(On);
    //open allelefreqoutputfile
    if(IsRandom() ){
      OpenOutputFile(options->getEtaOutputFilename());
    }
    
    // set which sampler will be used for allele freqs
    // current version uses conjugate sampler if annealing without thermo integration
    if( (options->getThermoIndicator() ) ||
	//using default allele freqs or CAF model
	(  !strlen(options->getPriorAlleleFreqFilename()) && !strlen(options->getInitialAlleleFreqFilename()) ) ) {
      FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
    } else {
      FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
    }
    
    InitialisePrior(Populations, NumberOfCompositeLoci, options, Log );

    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      if(RandomAlleleFreqs){
	if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER){
	  //set up samplers for allelefreqs
	  FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
						      &(DirichletParams[i]), true));
	}
      }
      //set AlleleProbs pointers in CompositeLocus objects to point to Freqs
      //initialise AlleleProbsMAP pointer to 0
      //allocate HapPairProbs and calculate them using AlleleProbs
      (*Loci)(i)->InitialiseHapPairProbs(Freqs[i]);

    }//end comp locus loop
    

  }//end if is freqsampler
  if(Comms::isFreqSampler() || Comms::isWorker()){
    AllocateAlleleCountArrays(options->getPopulations());
#ifdef PARALLEL
    //broadcast initial values of freqs
    BroadcastAlleleFreqs();
#endif
  }
}

void HapMixFreqs::AllocateDiploidGenotypeProbs(){
#ifdef ARRAY2D
  DiploidGenotypeProbs.array = new double*[NumberOfCompositeLoci];
  for(int i = 0; i < NumberOfCompositeLoci; ++i)
     DiploidGenotypeProbs.array[i] = new double[3*Populations*Populations];
#else
  DiploidGenotypeProbs.array = new double[NumberOfCompositeLoci*3*Populations*Populations];
#endif
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

void HapMixFreqs::InitialisePrior(unsigned Populations, unsigned L, const HapMixOptions* const  options, LogWriter& Log ){
  NumberOfCompositeLoci = L;
  //allocate prior arrays
  DirichletParams = new double[NumberOfCompositeLoci];//1D array of prior params for hapmixmodel
  Eta = new double[NumberOfCompositeLoci];
  if(options->OutputAlleleFreqPrior()){
    accumulateEta = true;
    SumEta = new double[NumberOfCompositeLoci];
  }
  
  EtaSampler = new StepSizeTuner[NumberOfCompositeLoci];

  //set parameters of prior on frequency Dirichlet prior params
  const std::vector<double> &params = options->getAlleleFreqPriorParams();
  etaHierModel = options->isFreqDispersionHierModel();
  if(etaHierModel && params.size()!=3)
    Log << On << "ERROR: allelefreqprior should have 3 elements.\n Using default instead.\n";

  if(etaHierModel){
    if(params.size()==3) {
      EtaShape = params[0];
      EtaRatePriorShape = params[1];
      EtaRatePriorRate = params[2];
    }
    else
      {//set defaults
        EtaShape = 1.0;
        EtaRatePriorShape = 0.2;
        EtaRatePriorRate = 1.0;
      }
  }
  else{
      EtaRatePriorRate = 1.0;    
      if(params.size()>=2) {
        EtaShape = params[0];
        EtaRatePriorShape = params[1];
      }
      else
        {//set defaults
          EtaShape = 0.2;
          EtaRatePriorShape = 1.0;
        }
  }

  EtaRate = EtaRatePriorShape / EtaRatePriorRate;
  HapMixMuArgs.K = Populations;
  MuSampler.Initialise(true, true, 1.0, 0.0, fmu_hapmix, dfmu_hapmix );

  std::ifstream initialvaluefile;
  const char* initialvaluefilename = options->getInitialFreqPriorFilename();
  if(strlen(initialvaluefilename)){
    Log << Quiet << "Reading initial values of allele freq prior from " << initialvaluefilename << "\n";
    initialvaluefile.open(initialvaluefilename);
  }

  //TODO: check initial value file has right format, number of entries;
  //TODO??: check prior is consistent with initial values
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    if(initialvaluefile.is_open()){//read values from file
      initialvaluefile >> Eta[i];
      initialvaluefile >> DirichletParams[i];//read mu
      DirichletParams[i] *= Eta[i];         //multiply by eta to get Dir params
    }
    else{
      //set dispersion to prior mean
      Eta[i] = (EtaShape / EtaRate);
      //set proportions to 0.5
      DirichletParams[i] =  Eta[i] * 0.5;
    }
    EtaSampler[i].SetParameters(0.1, 0.00001, 100.0, 0.26);
    if(accumulateEta)SumEta[i] = 0.0;
  }
  if(initialvaluefile.is_open())initialvaluefile.close();
}

void HapMixFreqs::PrintPrior(LogWriter& Log)const{
    Log << "Dirichlet prior on allele frequencies. ";
    Log << "Gamma prior on Dirichlet parameters with shape " << EtaShape << " and rate " << EtaRate << ".\n";
    //" and Gamma( " << EtaRatePriorShape << ", " << EtaRatePriorRate << " ) prior on rate.\n"; 
}

void HapMixFreqs::LoadAlleleFreqs(HapMixOptions* const options, InputData* const data_, LogWriter &Log)
{
  int newrow;
  int row = 0;
  const Matrix_s* temporary = 0;

  //allocate frequency arrays
#ifdef ARRAY2D
  Freqs.array = new double*[NumberOfCompositeLoci];
#else
  Freqs.array = new double[NumberOfCompositeLoci*Populations*2];
  Freqs.stride = Populations*2;
  AlleleFreqsMAP.stride = Populations*2;
#endif
  AlleleFreqsMAP.array = Freqs.array;

  //set static members of CompositeLocus
  CompositeLocus::SetRandomAlleleFreqs(RandomAlleleFreqs);
  CompositeLocus::SetNumberOfPopulations(Populations);

  //read initial values from file
  bool useinitfile = false;
  if(strlen( options->getInitialAlleleFreqFilename() )){
    useinitfile=true;
    LoadInitialAlleleFreqs(options->getInitialAlleleFreqFilename(), Log);
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
#ifdef ARRAY2D
      Freqs.array[i] = new double[Loci->GetNumberOfStates(i)* Populations];
#endif
      
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
void HapMixFreqs::OpenOutputFile(const char* filename){
  if(strlen(filename)){
    allelefreqprioroutput.open(filename);
    allelefreqprioroutput << "eta.Mean\teta.Var";
    if(etaHierModel)allelefreqprioroutput << "\teta.Rate";
    allelefreqprioroutput << std::endl;
  }
}

/// samples allele frequencies and prior parameters.
void HapMixFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  
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
      sumlogfreqs1 += mylog(Freqs[i][k*NumberOfStates]);//allele 1
      sumlogfreqs2 += mylog(Freqs[i][k*NumberOfStates + 1]);//allele 2
    }
    SamplePriorDispersion(i, Populations, sumlogfreqs1, sumlogfreqs2);
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
#ifndef PARALLEL
    //no need to update alleleprobs, they are the same as Freqs
    //set HapPair probs using updated alleleprobs
    (*Loci)(i)->SetHapPairProbs();
#endif
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

    Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
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
void HapMixFreqs::SamplePriorDispersion(unsigned locus, unsigned Populations, double sumlogfreqs1, double sumlogfreqs2){
  try{
    //for( unsigned locus = 0; locus < NumberOfCompositeLoci; ++locus ){
      double LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0;
      const double step = EtaSampler[locus].getStepSize();
      const double current = Eta[locus];
      const double logcurrent = mylog(current);
      const double proposal = exp(Rand::gennor(logcurrent, step));
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
      if( log(Rand::myrand()) < LogAcceptratio ){
	Eta[locus] = proposal;
	DirichletParams[locus] *= proposal / current;
	EtaSampler[locus].UpdateStepSize( exp(LogAcceptratio)  );//update stepsize
      }
      //}
      else EtaSampler[locus].UpdateStepSize( 0.0  );//update stepsize
      //    }//end locus loop

  }
  catch(string s){
    throw string ("Error encountered while sampling frequency prior dispersion: " + s);
  }
}

void HapMixFreqs::SampleEtaRate(bool afterBurnIn, double sum){
    //sample rate parameter of prior on dispersion params
    //cout << "Gamma (" << EtaRatePriorShape + (double)NumberOfCompositeLoci * EtaShape << ", " << EtaRatePriorRate + sum << ")" << std::endl;
  EtaRate = Rand::gengam( EtaRatePriorShape + (double)NumberOfCompositeLoci * EtaShape, 
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
    const hapmixmuargs* const args = (const hapmixmuargs* const)vargs;
    double f = 0.0;
    const double eta = args->eta;
//loglikelihood
    f += eta*(mu*args->sumlogfreqs1 + (1.0-mu)*args->sumlogfreqs2);
    f -= args->K * (lngamma(eta*mu) + lngamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f += lngamma(4) -2*lngamma(2) + log(mu) + log(1.0-mu);

    return f;
}
double HapMixFreqs::dfmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (const hapmixmuargs* const)vargs;
    double f = 0.0;
    const double eta = args->eta;
//derivative of loglikelihood
    f += eta * (args->sumlogfreqs1 - args->sumlogfreqs2);
    f -= args->K * eta * (digamma(eta*mu) - digamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f += 1/(mu) - 1/(1.0-mu);

    return f;
}
double HapMixFreqs::d2fmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (const hapmixmuargs* const)vargs;
    double f = 0.0;
    const double eta = args->eta;
//2nd derivative of loglikelihood
    f -= args->K * eta * eta * (trigamma(eta*mu) + trigamma(eta*(1.0-mu)));

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

void HapMixFreqs::OutputPriorParams(){
  if(allelefreqprioroutput.is_open()){
    OutputPriorParams(allelefreqprioroutput, false);
  }
}

void HapMixFreqs::OutputPriorParams(std::ostream& os, bool tofile){
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

    os << meaneta << "\t" << vareta; //<< "\t" << meanmu << "\t" << varmu 
    if(etaHierModel)
      os << "\t" << EtaRate;
    os << std::endl;

    if(tofile && allelefreqprioroutput.is_open()){
      allelefreqprioroutput << meaneta << "\t" << vareta;
      // << "\t" << meanmu << "\t" << varmu 
      if(etaHierModel)
        allelefreqprioroutput << "\t" << EtaRate;
      allelefreqprioroutput << std::endl;

    }
    //os << sumobs / (double)L << "\t" << sumexp / (double)L << std::endl;
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
void HapMixFreqs::OutputPosteriorMeans(const char* filename, LogWriter& Log)const{
  if(IsRandom() && strlen(filename)){
    std::ofstream outfile(filename);
    if(outfile.is_open()){
      Log << Quiet << "Writing posterior means of allele freq dispersion to " << filename << "\n"; 
      
      for(int j = 0; j < NumberOfCompositeLoci; ++j){
	outfile << SumEta[j] / (double)(NumEtaUpdates) << " ";
      }
      outfile.close();
    }
    else{
      //throw string("Error: cannot open " + filename);
      Log << On << "Error: cannot open " <<  filename
	  << ". Not writing freq dispersion posterior means.\n";
    }
  }
}
///Outputs final values of prior params to file
void HapMixFreqs::OutputFinalValues(const char* filename, LogWriter& Log)const{
  if(IsRandom() && strlen(filename)){
    std::ofstream outfile(filename);
    if(outfile.is_open()){
      Log << Quiet << "Writing final values of allele freq prior to " << filename << "\n"; 

      for(int j = 0; j < NumberOfCompositeLoci; ++j){
        outfile << Eta[j] << " " << DirichletParams[j] / Eta[j] << " ";
      }
      outfile.close();
    }
    else{
      //throw string("Error: cannot open " + filename);
      Log << On << "Error: cannot open " <<  filename
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
