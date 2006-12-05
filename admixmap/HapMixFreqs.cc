#include "HapMixFreqs.h"
#include "Genome.h"
#include "utils/misc.h"

HapMixFreqs::HapMixFreqs(){
  HapMixPriorParams = 0;
  HapMixPriorEta = 0;
  HapMixPriorEtaSampler = 0;
  SumLambda = 0;
}
HapMixFreqs::~HapMixFreqs(){
  delete[] HapMixPriorParams;
  delete[] HapMixPriorEta;
  delete[] HapMixPriorEtaSampler;
  if( allelefreqprioroutput.is_open()) allelefreqprioroutput.close();
}

void HapMixFreqs::Initialise(unsigned Populations, unsigned L, const std::vector<double> &params){
  NumberOfCompositeLoci = L;
  //allocate prior arrays
  HapMixPriorParams = new double[NumberOfCompositeLoci];//1D array of prior params for hapmixmodel
  HapMixPriorEta = new double[NumberOfCompositeLoci];
  
  HapMixPriorEtaSampler = new StepSizeTuner[NumberOfCompositeLoci];

  //set parameters of prior on frequency Dirichlet prior params
  if(params.size()==3) {
    HapMixPriorShape = params[0];
    HapMixPriorRatePriorShape = params[1];
    HapMixPriorRatePriorRate = params[2];
  }
  else
    {//set defaults
      //TODO: decide on sensible defaults
      HapMixPriorShape = 1.0;
      //      HapMixPriorRatePriorShape = 100.0;
      //HapMixPriorRatePriorRate = 1.0;
    }
  HapMixPriorRate = HapMixPriorRatePriorShape / HapMixPriorRatePriorRate;
  HapMixMuArgs.K = Populations;
  HapMixPriorMuSampler.Initialise(true, true, 1.0, 0.0, fmu_hapmix, dfmu_hapmix );

  for( unsigned i = 0; i < NumberOfCompositeLoci; i++ ){
    HapMixPriorEta[i] = (HapMixPriorShape / HapMixPriorRate);//set dispersion to prior mean
    HapMixPriorParams[i] =  HapMixPriorEta[i] * 0.5;//fixed proportions of 0.5
    HapMixPriorEtaSampler[i].SetParameters(0.1, 0.00001, 100.0, 0.26);
  }
}

void HapMixFreqs::PrintPrior(LogWriter& Log)const{
    Log << "Dirichlet prior on allele frequencies. ";
    Log << "Gamma prior on Dirichlet parameters with shape " << HapMixPriorShape << " and rate " << HapMixPriorRate << ".\n";
    //" and Gamma( " << HapMixPriorRatePriorShape << ", " << HapMixPriorRatePriorRate << " ) prior on rate.\n"; 
}

void HapMixFreqs::OpenOutputFile(const char* filename){
  if(strlen(filename)){
    allelefreqprioroutput.open(filename);
    // allelefreqprioroutput << "eta.Mean\teta.Var\tlambda" << std::endl;
    allelefreqprioroutput << "eta.Mean\teta.Var" << std::endl;
  }
}

///sample prior params using random walk
//NOTE: assuming 2 (allelic) states
void HapMixFreqs::SamplePriorDispersion(unsigned locus, unsigned Populations, double sumlogfreqs1, double sumlogfreqs2){
  //double sum = 0.0;
  try{
    //for( unsigned locus = 0; locus < NumberOfCompositeLoci; ++locus ){
      double LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0;
      const double step = HapMixPriorEtaSampler[locus].getStepSize();
      const double current = HapMixPriorEta[locus];
      const double logcurrent = mylog(current);
      const double proposal = exp(Rand::gennor(logcurrent, step));
      //const int NumberOfStates = Loci->GetNumberOfStates(locus);
      //    const double S = (double)NumberOfStates;
      
      //if(proposal > 0.0){
      LogPriorRatio = (HapMixPriorShape - 1.0) *( log(proposal) - log(current) )
	- HapMixPriorRate *(proposal - current);
      //        LogLikelihoodRatio += Populations *  ( gsl_sf_lngamma(proposal) - gsl_sf_lngamma(current) );
      
      //        LogLikelihoodRatio -= Populations * S * ( gsl_sf_lngamma(proposal / S ) 
      //  							     - gsl_sf_lngamma(current / S ) ) ;
      //        double sumlogfreqs = 0.0;
      //        for(int k = 0; k < Populations; ++k)
      //  	for(int s = 0; s < NumberOfStates; ++s)
      //  	  sumlogfreqs += mylog(Freqs[locus][k*NumberOfStates + s]);
      
      //        LogLikelihoodRatio += (sumlogfreqs / S) * (proposal - current );

      const double mu = HapMixPriorParams[locus]/current;
      LogLikelihoodRatio = (proposal - current)*( mu*sumlogfreqs1 + (1.0-mu)*sumlogfreqs2 )
	+ Populations*(lngamma(proposal) - lngamma(current) -lngamma(mu*proposal) 
	     - lngamma(proposal*(1.0-mu)) + lngamma(mu*current) + lngamma(current*(1.0-mu)));
      
      const double LogAcceptratio = LogLikelihoodRatio + LogPriorRatio;
      if( log(Rand::myrand()) < LogAcceptratio ){
	HapMixPriorEta[locus] = proposal;
	HapMixPriorParams[locus] *= proposal / current;
	HapMixPriorEtaSampler[locus].UpdateStepSize( exp(LogAcceptratio)  );//update stepsize
      }
      //}
      else HapMixPriorEtaSampler[locus].UpdateStepSize( 0.0  );//update stepsize
    // sum += HapMixPriorParams[locus];
      //    }//end locus loop
    //sample rate parameter of prior on prior params
    //cout << "Gamma (" << HapMixPriorRatePriorShape + (double)NumberOfCompositeLoci * HapMixPriorShape << ", " << HapMixPriorRatePriorRate + sum << ")" << std::endl;
    //HapMixPriorRate = Rand::gengam( HapMixPriorRatePriorShape + (double)NumberOfCompositeLoci * HapMixPriorShape, 
    //			  HapMixPriorRatePriorRate + sum);
    //if(afterBurnIn) SumLambda += HapMixPriorRate;
  }
  catch(string s){
    throw string ("Error encountered while sampling frequency prior dispersion: " + s);
  }
}

void HapMixFreqs::SamplePriorProportions(unsigned locus, double sumlogfreqs1, double sumlogfreqs2){
  try{
    //for( unsigned locus = 0; locus < NumberOfCompositeLoci; ++locus ){
    //const int NumberOfStates = Loci->GetNumberOfStates(locus);
      HapMixMuArgs.eta = HapMixPriorEta[locus];
      HapMixMuArgs.sumlogfreqs1 = sumlogfreqs1, HapMixMuArgs.sumlogfreqs2 = sumlogfreqs2;
      double mu = HapMixPriorMuSampler.Sample(&HapMixMuArgs, d2fmu_hapmix);
      HapMixPriorParams[locus] = mu*HapMixPriorEta[locus];
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
  if(HapMixPriorParams){
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
	//double alpha = HapMixPriorParams[j] ;
	//double beta = HapMixPriorEta[j] - alpha;
	//double expvar = alpha*beta / ((alpha +  beta)*(alpha+beta)*(alpha+beta+ 1.0));

	sumeta += HapMixPriorEta[j];
	sumetasq += HapMixPriorEta[j] *HapMixPriorEta[j];
        //summu += HapMixPriorParams[j]/HapMixPriorEta[j];
        //summusq += (HapMixPriorParams[j]*HapMixPriorParams[j])/(HapMixPriorEta[j]*HapMixPriorEta[j]);
	//sumobs += obsvar;
	//sumexp += expvar;
    }
   
    double meaneta = sumeta / (double) L;
    double vareta = sumetasq / (double)L - meaneta*meaneta;
    //double meanmu = summu / (double) L;
    //double varmu = summusq/ (double) L - meanmu*meanmu;
    os << meaneta << "\t" << vareta //<< "\t" << meanmu << "\t" << varmu 
       << /*"\t" << HapMixPriorRate <<*/ std::endl;
    if(tofile && allelefreqprioroutput.is_open())
	allelefreqprioroutput << meaneta << "\t" << vareta// << "\t" << meanmu << "\t" << varmu 
			      << "\t" /*<< HapMixPriorRate << "\t"*/ << std::endl;
    //os << sumobs / (double)L << "\t" << sumexp / (double)L << std::endl;
    //if(tofile && allelefreqprioroutput.is_open())
    //	allelefreqprioroutput <<sumobs / (double)L << "\t" << sumexp / (double)L << std::endl;
  }
}
float HapMixFreqs::getAcceptanceRate()const{
  float sum = 0.0;
  for(unsigned j = 0; j < NumberOfCompositeLoci; ++j)
    sum += HapMixPriorEtaSampler[j].getExpectedAcceptanceRate();
  return sum / (float)NumberOfCompositeLoci;
}
float HapMixFreqs::getStepSize()const{
  float sum = 0.0;
  for(unsigned j = 0; j < NumberOfCompositeLoci; ++j)
    sum += HapMixPriorEtaSampler[j].getStepSize();
  return sum / (float)NumberOfCompositeLoci;

}

double HapMixFreqs::getParams(unsigned locus)const{
  return HapMixPriorParams[locus];
}
