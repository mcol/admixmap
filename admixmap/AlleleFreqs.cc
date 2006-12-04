/** 
 *   ADMIXMAP
 *   AlleleFreqs.cc 
 *   Class to hold and update allele frequencies, their prior parameters, allele counts and sums. Also holds and updates dispersion
 *   parameter eta and its prior parameters, for a dispersion model. Also computes Fst if required.
 *   Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AlleleFreqs.h"
#include "samplers/AdaptiveRejection.h"
#include "utils/misc.h"
#include "samplers/MuSampler.h"
#include <math.h>
#include <numeric>
#include "Comms.h"
#ifdef PARALLEL
#include <mpe.h>
#endif

//#define DEBUGETA 1

static double convertValueFromFile(const string s){
  double d = atof(s.c_str());
  if(d < 0.000001) d = 0.000001;//values must be strictly positive
  return d;
}

AlleleFreqs::AlleleFreqs(){
  eta = 0;
  psi = 0;
  tau = 0; 
  SumEta = 0;
  psi0 = 0.0;
  RandomAlleleFreqs = false;
  IsHistoricAlleleFreq = false;
  CorrelatedAlleleFreqs = false;
  AlleleCounts.array = 0;
  hetCounts.array = 0;
  Freqs.array = 0;
  AlleleFreqsMAP.array = 0;
  HistoricAlleleFreqs = 0;
  HistoricAlleleCounts = 0;
  PriorParams = 0;
  HapMixPriorParams = 0;
  HapMixPriorEta = 0;
  HapMixPriorEtaSampler = 0;
  SumLambda = 0;
  Fst = 0;
  SumFst = 0;
  calculateFST = false;
  MuProposal = 0;
  TuneEtaSampler = 0;
  w = 1; // frequency of tuning sampler
  NumberOfEtaUpdates  = 0;
  etastep = 0;
#ifdef PARALLEL
  globalAlleleCounts = 0;
  globalHetCounts = 0;
#endif
}

AlleleFreqs::~AlleleFreqs(){
  if(AlleleFreqsMAP.array != Freqs.array) {
    AlleleFreqsMAP.dealloc(NumberOfCompositeLoci);
  }
  Freqs.dealloc(NumberOfCompositeLoci);
  AlleleCounts.dealloc(NumberOfCompositeLoci);
  hetCounts.dealloc(NumberOfCompositeLoci);

  if(PriorParams) {
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      if(PriorParams[i])delete[] PriorParams[i];
    }
    delete[] PriorParams;
  }
  delete[] HapMixPriorParams;
  delete[] HapMixPriorEta;
  delete[] HapMixPriorEtaSampler;

  if(IsHistoricAlleleFreq){
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      delete[] HistoricAlleleCounts[i];
      delete[] HistoricAlleleFreqs[i];
    }
    delete[] HistoricAlleleCounts;
    delete[] HistoricAlleleFreqs;
  }

  delete[] Fst;
  delete[] SumFst;
  delete[] psi;
  delete[] tau;
  delete[] SumEta;
  delete[] MuProposal;

#if ETASAMPLER == 1 
if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs ) {
  delete[] etastep;
  delete[] TuneEtaSampler;
}
#elif ETASAMPLER == 2
if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs ) {
  delete[] etastep;
  delete[] EtaSampler;
}
#endif

#ifdef PARALLEL
  delete[] globalAlleleCounts;
  delete[] globalHetCounts;
  // AlleleFreqArrayType.Free();
#endif
  if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER)
    for(vector<AlleleFreqSampler*>::const_iterator i = FreqSampler.begin(); i !=FreqSampler.end(); ++i)
      delete *i;
  FreqSampler.clear();

  if( allelefreqoutput.is_open()) allelefreqoutput.close();
  if( allelefreqprioroutput.is_open()) allelefreqprioroutput.close();
}

// ************** Initialisation and loading of data  *******************
void AlleleFreqs::Initialise(AdmixOptions* const options, InputData* const data, Genome *pLoci, LogWriter &Log ){
  //initialise Freqs, PriorAlleleFreqs, HistoricAlleleFreqs etc
  Loci = pLoci;
  Populations = options->getPopulations();
  hapmixmodel = options->getHapMixModelIndicator();

  if(Comms::isFreqSampler()){
  LoadAlleleFreqs(options, data);
  Log.setDisplayMode(On);
  //open allelefreqoutputfile
  if(IsRandom() ){
    OpenOutputFile(options);
  }

  // set which sampler will be used for allele freqs
  // current version uses conjugate sampler if annealing without thermo integration
  if( hapmixmodel ||
      (options->getThermoIndicator() && !options->getTestOneIndivIndicator()) ||
      //using default allele freqs or CAF model
      ( !strlen(options->getAlleleFreqFilename()) &&
	!strlen(options->getHistoricalAlleleFreqFilename()) && 
	!strlen(options->getPriorAlleleFreqFilename()) && 
	!options->getCorrelatedAlleleFreqs() ) ) {
    FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
  } else {
    FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
  }

  if(hapmixmodel) {
    //set parameters of prior on frequency Dirichlet prior params
    const vector<double> &params = options->getAlleleFreqPriorParams();
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
  }
  
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    if(RandomAlleleFreqs){
      if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER){
	//set up samplers for allelefreqs
	if(hapmixmodel){
	  FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
						      &(HapMixPriorParams[i]), true));
	}
	else
	  FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
						      PriorParams[i], false));
      }
      if(hapmixmodel){
	HapMixPriorEta[i] = (HapMixPriorShape / HapMixPriorRate);//set dispersion to prior mean
	HapMixPriorParams[i] =  HapMixPriorEta[i] * 0.5;//fixed proportions of 0.5
	HapMixPriorEtaSampler[i].SetParameters(0.1, 0.00001, 100.0, 0.26);
      }

    }
    //set AlleleProbs pointers in CompositeLocus objects to point to Freqs
    //initialise AlleleProbsMAP pointer to 0
    //allocate HapPairProbs and calculate them using AlleleProbs
    (*Loci)(i)->InitialiseHapPairProbs(Freqs[i]);
    //If using Chib algorithm, allocate HapPairProbsMAP and copy values in HapPairProbs
    if(options->getChibIndicator()) {
      // allocate AlleleFreqsMAP and set AlleleProbsMAP in Composite Locus to point to it
      setAlleleFreqsMAP();
      //allocate HapPairProbsMAP
      (*Loci)(i)->InitialiseHapPairProbsMAP();
    }
  }//end comp locus loop
  if(options->getChibIndicator())
      setAlleleFreqsMAP();
  

  // ** settings for sampling of dispersion parameter **
  if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs) {
    unsigned dim = 1;
    if(IsHistoricAlleleFreq) {
      dim = Populations;
    } 
    // ** dispersion parameter(s) and priors **
    eta = new double[ dim ];//dispersion parameters
    psi = new double[ dim ];//gamma prior shape parameter
    tau = new double[ dim ];//gamma prior rate parameter
    SumEta = new double[ dim ];//running sums
    
    //NOTE: if DispersionSampler is used to sample eta in both historic and correlated allele freq models, we don't need to
    //store psi and tau here, just pass the values to the sampler
    
    
    // ** set eta priors **
    if( strlen(options->getEtaPriorFilename()) ){
      //specified by user in file
      Log << "Loading gamma prior parameters for allele frequency dispersion from "
	  << options->getEtaPriorFilename() << ".\n";
      const DataMatrix& etaprior = data->getEtaPriorMatrix();
      
      for( unsigned k = 0; k < dim; k++ ){
	psi[k] = etaprior.get( k, 0 );
	tau[k] = etaprior.get( k, 1 );
      }
    }
    else{//specified by user, otherwise default of psi = 3.0, tau = 0.01
      double rate = options->getEtaMean() / options->getEtaVar();
      double shape = rate *options->getEtaMean();
      
      fill(psi, psi + dim, shape);
      fill(tau, tau + dim, rate);
    }

    //w = 1;
    
    if(IsHistoricAlleleFreq) {
      // ** settings for sampling of proportion vector mu of prior on allele freqs
      muSampler = new MuSampler[dim*NumberOfCompositeLoci]; // 
      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	for(int k = 0; k < Populations; ++k)
	  muSampler[i*Populations+k].setDimensions(2, Loci->GetNumberOfStates(i), 0.001, 0.000001, 10.0, 0.9);
    } else {//correlated allele freq model
      muSampler = new MuSampler[NumberOfCompositeLoci];

      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	muSampler[i].setDimensions(Populations, Loci->GetNumberOfStates(i), 0.002, 0.00001, 10.0, 0.9);

      if(ETASAMPLER==2) {
	EtaSampler = new DispersionSampler[dim];

	int* NumStates;
	NumStates = new int[NumberOfCompositeLoci];
	for(int i = 0; i < NumberOfCompositeLoci; ++i) {
	  NumStates[i] = Loci->GetNumberOfStates(i);
	}
	DispersionSampler::setDimensions(NumberOfCompositeLoci, Populations, NumStates);
	delete[] NumStates;

	EtaSampler[0].Initialise( 0.05, 0.01, 1000.0, 0.9);
	EtaSampler[0].setEtaPrior(psi[0], tau[0]); 

      }
    }
    
    for( unsigned k = 0; k < dim; k++ ){
      //Initialise eta at its prior expectation
      eta[k] = psi[k]/tau[k];
      //Rescale priorallelefreqs so the columns sum to eta 
      for(int j = 0; j < NumberOfCompositeLoci; j++ ){
	int NumberOfStates = Loci->GetNumberOfStates(j);
	double sum = 0.0;
	for(int s = 0; s < NumberOfStates; ++s)sum+= PriorParams[j][k*NumberOfStates+s];
	for(int s = 0; s < NumberOfStates; ++s)PriorParams[j][k*NumberOfStates+s]*= eta[k] / sum;
      }
    }
    
    if(ETASAMPLER==1) {
      etastep0 = 0.1; // sd of proposal distribution for log eta
      etastep = new double[ dim ];
      for(unsigned k = 0; k < dim; ++k) etastep[k] = etastep0;
      NumberOfEtaUpdates = 0;
      TuneEtaSampler = new StepSizeTuner[ dim ];
      for( unsigned k = 0; k < dim; k++ )
	TuneEtaSampler[k].SetParameters( etastep0, 0.01, 10, 0.44 );
    }
    
    // ** Open output file for eta **
    if ( options->getIndAdmixHierIndicator()){
      if (strlen( options->getEtaOutputFilename() ) ){
	InitializeEtaOutputFile(options, data->GetPopLabels(), Log); 
      }
      else{
	Log << "No dispparamfile given\n";
	//exit(1);
      }
    }
    // ** open fst output file if specified **
    if( options->getOutputFST() ){
      OpenFSTFile(options,Log);
    }
  } //end if dispersion parameter
  }
}

void AlleleFreqs::PrintPrior(const Vector_s& PopLabels, LogWriter& Log)const{
  if(IsHistoricAlleleFreq){//historic allele freq (dispersion) model
    Log << "Gamma prior on dispersion parameters with means and variances:\n";
    for( int k = 0; k < Populations; k++ ){
      Log << PopLabels[k] << ": "
	  << psi[k]/tau[k] << "  " << psi[k]/(tau[k]*tau[k]) << "\n";
    }
    Log << "\n";
  }
  else if( CorrelatedAlleleFreqs){//correlated allele freq model
    Log << "Gamma prior on allele frequency dispersion parameter with mean and variance:\n"
	<< psi[0]/tau[0] << "  " << psi[0]/(tau[0]*tau[0]) << "\n";
  }
  else if(hapmixmodel){
    Log << "Dirichlet prior on allele frequencies. ";
    Log << "Gamma prior on Dirichlet parameters with shape " << HapMixPriorShape << " and rate " << HapMixPriorRate << ".\n";
    //" and Gamma( " << HapMixPriorRatePriorShape << ", " << HapMixPriorRatePriorRate << " ) prior on rate.\n"; 
  }
}

void AlleleFreqs::LoadAlleleFreqs(AdmixOptions* const options, InputData* const data_)
{
  int newrow;
  int row = 0;
  const Matrix_s* temporary = 0;
  //set model indicators
  RandomAlleleFreqs = !options->getFixedAlleleFreqs();
  CorrelatedAlleleFreqs = options->getCorrelatedAlleleFreqs();
  //set dimensions

  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();

  //allocate frequency arrays
#ifdef ARRAY2D
  Freqs.array = new double*[NumberOfCompositeLoci];
#else
  Freqs.array = new double[NumberOfCompositeLoci*Populations*2];
  Freqs.stride = Populations*2;
  AlleleFreqsMAP.stride = Populations*2;
#endif
  AlleleFreqsMAP.array = Freqs.array;

  int offset = 0;
  bool oldformat = false;//flag for old format allelefreqfile
  bool file = false;//flag to indicate if priors have been supplied in a file

  //old format allelefreqfile
  if( strlen( options->getAlleleFreqFilename() ) ){
    temporary = &(data_->getAlleleFreqData());

    offset = 1;
    oldformat = true;
    file = true;
  }
  else if( strlen( options->getHistoricalAlleleFreqFilename() ) || strlen( options->getPriorAlleleFreqFilename() ) ){
    offset = 0;
    oldformat = false;
    file = true;
    //Historic AlleleFreqs
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      HistoricAlleleFreqs = new double*[NumberOfCompositeLoci];
      HistoricAlleleCounts = new double*[NumberOfCompositeLoci];
      if( options->getOutputFST() ){
	calculateFST = true;
	Fst = alloc2D_d(NumberOfCompositeLoci, Populations);
	SumFst = alloc2D_d(NumberOfCompositeLoci, Populations);
      }
      MuProposal = new std::vector<StepSizeTuner>[NumberOfCompositeLoci];
      temporary = &(data_->getHistoricalAlleleFreqData());
      IsHistoricAlleleFreq = true;
      
    } else {
      //Prior on AlleleFreqs
      temporary = &(data_->getPriorAlleleFreqData());
      IsHistoricAlleleFreq = false;
    }
  }

  bool useinitfile = false;
  if(RandomAlleleFreqs){
    //allocate prior arrays
    if(hapmixmodel){
      HapMixPriorParams = new double[NumberOfCompositeLoci];//1D array of prior params for hapmixmodel
      HapMixPriorEta = new double[NumberOfCompositeLoci];

      HapMixPriorEtaSampler = new StepSizeTuner[NumberOfCompositeLoci];
      file = false;
      if(strlen( options->getAlleleFreqFilename() )){
	  useinitfile=true;
	  LoadInitialAlleleFreqs(options->getAlleleFreqFilename());
      }
    }
    else
      PriorParams = new double*[NumberOfCompositeLoci];//2D array    "      "     otherwise
  }

  //set static members of CompositeLocus
  CompositeLocus::SetRandomAlleleFreqs(RandomAlleleFreqs);
  CompositeLocus::SetNumberOfPopulations(Populations);

if(!useinitfile)
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
#ifdef ARRAY2D
    Freqs.array[i] = new double[Loci->GetNumberOfStates(i)* Populations];
#endif

    if(file){//read allele freqs from file
	  newrow = row + (*Loci)(i)->GetNumberOfStates() - offset;
	  LoadAlleleFreqs( *temporary, i, row+1, oldformat);//row+1 is the first row for this locus (+1 for the header)
	  row = newrow;
    }
    else {  //set default Allele Freqs
      SetDefaultAlleleFreqs(i);
      if(!hapmixmodel && RandomAlleleFreqs){
//prior for hapmix model is set later in Initialise
	// reference prior on allele freqs: all elements of parameter vector set to 0.5
	// this is unrealistic for large haplotypes - should set all elements to sum to 1
	double defaultpriorparams = 0.5;
	SetDefaultPriorParams(i, defaultpriorparams);
      }
    }
  }
}

void AlleleFreqs::AllocateAlleleCountArrays(unsigned K){
    const int L = Loci->GetNumberOfCompositeLoci();
#ifdef ARRAY2D
  AlleleCounts.array = new int*[L];
  hetCounts.array = new int*[L];

  for( int i = 0; i < L; i++ ){
    // set size of allele counts array
    // allele counts array has NumberOfStates elements for each population 
    AlleleCounts.array[i] = new int[K*Loci->GetNumberOfStates(i)];
    //fill(AlleleCounts[i], AlleleCounts[i]+ Loci->GetNumberOfStates(i)*Populations, 0);
    if(//options->getThermoIndicator() && 
       Loci->GetNumberOfStates(i)==2){//fill hetCounts for SNPs
      hetCounts.array[i] = new int[K * K];
      fill(hetCounts[i], hetCounts[i]+ K*K, 0);
    }
    else hetCounts.array[i] = 0;//set to null pointer for safety
  }

#else
  AlleleCounts.array = new int[L*K*2];
  AlleleCounts.stride = K*2;
  //if(options->getThermoIndicator()){
  hetCounts.array = new int[L*K*K];
  hetCounts.stride = K*K;
  //}
#ifdef PARALLEL
  if(MPI::COMM_WORLD.Get_rank() == 1){
    globalAlleleCounts = new int[L * K*2];
    globalHetCounts = new int[L*K*K];
  }
#endif
#endif
#ifdef PARALLEL
  if(!Freqs.array){//workers need to allocate Freqs
    Freqs.array = new double[L*K*2];
    Freqs.stride = K*2;
  }

#endif
}

///reads initial values of allele freqs from file
///note: no checking of this file is done; ok if produced by function OutputFinalAlleleFreqs
void AlleleFreqs::LoadInitialAlleleFreqs(const char*filename){
    ifstream infile(filename);
    if(infile.is_open()){
	cout << "Reading initial values of allele freqs from " << filename << endl;
    double phi = 0.0;
    for( int locus = 0; locus < NumberOfCompositeLoci; locus++ ){
	const int NumStates = Loci->GetNumberOfStates(locus);
#ifdef ARRAY2D
    Freqs.array[locus] = new double[NumStates* Populations];
#endif
	for( int pop = 0; pop < Populations; pop++ ){
	    Freqs[locus][NumStates-1 + pop*NumStates] = 1.0;
	    for( int state = 0; state < NumStates-1; state++ ){
	    infile >> phi;
	    if(phi == 1.0)phi = 0.999;
	    if(phi == 0.0)phi = 0.001;
	    Freqs[locus][state + pop*Loci->GetNumberOfStates(locus)] = phi;
	    Freqs[locus][NumStates-1 + pop*NumStates] -= phi;
	}
      }
    }
    infile.close();
    }
    else{
	cout << "Error: cannot open " << filename << endl;
	exit(1);
    }
}
/**
 * Initialises the frequencies of each haplotype in the ith
 * composite locus, given Dirichlet priors (oldformat=false) or frequencies (oldformat=true)in matrix New.  
 * If oldformat=false, Allele freqs are set to their prior expectation. 
 * If fixed, allele freqs will be fixed at their prior expectations   
 *
 * New - a matrix containing either frequencies or
 *   parameters for the Dirichlet prior distribution of the allele frequencies. The second dimension is the allele number, 
 *   being in the range of zero to one less than the number of states.
 *   The sum of the prior parameters over all alleles in a population 
 *   (sumalpha) can be interpreted as 
 *   the "prior sample size". The first dimension is the population. Thus, for a 
 *   composite locus with four states and European and African 
 *   populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 9.0 | 3.0 |
 *   State  1 | 3.0 | 4.0 |
 *          2 | 1.0 | 8.0 |
 *          3 | 2.0 | 1.0 |
 *
 * It is also permissable to have one column, in which case the parameters are (initially) the same across all
 * populations. This is required for a correlated allelefreqs model.
 * If Historic is true, sets "historical allele frequencies", where the model has been specified to allow the 
 * allele freqs in the admixed population 
 * to vary from the historical allele frequencies in the unadmixed ancestral populations that have 
 * been sampled. 
 */
void AlleleFreqs::LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool oldformat){
  // set size of allele freqs array for this locus
  // Freqs array has NumberOfStates - 1 elements for each population
  unsigned NumberOfStates = Loci->GetNumberOfStates(i);
  double d = 0.0;

  if(oldformat){//old format allelefreqfile
    for(int k = 0; k < Populations; ++k){
      Freqs[i][(NumberOfStates-1) + k*NumberOfStates] = 1.0;
      for(unsigned j = 0; j < NumberOfStates-1; ++j){
	d = convertValueFromFile(New[row0+j][k+1]);
	Freqs[i][j+ k*NumberOfStates] = d;//New.get(row0+j,k+1);
	Freqs[i][(NumberOfStates-1) + k*NumberOfStates] -= d;//New.get(row0+j,k+1);
      }
    }
  }
  else{
    double sumalpha;
    // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
    // by dividing each prior parameter by the sum of the parameters.
    for( int j = 0; j < Populations; j++ ){
      //vector<double> NewCol = New.getCol(j);
      //sumalpha = accumulate( NewCol.begin(), NewCol.end(), 0.0, std::plus<double>() );
      sumalpha = 0.0;
      for( unsigned k = 0; k < NumberOfStates ; k++ )sumalpha+= convertValueFromFile(New[row0+k][j+1]);//New.get(row0+k, j+1);
      for( unsigned k = 0; k < NumberOfStates ; k++ )
	Freqs[i][ k + j*NumberOfStates ] = ( convertValueFromFile(New[row0+k][(!CorrelatedAlleleFreqs)*j+1] )
	  /*New.get( row0+k, (!CorrelatedAlleleFreqs)* j +1)*/ ) / sumalpha;
    }
  }
  
  if(RandomAlleleFreqs){//no need to allocate remaining arrays if fixed allele freqs model
    PriorParams[i] = new double[NumberOfStates* Populations];
    if(IsHistoricAlleleFreq){
      HistoricAlleleFreqs[i] = new double[(NumberOfStates - 1)* Populations];
      fill(HistoricAlleleFreqs[i],HistoricAlleleFreqs[i] + (NumberOfStates - 1)* Populations, 0.0);
      HistoricAlleleCounts[i] = new double[NumberOfStates* Populations];
      
      for(unsigned row = 0; row < NumberOfStates; ++row)
	for(int col = 0; col < Populations; ++col){
	  d = convertValueFromFile(New[row0+row][col+1]);
	  HistoricAlleleCounts[i][row*Populations +col] = d;//New.get(row0+row, col+1);
	  PriorParams[i][col*NumberOfStates + row] = d/*New.get(row0+row, col+1)*/ + 0.501; // why add 0.501? 
	}
      // set size of vector MuProposal
      if( NumberOfStates > 2 ){
	MuProposal[i].resize( Populations );
	for( int k = 0; k < Populations; k++ ){
	  MuProposal[i][k].SetParameters( 0.01, 0.001, 10.0, 0.23 );
	}
      }
    }
    else{ // priorallelefreqs model, with or without correlated allelefreqs
      RandomAlleleFreqs = true;
      for(unsigned row = 0; row < NumberOfStates; ++row)
	for(int col = 0; col < Populations; ++col){
	  PriorParams[i][col*NumberOfStates + row] = convertValueFromFile(New[row0+row][col+1]);//New.get(row0+row, col+1); 
	}
    }
  }

}

/**
 * Given the number of ancestral populations, sets default values for
 * prior allele freq params.
 */
void AlleleFreqs::SetDefaultPriorParams(int i, double defaultpriorparams){
  int NumberOfStates = Loci->GetNumberOfStates(i);
  if(CorrelatedAlleleFreqs){
    PriorParams[i] = new double[NumberOfStates];
    // reference prior on allele freqs: all elements of parameter vector set to 0.5
    // this is unrealistic for large haplotypes - should set all elements to sum to 1
    fill(PriorParams[i], PriorParams[i] + NumberOfStates, 0.5);
  }
  else{
    PriorParams[i] = new double[NumberOfStates* Populations];
    fill(PriorParams[i], PriorParams[i] + NumberOfStates* Populations, defaultpriorparams);
  }
}

void AlleleFreqs::SetDefaultAlleleFreqs(int i){
  int NumberOfStates = Loci->GetNumberOfStates(i);
  // initialize frequencies as equal for all alleles at locus
  for(int j = 0; j < NumberOfStates*Populations; ++j)
    Freqs[i][j] = 1.0/NumberOfStates;
}

// ************************** Sampling and Updating *****************************************
/// samples allele frequencies and prior parameters.
void AlleleFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness, bool hapmixmodel){
  // Sample prior parameters
  if(HapMixPriorParams){
    SampleHapMixPriorParams();
    SampleHapMixPriorProportions();
    if(afterBurnIn) SumLambda += HapMixPriorRate;
  }
  else
    if(IsHistoricAlleleFreq ){ // TODO: use class DirichletParamSampler 
      for( int i = 0; i < NumberOfCompositeLoci; i++ ){
	if( Loci->GetNumberOfStates(i) == 2 )
	  SampleDirichletParams1D( i);
	else
	  SampleDirichletParamsMultiDim( i);
      }
    }
    else if(CorrelatedAlleleFreqs){
      SampleDirichletParams();
    }
  
  // Sample allele frequencies conditional on Dirichlet priors 
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    const unsigned NumberOfStates = Loci->GetNumberOfStates(i);

    if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER && (NumberOfStates > 2 ||
	accumulate(hetCounts[i], hetCounts[i]+Populations*Populations, 0, std::plus<int>()) > 0 )) {
      if(NumberOfStates==2) //shortcut for SNPs
	FreqSampler[i]->SampleSNPFreqs(Freqs[i], AlleleCounts[i], hetCounts[i], i, Populations, 
				       coolness);
      else FreqSampler[i]->SampleAlleleFreqs(Freqs[i], IC, i, NumberOfStates, Populations, 
					     coolness);
    }
    else //if (FREQSAMPLER==FREQ_CONJUGATE_SAMPLER)
      SampleAlleleFreqs(i, coolness, hapmixmodel);

    if(afterBurnIn)
      (*Loci)(i)->AccumulateAlleleProbs();
#ifndef PARALLEL
    //no need to update alleleprobs, they are the same as Freqs
    //set HapPair probs using updated alleleprobs
    (*Loci)(i)->SetHapPairProbs();
#endif
  }
  
  // Sample for allele frequency dispersion parameters, eta, conditional on allelefreqs using
  // Metropolis random-walk.
  if(  IsHistoricAlleleFreq){ 
    NumberOfEtaUpdates++;
    for( int k = 0; k < Populations; k++ ){
      SampleEtaWithRandomWalk(k, afterBurnIn);
    }
  }
  else if(CorrelatedAlleleFreqs){
    NumberOfEtaUpdates++;
    if(ETASAMPLER == 1) { // random walk metropolis
    SampleEtaWithRandomWalk(0, afterBurnIn);
    } // otherwise eta will be sampled in sampleDirichletParams
      
  }
  if( calculateFST && afterBurnIn && IsHistoricAlleleFreq ){
    UpdateFst();
  }
}

/// resets all counts to 0
void AlleleFreqs::ResetAlleleCounts(unsigned K) { 
    const int L = Loci->GetNumberOfCompositeLoci();
  for( int i = 0; i < L; i++ ){
#ifdef PARALLEL
      int NumberOfStates = 2;
#else
      int NumberOfStates = Loci->GetNumberOfStates(i);
#endif
      if(AlleleCounts.array)fill(AlleleCounts[i], AlleleCounts[i] + NumberOfStates*K, 0);
    //TODO: only do next line if thermo = 1 and testoneindiv = 0 and annealing
      if(hetCounts.array && NumberOfStates==2)
      fill(hetCounts[i], hetCounts[i]+K*K, 0);
  }
#ifdef PARALLEL
      if(globalAlleleCounts)
	  fill(globalAlleleCounts, globalAlleleCounts + L*2*K, 0);
    //TODO: only do next line if thermo = 1 and testoneindiv = 0 and annealing
      if(globalHetCounts)
	  fill(globalHetCounts, globalHetCounts+L*K*K, 0);
#endif
}

/**
   Updates the counts of alleles observed in each state of ancestry.
   Given a haplotype pair, h, and the ordered ancestry states at a locus.
   * should use hap pairs stored in Individual object
   */
void AlleleFreqs::UpdateAlleleCounts(const int locus, const int h[2], const int ancestry[2], const bool diploid, 
				     const bool anneal) {
  if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER ) { 
    if(diploid) {
      if(Loci->GetNumberOfStates(locus)==2) { //diallelic
	if( (h[0] != h[1]) && (ancestry[0] !=ancestry[1])) { //heterozygous with distinct ancestry states
	  ++hetCounts[locus][ancestry[0]*Populations + ancestry[1]];
	} else {
	  ++AlleleCounts[locus][h[0]*Populations + ancestry[0]];
	  ++AlleleCounts[locus][h[1]*Populations + ancestry[1]];
	}
      } else { // > 2 alleles: update allele counts only if not annealed run 
	if(!anneal) { 
	  ++AlleleCounts[locus][ h[0]*Populations + ancestry[0] ];
	  ++AlleleCounts[locus][ h[1]*Populations + ancestry[1] ];
	}
      }
    } else { // haploid: ignore h[1] and ancestry[1]
      if( Loci->GetNumberOfStates(locus)==2 || !anneal) {
	++AlleleCounts[locus][ h[0]*Populations + ancestry[0] ];
      } 
    }
  } else { // conjugate sampler
    ++AlleleCounts[locus][ h[0]*Populations + ancestry[0] ];
    if(diploid) ++AlleleCounts[locus][ h[1]*Populations + ancestry[1] ];
  }
}

void AlleleFreqs::resetStepSizeApproximator(int k) {
  if (RandomAlleleFreqs && FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER) {
    for( int i = 0; i < NumberOfCompositeLoci; ++i ){
     FreqSampler[i]->resetStepSizeApproximator(k);
    }
  }
}


#ifdef PARALLEL
void AlleleFreqs::SumAlleleCountsOverProcesses(unsigned K){
  const int L = Loci->GetNumberOfCompositeLoci();
  Comms::reduceAlleleCounts(AlleleCounts.array, globalAlleleCounts, L*K*2);
  Comms::reduceAlleleCounts(hetCounts.array, globalHetCounts, L*K*K);
}
void AlleleFreqs::BroadcastAlleleFreqs(){
  Comms::BroadcastAlleleFreqs(Freqs.array, (Loci->GetNumberOfCompositeLoci()) * (Freqs.stride));
}
#endif

/** samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
 and stores result in array Freqs 
*/
void AlleleFreqs::SampleAlleleFreqs(int i, double coolness, bool hapmixmodel)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  int c = CorrelatedAlleleFreqs? 0 : 1; //indicates whether correlated allelefreqmodel
  //if there is, the Dirichlet params are common across populations
  for( int j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    if(hapmixmodel)
	for(unsigned s = 0; s < NumStates; ++s)
	    temp[s] = HapMixPriorParams[i] + coolness*AlleleCounts[i][s*Populations +j];
    else
	for(unsigned s = 0; s < NumStates; ++s)
	    temp[s] = PriorParams[i][c*j*NumStates + s] + coolness*AlleleCounts[i][s*Populations +j];

    Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
    for(unsigned s = 0; s < NumStates; ++s){
	if(Freqs[i][j*NumStates+s]==0.0) Freqs[i][j*NumStates+s] = 0.000001;
	if(Freqs[i][j*NumStates+s]==1.0) Freqs[i][j*NumStates+s] = 0.999999;
    }
    // sample HistoricAlleleFreqs as conjugate Dirichlet update with prior specified by PriorParams
    if( IsHistoricAlleleFreq ){
      for(unsigned s = 0; s < NumStates; ++s)
	temp[s] = PriorParams[i][c*j*NumStates + s] + HistoricAlleleCounts[i][s*Populations+j];
      Rand::gendirichlet(NumStates, temp, freqs);
      for(unsigned s = 0; s < NumStates-1; ++s)HistoricAlleleFreqs[i][s*Populations+j] = freqs[s];
    }
  }
  delete[] freqs;  
  delete[] temp;
}

///sample prior params using random walk
//NOTE: assuming 2 (allelic) states
void AlleleFreqs::SampleHapMixPriorParams(){
  //double sum = 0.0;
  try{
    const double K = (double)Populations;
    for( int locus = 0; locus < NumberOfCompositeLoci; ++locus ){
      double LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0;
      const double step = HapMixPriorEtaSampler[locus].getStepSize();
      const double current = HapMixPriorEta[locus];
      const double logcurrent = mylog(current);
    const double proposal = exp(Rand::gennor(logcurrent, step));
    const int NumberOfStates = Loci->GetNumberOfStates(locus);
    //    const double S = (double)NumberOfStates;
    
    //if(proposal > 0.0){
    LogPriorRatio = (HapMixPriorShape - 1.0) *( log(proposal) - log(current) )
      - HapMixPriorRate *(proposal - current);
    //        LogLikelihoodRatio += K *  ( gsl_sf_lngamma(proposal) - gsl_sf_lngamma(current) );
    
    //        LogLikelihoodRatio -= K * S * ( gsl_sf_lngamma(proposal / S ) 
    //  							     - gsl_sf_lngamma(current / S ) ) ;
    //        double sumlogfreqs = 0.0;
    //        for(int k = 0; k < Populations; ++k)
    //  	for(int s = 0; s < NumberOfStates; ++s)
    //  	  sumlogfreqs += mylog(Freqs[locus][k*NumberOfStates + s]);
    
    //        LogLikelihoodRatio += (sumlogfreqs / S) * (proposal - current );
    double sumlogfreqs1 = 0.0,sumlogfreqs2 = 0.0;
    for(int k = 0; k < Populations; ++k){
      sumlogfreqs1 += mylog(Freqs[locus][k*NumberOfStates]);//allele 1
      sumlogfreqs2 += mylog(Freqs[locus][k*NumberOfStates + 1]);//allele 2
      }
    const double mu = HapMixPriorParams[locus]/current;
    LogLikelihoodRatio = (proposal - current)*( mu*sumlogfreqs1 + (1.0-mu)*sumlogfreqs2 )
      + K*(lngamma(proposal) - lngamma(current) -lngamma(mu*proposal) 
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
    }
    //sample rate parameter of prior on prior params
    //cout << "Gamma (" << HapMixPriorRatePriorShape + (double)NumberOfCompositeLoci * HapMixPriorShape << ", " << HapMixPriorRatePriorRate + sum << ")" << endl;
    //HapMixPriorRate = Rand::gengam( HapMixPriorRatePriorShape + (double)NumberOfCompositeLoci * HapMixPriorShape, 
    //			  HapMixPriorRatePriorRate + sum);
  }
  catch(string s){
    throw string ("Error encountered while sampling frequency prior dispersion: " + s);
  }
}

void AlleleFreqs::SampleHapMixPriorProportions(){
    //double sum1 = 0.0, sum2 = 0.0;
  try{
    for( int locus = 0; locus < NumberOfCompositeLoci; ++locus ){
      const int NumberOfStates = Loci->GetNumberOfStates(locus);
      HapMixMuArgs.eta = HapMixPriorEta[locus];
      HapMixMuArgs.sumlogfreqs1 = 0.0, HapMixMuArgs.sumlogfreqs2 = 0.0;
      for(int k = 0; k < Populations; ++k){
	HapMixMuArgs.sumlogfreqs1 += mylog(Freqs[locus][k*NumberOfStates]);//allele 1
	HapMixMuArgs.sumlogfreqs2 += mylog(Freqs[locus][k*NumberOfStates + 1]);//allele 2
      }
      //    sum1 += HapMixMuArgs.sumlogfreqs1;
      //sum2 += HapMixMuArgs.sumlogfreqs2;
      double mu = HapMixPriorMuSampler.Sample(&HapMixMuArgs, d2fmu_hapmix);
      HapMixPriorParams[locus] = mu*HapMixPriorEta[locus];
    }
    //   cout << sum1 / (double)NumberOfCompositeLoci <<" " << sum2/ (double)NumberOfCompositeLoci << endl; 
  }
  catch(string s){
    throw string ("Error encountered while sampling frequency prior proportions:\n " + s);
  }
}

double AlleleFreqs::fmu_hapmix(double mu, const void* const vargs){
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
double AlleleFreqs::dfmu_hapmix(double mu, const void* const vargs){
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
double AlleleFreqs::d2fmu_hapmix(double mu, const void* const vargs){
    const hapmixmuargs* const args = (const hapmixmuargs* const)vargs;
    double f = 0.0;
    const double eta = args->eta;
//2nd derivative of loglikelihood
    f -= args->K * eta * eta * (trigamma(eta*mu) + trigamma(eta*(1.0-mu)));

//log beta(2,2) prior
    //f -= 1/(mu*mu) - 1/((1.0-mu)*(1.0-mu));

    return f;
}


/**
 * Used when model is specified to allow allele freqs in admixed
 * population to vary from the "historical" allele freqs in the
 * unadmixed population.  Dirichlet distribution for allele freqs at
 * locus with k alleles is specified with (k - 1) frequency parameters
 * (mu) and with a single dispersion parameter (eta) this method
 * samples mu and eta, and updates Dirichlet parameters for the allele
 * frequencies
 */
const double *AlleleFreqs::GetStatsForEta( int locus, int population)const
{
  int NumberOfStates = Loci->GetNumberOfStates(locus);
  double *stats = new double[ NumberOfStates ];
  fill(stats, stats+NumberOfStates, 0.0);
  if(IsHistoricAlleleFreq){
    // calculates sufficient stats for update of dispersion parameter
    double sumHistoric = 0.0;
    for( int i = 0; i < NumberOfStates - 1; i++ ){
      stats[ i ] = log( Freqs[locus][ i + population*NumberOfStates ] ) + log( HistoricAlleleFreqs[locus][ i*Populations + population ] );
      sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + population ];
    }
    stats[ NumberOfStates - 1 ] = log( Freqs[locus][NumberOfStates-1 + population*NumberOfStates] ) 
      + log( 1 - sumHistoric );
  }
  else{//correlated allelefreqs model, sum over populations
    for(int k = 0; k < Populations; ++k){
      for( int i = 0; i < NumberOfStates; i++ ){
	stats[ i ] += log( Freqs[locus][ i + k*NumberOfStates ] );
      }
    }
  }

  return stats;
}

///sets PriorParams to sum to eta after sampling of eta
void AlleleFreqs::UpdatePriorAlleleFreqs(int j, const vector<vector<double> >& mu)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      for(int h = 0; h < Loci->GetNumberOfStates(i); ++h)
	PriorParams[i][j*Loci->GetNumberOfStates(i)+h] = mu[i][h];
  }
}

// problem here is to sample the Dirichlet parameters of a multinomial-Dirichlet distribution
// for a multi-allelic locus, we sample the Dirichlet proportion parameters conditional on the  
// dispersion parameter and the allele counts in admixed and historic populations
// by a Metropolis random walk
void AlleleFreqs::SampleDirichletParamsMultiDim( int locus)
{
  int NumberOfStates = Loci->GetNumberOfStates(locus);
  double LogAccProb;
  double* mu1 = new double[NumberOfStates];
  double* mu2 = new double[NumberOfStates];// mu1 is current vector of proportion parameters, mu2 is proposal

  for( int j = 0; j < Populations; j++ ){
    double Proposal1=0, Proposal2=0, f1=0, f2=0;

    for(int s = 0; s< NumberOfStates; ++s)
      mu1[s] = PriorParams[locus][j*NumberOfStates +s] / (eta[j]*MuProposal[locus][j].getStepSize());
    // propose mu2 from Dirichlet distribution with vector of expectations given by mu1
    // step size parameter controls variance: small step size gives small variance

    Rand::gendirichlet(NumberOfStates, mu1, mu2);
 
        
    for( int i = 0; i < NumberOfStates; i++ ){
      mu1[i] *= MuProposal[locus][j].getStepSize();
      // priors on proportion parameters are apparently Dirichlet(1.1, ,,, 1,1) 
      // f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) ); 
      //f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );

      //using Di(1,...,1) prior
      Proposal1 += (eta[j] * mu2[i] - 1) * log( mu1[i] ) - gsl_sf_lngamma( eta[j] * mu2[i] );
      Proposal2 += (eta[j] * mu1[i] - 1) * log( mu2[i] ) - gsl_sf_lngamma( eta[j] * mu1[i] );
    }
      
    for( int k = 0; k < NumberOfStates; k++ ){
      f1 -= 2.0 * gsl_sf_lngamma( mu1[ k ]* eta[j] );
      f2 -= 2.0 * gsl_sf_lngamma( mu2[ k ]* eta[j] );
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
    }
    LogAccProb = f2-f1-Proposal2 + Proposal1;
    if(LogAccProb > 0.0) LogAccProb = 0.0;

    if( log(Rand::myrand()) < LogAccProb ){
      for(int s = 0; s < NumberOfStates; ++s)PriorParams[locus][j*NumberOfStates + s] = mu2[s]*eta[j];
    }
    MuProposal[locus][j].UpdateStepSize(exp(LogAccProb));
  }
  delete[] mu1;
  delete[] mu2;
}

void AlleleFreqs::SampleDirichletParams1D( int locus)
  // with a dispersion model, we sample PriorParams conditional on the observed counts 
  // (with the realized allele freqs integrated out) from a distribution that is 
  // proportional to the product of two binomial-beta likelihoods (for a diallelic locus)
  // we sample the proportion parameter of the beta distribution
  // conditional on the realized allele counts in the admixed population, 
  // the allele counts in the historic population, and the dispersion parameter
  // using an adaptive rejection sampler
  //Note: here NumberOfStates == 2  
{
  double lefttruncation = 0.1;//should be smaller
  MuSamplerArgs MuParameters;
  int* counts0 = new int[2 * Populations];
  double* counts1 = new double[2 * Populations];

  // Construct adaptive rejection sampler for mu.
  for(int i = 0; i < 2; ++i)//loop over the two states/alleles
    for(int j = 0; j < Populations; ++j) {
      counts0[i + j*2] = AlleleCounts[locus][i*Populations +j];
      counts1[i+ j*2] = HistoricAlleleCounts[locus][i*Populations +j];
    }

  AdaptiveRejection SampleMu;
  SampleMu.Initialise(true, true, 1.0, lefttruncation, fMu, dfMu );

  for( int j = 0; j < Populations; j++ ){

    MuParameters.eta = eta[j];
    MuParameters.K = j;
    MuParameters.counts = counts0;
    MuParameters.counts1 = counts1;

    SampleMu.setUpperBound( eta[j] - lefttruncation );//set upper limit for sampler
    //?? wrong  - should have upper limit of 1.0

    PriorParams[locus][ j*2 ] = SampleMu.Sample(&MuParameters, ddfMu); //first state/allele
    //?? wrong should be eta * SampleMu.Sample(...
    // Last (second) prior frequency parameter is determined by sum of mu's = eta.
    PriorParams[locus][ j*2 +1 ] = eta[j] - PriorParams[locus][ j*2 ];
  }
  delete[] counts0;
  delete[] counts1;
}

///sampling of Dirichlet parameters for allelefreqs in historic or correlated allele freq model
void AlleleFreqs::SampleDirichletParams() {
  if(IsHistoricAlleleFreq){
    for(int k = 0; k < Populations; ++k){
      for(int i = 0; i < NumberOfCompositeLoci; ++i){
	int NumberOfStates = Loci->GetNumberOfStates(i);
	int* counts = new int[2*NumberOfStates];
	for(int j = 0; j < NumberOfStates; ++j)
	  {
	    counts[j*NumberOfStates] = AlleleCounts[i][j*Populations +k];
	    counts[j*NumberOfStates+1] = (int)HistoricAlleleCounts[i][j*Populations+k];
	  }
	
	muSampler[i*Populations + k].Sample(PriorParams[i] + (k*NumberOfStates), eta[k], counts);

	//EtaSampler[k].addAlphas(i, PriorParams[i] + (k*NumberOfStates));
	//EtaSampler[k].addCounts(i, counts);
	delete[] counts;
      }

      //sample eta
      //eta[k] = EtaSampler[k].Sample();
      //SumEta[k]+=eta[k];
      //rescale PriorParams so they sum to eta
      //for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      //double sum = accumulate(PriorParams[i][k*NumberOfStates], PriorParams[i][k*NumberOfStates]+NumberOfStates, 
      //			0.0, std::plus<double>());
      //sum = eta[k] / sum;
      //for(int t = 0; t < NumberOfStates; ++t){
      //  PriorParams[i][t] *= sum;
      //}
      //}

    }
  }
  
  else if(CorrelatedAlleleFreqs) {
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      muSampler[i].Sample(PriorParams[i], eta[0], AlleleCounts[i]);
    }

    if(ETASAMPLER == 2) { // hamiltonian sampler
      for(int i = 0; i < NumberOfCompositeLoci; ++i) {
	EtaSampler[0].addAlphas(i, PriorParams[i]);
	EtaSampler[0].addCounts(i, AlleleCounts[i]);
      }
      eta[0] = EtaSampler[0].Sample();
      SumEta[0] += eta[0];
      //rescale PriorParams so they sum to eta
      for( int i = 0; i < NumberOfCompositeLoci; ++i ) {
	int NumberOfStates = Loci->GetNumberOfStates(i);
	double sum = accumulate(PriorParams[i], PriorParams[i]+NumberOfStates, 0.0, 
				std::plus<double>());
	sum = eta[0] / sum;
	for(int t = 0; t < NumberOfStates; ++t) {
	  PriorParams[i][t] *= sum;
	}
      }
    }

  }
}

void AlleleFreqs::SampleEtaWithRandomWalk(int k, bool updateSumEta) {
  double etanew, LogPostRatio = 0.0, LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0, AccProb = 0.0;
  double Denom = 0.0;
  double mineta = 0;
  vector< vector<double> > munew;
  // propose etanew from truncated log-normal distribution.
  do{
    etanew = exp( Rand::gennor( log( eta[k] ), etastep[k] ) );
  }while( etanew > 5000.0 );
  // Prior log-odds ratio (proposal ratio cancels with a part of the prior ratio)   
  LogPriorRatio = ( psi[k] - 1 ) * (log(etanew) - log(eta[k])) - tau[k] * ( etanew - eta[k] );
  // Log-likelihood ratio; numerator of integrating constant
  LogLikelihoodRatio += 2 * NumberOfCompositeLoci * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta[k] ) );
  for(int j = 0; j < NumberOfCompositeLoci; j++ ){
    std::vector<double> mu = GetPriorAlleleFreqs(j,k);
    vector<double>::const_iterator it = min_element(mu.begin(), mu.end());
    
    //mineta is a lower bound for proposal etanew
    if( mineta < 0.1 * eta[k] / *it )
      mineta = 0.1 * eta[k] / *it;
    
    munew.push_back( mu );
    //rescale munew so munew sums to etanew
    for( unsigned l = 0; l < mu.size(); l++ )
      munew[j][l] *= etanew / eta[k];
    
    const double *SumLogFreqs = GetStatsForEta(j,k);
    for( unsigned l = 0; (int)l < (*Loci)(j)->GetNumberOfStates(); l++ ){
      // Denominator of integrating constant
      Denom += 2*(gsl_sf_lngamma( mu[l] ) - gsl_sf_lngamma( munew[j][l] ));
      // SumLogFreqs = log phi_1 + log phi_2
      Denom += (munew[j][l] - mu[l])*SumLogFreqs[l];
    }
    delete[] SumLogFreqs;
  }
  LogPostRatio = LogPriorRatio + LogLikelihoodRatio + Denom;
  // Log acceptance probability = Log posterior ratio since the
  // proposal ratio (log-normal) cancels with prior.
  if(LogPostRatio > 0.0)LogPostRatio = 0.0;
  AccProb = 0.0; 
  if(mineta < etanew )AccProb = xexp(LogPostRatio);
  
#ifdef DEBUGETA
  cout<< "eta["<<k<<"] = "<< eta[k]<<" eta* = "<<etanew<< " stepsize = "<<etastep[k]<<endl;
  cout << "LogPriorRatio= "<< LogPriorRatio <<" LogLikelihoodRatio= " << LogLikelihoodRatio <<endl
       << "Denom = "<< Denom<< " AccProb= exp("<<LogPostRatio<<") = "<< AccProb<<endl<<endl;
#endif
  
  // Acceptance test.
  if( log( Rand::myrand() ) < LogPostRatio && mineta < etanew ){
    eta[k] = etanew;
    UpdatePriorAlleleFreqs( k, munew );
  }
  
  if( !( NumberOfEtaUpdates % w ) ){
    etastep[k] = TuneEtaSampler[k].UpdateStepSize( AccProb );
  }
  
  if( updateSumEta )
    SumEta[k]+=eta[k];
}

/** does three things: 
 1. allocates AlleleFreqsMAP, 
 2. copies current values from Freqs to AlleleFreqsMAP 
 3. sets AlleleProbsMAP in CompositeLocus objects to point to AlleleFreqsMAP

 Used in Chib algorithm which just requires a value near the posterior mode. 
*/
void AlleleFreqs::setAlleleFreqsMAP()
{
  bool allocate = false;
  if(!AlleleFreqsMAP.array || (AlleleFreqsMAP.array == Freqs.array) ){
    allocate = true;
#ifdef ARRAY2D
    AlleleFreqsMAP.array = new double*[NumberOfCompositeLoci];
#else
    AlleleFreqsMAP.array = new double[NumberOfCompositeLoci *Populations];
#endif
  }
  for(int i = 0; i < NumberOfCompositeLoci;++i) {
#ifdef ARRAY2D
    if(allocate) {
      AlleleFreqsMAP.array[i] = new double[Loci->GetNumberOfStates(i) * Populations];
    }
#endif
    int numstates = Loci->GetNumberOfStates(i);
    for(int k = 0; k < Populations; ++k) {
      for(int j = 0; j < numstates; ++j) {
	AlleleFreqsMAP[i][k*numstates + j] = Freqs[i][k*numstates + j]; //*Loci->GetNumberOfStates(i)];
      } 
    }
    (*Loci)(i)->setAlleleProbsMAP(AlleleFreqsMAP[i]);
  }
}

// ************ Accessors **********************************************************

/// Indicates whether allele frequencies are fixed or random.
/// returns a boolean true if allele frequencies are random, false otherwise.
bool AlleleFreqs::IsRandom()const
{
  return( RandomAlleleFreqs );
}

/**
 * Returns Dirichlet parameters for allele frequencies for a particular population and locus.
 * 
 * population - the number of the population (zero based)
 * locus - the number of the locus
 * returns:
 * a vector containing Dirichlet parameters for frequencies of each allele at the locus. 
 * Expected frequencies are calculated by dividing each parameter by the sum of parameters
 */
std::vector<double> AlleleFreqs::GetPriorAlleleFreqs( int locus, int population )const
{
  std::vector<double> counts(Loci->GetNumberOfStates(locus));
  for(int s = 0; s < Loci->GetNumberOfStates(locus); ++s)
    counts[s] = PriorParams[locus][population*Loci->GetNumberOfStates(locus) + s];
  return counts;
}

// AlleleCounts is a 2D array: 1st dimension is locus, 2nd dimension is state*population
// this function returns counts for all populations 
// const int *AlleleFreqs::GetAlleleCounts(int locus)const
// {
//   return( AlleleCounts[locus] );
// }

/// returns counts for the population specified
std::vector<int> AlleleFreqs::GetAlleleCounts( int locus, int population )const
{
  std::vector<int> counts(Loci->GetNumberOfStates(locus));
  for(int s = 0; s < Loci->GetNumberOfStates(locus); ++s)
    counts[s] = AlleleCounts[locus][s*Populations + population];
  return counts;
}

std::vector<double> AlleleFreqs::getAlleleFreqsMAP( int locus, int population )const
{
  int numstates = Loci->GetNumberOfStates(locus);
  std::vector<double> A( numstates );
  for(int i = 0; i < numstates; ++i)
    A[i] = AlleleFreqsMAP[locus][population*numstates + i];
  return A;
}

/**
 * Gets the frequencies of each haplotype in a composite locus.
 *
 * returns:
 * array containing the frequencies of each allele
 * one row per allele state, one column per population
 */
const double* AlleleFreqs::GetAlleleFreqs(int locus)const
{
  return( Freqs[locus] );
}
vector<double> AlleleFreqs::GetAlleleFreqs( int locus, int population )const
{
  vector<double> A(Loci->GetNumberOfStates(locus));
  for(int i = 0; i < Loci->GetNumberOfStates(locus); ++i)
    A[i] = Freqs[locus][i + population*Loci->GetNumberOfStates(locus)];
  return A;
}
const array_of_allelefreqs& AlleleFreqs::GetAlleleFreqs()const{
  return Freqs;
}
// get posterior mode of frequency of allele x, given locus and subpopulation
// double AlleleFreqs::GetAlleleProbsMAP( int x, int ancestry , int locus)const
// {
//   double P;
//   if( x < (*Loci)(locus)->GetNumberOfAllelesOfLocus(0) - 1 )
//     // calculate posterior mode
//     P = AlleleFreqsMAP[locus][ x*Populations+ ancestry ];
//   else // frequency of last allele is set by subtracting sum of posterior modes of other alleles from 1 
//     {
//       P = 1;
//       for( int j = 0; j < (*Loci)(locus)->GetNumberOfAllelesOfLocus(0) - 1; j++ )
// 	P -= AlleleFreqsMAP[locus][ j *Populations+ ancestry ];
//     }
//   return P;
// }

// ******************** Output **************************
void AlleleFreqs::OpenOutputFile(const AdmixOptions* const options)
{
  if(hapmixmodel){
    const char* s = options->getAlleleFreqPriorOutputFilename();
    if(strlen(s)){
      allelefreqprioroutput.open(s, ios::out);
      // allelefreqprioroutput << "eta.Mean\teta.Var\tlambda" << endl;
      allelefreqprioroutput << "eta.Mean\teta.Var" << endl;
    }
  }

  else{
//open file to output freqs as R object
      const char* s = options->getAlleleFreqOutputFilename();
      if(strlen(s) && options->getIndAdmixHierIndicator()){
	  allelefreqoutput.open(s, ios::out );
	  if( !allelefreqoutput){
	      cerr << "Warning: Couldn't open allelefreqoutputfile: " << options->getAlleleFreqOutputFilename() << endl;
	      exit( 1 );
	  }
	  allelefreqoutput << "structure(.Data=c(" << endl;
      }
  }
}
void AlleleFreqs::OutputAlleleFreqs()
{
  if( IsRandom() ){
    for( int locus = 0; locus < NumberOfCompositeLoci; locus++ ){
      for( int state = 0; state < Loci->GetNumberOfStates(locus)-1; state++ ){
	allelefreqoutput << "\"" << (*Loci)(locus)->GetLabel(0) << "\",";
	for( int pop = 0; pop < Populations; pop++ ){
	  allelefreqoutput <<  Freqs[locus][state + pop*Loci->GetNumberOfStates(locus)] << ",";
	  // allelefreqoutput << "count " << AlleleCounts[locus][state*Populations + pop] <<"; ";
	}
	allelefreqoutput << endl;
      }
      allelefreqoutput << endl;
    }
    allelefreqoutput << endl;
  }
}

void AlleleFreqs::OutputAlleleFreqs(const char* filename)
{
    ofstream outfile(filename);
    if(outfile.is_open()){
	cout << "Writing final values of allele freqs to " << filename << endl;
	//if( IsRandom() ){
    for( int locus = 0; locus < NumberOfCompositeLoci; locus++ ){
	for( int pop = 0; pop < Populations; pop++ ){
	    for( int state = 0; state < Loci->GetNumberOfStates(locus)-1; state++ ){
		outfile <<  Freqs[locus][state + pop*Loci->GetNumberOfStates(locus)] << "\t";
//			<< AlleleCounts[locus][state*Populations + pop] << "\t"
		    //		<< AlleleCounts[locus][(state+1)*Populations + pop] << endl;
	}
      }
    }
    //}
  outfile.close();
    }
    else{
	cout << "Error: cannot open " << filename << ", not writing allele freqs." << endl;
    }
}

void AlleleFreqs::OutputPriorParams(){
  //to be used only in hapmixmodel
  if(allelefreqprioroutput.is_open()){
    OutputPriorParams(allelefreqprioroutput, false);
  }
}

void AlleleFreqs::OutputPriorParams(ostream& os, bool tofile){
  //to be used only in hapmixmodel
  if(HapMixPriorParams){
    unsigned L = Loci->GetNumberOfCompositeLoci();
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
       << /*"\t" << HapMixPriorRate <<*/ endl;
    if(tofile && allelefreqprioroutput.is_open())
	allelefreqprioroutput << meaneta << "\t" << vareta// << "\t" << meanmu << "\t" << varmu 
			      << "\t" /*<< HapMixPriorRate << "\t"*/ << endl;
    //os << sumobs / (double)L << "\t" << sumexp / (double)L << endl;
    //if(tofile && allelefreqprioroutput.is_open())
    //	allelefreqprioroutput <<sumobs / (double)L << "\t" << sumexp / (double)L << endl;
  }
}
void AlleleFreqs::CloseOutputFile(int iterations, const Vector_s& PopulationLabels)
{
  int nrows = 0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    nrows += Loci->GetNumberOfStates(j)-1;
  allelefreqoutput << ")," << endl;
  allelefreqoutput << ".Dim = c(";
  allelefreqoutput << Populations+1 << ",";
  allelefreqoutput << nrows << ",";
  allelefreqoutput << iterations;
  allelefreqoutput << ")," << endl;
  allelefreqoutput << ".Dimnames=list(c(\"Locus\",";
  for (int i = 0; i < Populations; i++){
    allelefreqoutput << "\""<<PopulationLabels[i]<<"\"";
    if(i < Populations-1){
      allelefreqoutput << ",";
    }
  }
  allelefreqoutput << "), character(0), character(0)))" << endl;
  allelefreqoutput.close();
} 

void AlleleFreqs::InitializeEtaOutputFile(const AdmixOptions* const options, const Vector_s& PopulationLabels, LogWriter &Log)
{
  Log.setDisplayMode(On);
  outputstream.open( options->getEtaOutputFilename(), ios::out );
  if( !outputstream )
    {
      Log <<  "ERROR: Couldn't open dispparamfile\n";
      exit( 1 );
    }
  else{
    Log << "Writing dispersion parameters to " << options->getEtaOutputFilename() << "\n";
    //Dispersion parameters (eta)
    if(IsHistoricAlleleFreq){
      for( int k = 0; k < Populations; k++ ){
	outputstream << "\"eta." << PopulationLabels[k]<< "\"\t"; 
      }
      outputstream << endl;
    }
    else if(CorrelatedAlleleFreqs){
      outputstream << "\"eta\"" <<endl; 
    }
  }
}

void AlleleFreqs::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  if( IsHistoricAlleleFreq ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta[j] / samples << "\t";
    }
  }
  else if(CorrelatedAlleleFreqs){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumEta[0] / samples << "\t";
  }
  else if (SumLambda >0){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumLambda / samples << "\t";
  }
}

void AlleleFreqs::OutputEta(int iteration, const AdmixOptions *options, LogWriter &Log){
  if( IsHistoricAlleleFreq ){
    //output to logfile
    if( iteration == -1 )
      {
	Log.setDisplayMode(Off);
	Log.setPrecision(6);
	for( int j = 0; j < Populations; j++ ){
	  //Log.width(9);
	  Log << eta[j];
	}
      }
    //output to screen
    if( options->getDisplayLevel()>2 )
      {
	for( int j = 0; j < Populations; j++ ){
	  cout.width(9);
	  cout << setprecision(6) << eta[j] << " ";
	}
      }
    //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      outputstream << eta[0];
      for( int j = 1; j < Populations; j++ ){
	outputstream << "\t" << eta[j];
      }

      outputstream << endl;
    }
  }
  else if(CorrelatedAlleleFreqs){
    
    //output to logfile
    if( iteration == -1 )
      {
	Log.setDisplayMode(Off);
	Log.setPrecision(6);
	Log.width(9);
	Log << eta[0];
      }
    //output to screen
    if( options->getDisplayLevel()>2 )
      {
	cout.width(9);
	cout << setprecision(6) << eta[0] << " ";
      }
    //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      outputstream << eta[0]<< endl;
    }
  }
}

// *** FST functions ****************************
void AlleleFreqs::OpenFSTFile(const AdmixOptions* const options, LogWriter &Log){
  Log.setDisplayMode(On);
  Log << "Writing ergodic averages of FSTs to: " << options->getFSTOutputFilename() << "\n";
  fstoutputstream.open( options->getFSTOutputFilename(), ios::out );
  if( !fstoutputstream ){
    Log << "ERROR: Couldn't open fstoutputfile\n";
    exit( 1 );
  }
}

void AlleleFreqs::UpdateFst()
{
  for(int locus = 0; locus < NumberOfCompositeLoci; ++locus){
    int NumberOfStates = Loci->GetNumberOfStates(locus);
    double q_admix,q_parental,f,H_admix, H_parental, H_combined, pbar;
    for( int k = 0; k < Populations; k++ ){
      H_admix = 0;
      H_parental = 0;
      H_combined = 0;
      q_admix=1.0;
      
      for( int i = 0; i < NumberOfStates - 1; i++ ){
	H_admix += Freqs[locus][ i + k*NumberOfStates ] * Freqs[locus][ i + k*NumberOfStates ];
	H_parental += HistoricAlleleFreqs[locus][ i*Populations+ k ] * HistoricAlleleFreqs[locus][ i*Populations+ k ];
	pbar = 0.5 * ( Freqs[locus][ i + k*NumberOfStates ] + HistoricAlleleFreqs[locus][ i*Populations+ k ] );
	H_combined += pbar * pbar;
	q_admix -= Freqs[locus][i + k*NumberOfStates];
      }
      
      H_admix += q_admix * q_admix;
      double sumHistoric = 0.0;
      for( int i = 0; i < NumberOfStates - 1; i++ ){
	sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + k ];
      }
      q_parental = 1 - sumHistoric;
      H_parental += q_parental * q_parental;
      pbar = 0.5 * ( q_admix + q_parental );
      H_combined += pbar * pbar;
      
      H_combined = 1 - H_combined;
      H_admix = 1 - H_admix;
      H_parental = 1 - H_parental;
      f = ( H_combined - 0.5 * ( H_admix + H_parental ) ) / H_combined;
      Fst[locus][k] = 2*f / ( 1 + f );
    }
    for(int i=0;i < Populations;++i)
      SumFst[locus][i] += Fst[locus][i];
  }
}
void AlleleFreqs::OutputFST(){ 
  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    fstoutputstream << (*Loci)(j)->GetLabel(0);
    for(int k=0; k<Populations; ++k)fstoutputstream << " " << Fst[j][k];
    fstoutputstream << endl;
  }
}

// ********** log density and derivatives for Adaptive Rejection sampling of PriorParams  ***********

double fMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  //counts 0 are counts for admixed pop, counts1 for historic/unadmixed pop
  double f = 0.0, mu = alpha / eta;
  try{
    double logprior = 0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1,1) prior
    f = logprior - 2 * lngamma( alpha ) - 2 * lngamma( eta - alpha );
    
    f += lngamma( alpha+counts0[pop*2] ) + lngamma( eta-alpha+counts0[1+pop*2] );//state 1 + state 2, admixed pop
    f += lngamma( alpha+counts1[pop*2] ) + lngamma( eta-alpha+counts1[1+pop*2] );//historic pop
  }
  catch (string s){
    throw string("Error in ALleleFreqs::fMu - " + s);
  }
  return f;
}

double dfMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x;
  double logprior = 0.1 / alpha - 0.1 / ( eta - alpha );//Beta(1,1) prior
  double f = logprior;
  x = eta - alpha;
  try{
    if(alpha < 0)throw string("arg mu to digamma is negative\n"); 
    if(x < 0)throw string("arg x to digamma is negative\n"); 
    f += 2 * ( digamma(x) - digamma(alpha) );
    
    //admixed pop
    //first state/allele
    x = alpha + counts0[pop*2];
    f += digamma(x);
    //second state/allele
    x = eta - alpha + counts0[1+pop*2];
    f -= digamma(x);
    
    //unadmixed pop
    x = alpha + counts1[pop*2];
    f += digamma(x);
    x = eta - alpha + counts1[1+pop*2];
    f -= digamma(x);
  }
  catch(string s){
    throw string("Error in dfMu in allelefreqs.cc - " + s);
  }
  return f;
}

double ddfMu( double alpha, const void* const args ) {
  // 
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x;
  double prior = -0.1 / (alpha*alpha) - 0.1 / (( eta - alpha ) * ( eta - alpha ) );
  double f = prior;
  try{
    x = eta - alpha;
    f -= 2 * ( trigamma(x) + trigamma(alpha) );
    
    x = alpha + counts0[pop*2];
    f += trigamma(x);
    x = eta - alpha + counts0[1+pop*2];
    f += trigamma(x);
    
    x = alpha + counts1[pop*2];
    f += trigamma(x);
    x = eta - alpha + counts1[1+pop*2];
    f += trigamma(x);
  }
  catch(string s){
    throw string("Error in ddfMu in allelefreqs.cc - " + s);
  }
  return f;
}

float AlleleFreqs::getEtaSamplerAcceptanceRate(int k)const{
  if(ETASAMPLER == 1) {
    return TuneEtaSampler[k].getExpectedAcceptanceRate();
  } else {
    return EtaSampler[k].getAcceptanceRate();
  }
}

float AlleleFreqs::getEtaSamplerStepsize(int k)const{
  if(ETASAMPLER==1) {
    return etastep[k];
  } else {
    return EtaSampler[k].getStepsize();
  }
}

float AlleleFreqs::getAlphaSamplerAcceptanceRate(int i)const{
  return muSampler[i].getAcceptanceRate();
}

float AlleleFreqs::getAlphaSamplerStepsize(int i)const{
  return muSampler[i].getStepsize();
}

float AlleleFreqs::getHapMixPriorSamplerAcceptanceRate()const{
  float sum = 0.0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    sum += HapMixPriorEtaSampler[j].getExpectedAcceptanceRate();
  return sum / (float)NumberOfCompositeLoci;
}
float AlleleFreqs::getHapMixPriorSamplerStepSize()const{
  float sum = 0.0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    sum += HapMixPriorEtaSampler[j].getStepSize();
  return sum / (float)NumberOfCompositeLoci;

}

void AlleleFreqs::OutputAlleleFreqSamplerAcceptanceRates(const char* filename){
  if(FreqSampler.size()){
    ofstream file(filename);
    for(vector<AlleleFreqSampler*>::const_iterator i = FreqSampler.begin(); i !=FreqSampler.end(); ++i)
      file << (*i)->getAcceptanceRate() << endl;
    file.close();
  }
}
