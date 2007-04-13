/** 
 *   AlleleFreqs.cc
 *   Class to hold and update allele frequencies, their prior parameters, allele counts and sums. Also holds and updates dispersion
 *   parameter eta and its prior parameters, for a dispersion model. Also computes Fst if required.
 *   Copyright (c) 2005 - 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AlleleFreqs.h"
//#include "bcppcl/AdaptiveRejection.h"
//#include "bcppcl/misc.h"
//#include <math.h>
#include <numeric>
#include "Comms.h"
#include "bcppcl/LogWriter.h"

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
  RandomAlleleFreqs = false;
  AlleleCounts.array = 0;
  hetCounts.array = 0;
  Freqs.array = 0;
  AlleleFreqsMAP.array = 0;
  PriorParams = 0;

#ifdef PARALLEL
  globalAlleleCounts = 0;
  globalHetCounts = 0;
#endif
}

AlleleFreqs::~AlleleFreqs(){
//  if(AlleleFreqsMAP.array != Freqs.array) {
//    AlleleFreqsMAP.dealloc(NumberOfCompositeLoci);
//  }
//  Freqs.dealloc(NumberOfCompositeLoci);
  AlleleCounts.dealloc(NumberOfCompositeLoci);
  hetCounts.dealloc(NumberOfCompositeLoci);

  if(PriorParams) {
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      if(PriorParams[i])delete[] PriorParams[i];
    }
    delete[] PriorParams;
  }

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
}

// ************** Initialisation and loading of data  *******************
void AlleleFreqs::Initialise(Options* const options, InputData* const data, Genome *pLoci, LogWriter &Log, bool MAP ){
  //initialise Freqs, PriorAlleleFreqs
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
      const char* s = options->getAlleleFreqOutputFilename();
      if(strlen(s))
        OpenOutputFile(s);
    }
    
    // set which sampler will be used for allele freqs
    // current version uses conjugate sampler if annealing without thermo integration
    if( hapmixmodel ||
	(options->getThermoIndicator() && !options->getTestOneIndivIndicator()) ||
	//using default allele freqs or CAF model
	( !strlen(options->getAlleleFreqFilename()) &&
  	  !strlen(options->getPriorAlleleFreqFilename())) ) {
      FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
    } else {
      FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
    }

    if(MAP)
      setAlleleFreqsMAP();    
    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      if(RandomAlleleFreqs){
        if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER){
	  //set up samplers for allelefreqs
          FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
                                                      PriorParams[i], false));
        }
      }
      //set AlleleProbs pointers in CompositeLocus objects to point to Freqs
      //initialise AlleleProbsMAP pointer to 0
      //allocate HapPairProbs and calculate them using AlleleProbs
      (*Loci)(i)->InitialiseHapPairProbs(Freqs[i]);
      //If using Chib algorithm, allocate HapPairProbsMAP and copy values in HapPairProbs
      if(MAP) {
        //set AlleleProbsMAP in Composite Locus to point to AlleleFreqsMAP
        (*Loci)(i)->setAlleleProbsMAP(AlleleFreqsMAP[i]);
    	(*Loci)(i)->InitialiseHapPairProbsMAP();
      }
    }//end comp locus loop

  }//end if isfreqsampler
  if(Comms::isFreqSampler() || Comms::isWorker()){
    AllocateAlleleCountArrays(options->getPopulations());
#ifdef PARALLEL
    //broadcast initial values of freqs
    BroadcastAlleleFreqs();
#endif
  }
}

void AlleleFreqs::PrintPrior(const Vector_s& , LogWriter& )const{
}

void AlleleFreqs::LoadAlleleFreqs(Options* const options, InputData* const data_, LogWriter &Log)
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
  else if( strlen( options->getPriorAlleleFreqFilename() ) ){
    offset = 0;
    oldformat = false;
    file = true;
    //Prior on AlleleFreqs
    temporary = &(data_->getPriorAlleleFreqData());
  }

  bool useinitfile = false;
  if(RandomAlleleFreqs){

    if(hapmixmodel){
 
      file = false;
      if(strlen( options->getAlleleFreqFilename() )){
        useinitfile=true;
        LoadInitialAlleleFreqs(options->getAlleleFreqFilename(), Log);
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
void AlleleFreqs::LoadInitialAlleleFreqs(const char*filename, LogWriter &Log){
  ifstream infile(filename);
  if(infile.is_open()){
    Log << Quiet  << "Reading initial values of allele freqs from " << filename << "\n";
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
    Log << On << "Error: cannot open " << filename << "\n";
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
	Freqs[i][ k + j*NumberOfStates ] = ( convertValueFromFile(New[row0+k][j+1] )
                                             /*New.get( row0+k, (!CorrelatedAlleleFreqs)* j +1)*/ ) / sumalpha;
    }
  }
  
  if(RandomAlleleFreqs){//no need to allocate remaining arrays if fixed allele freqs model
    PriorParams[i] = new double[NumberOfStates* Populations];
    // priorallelefreqs model, with or without correlated allelefreqs
    RandomAlleleFreqs = true;
    for(unsigned row = 0; row < NumberOfStates; ++row)
      for(int col = 0; col < Populations; ++col){
        PriorParams[i][col*NumberOfStates + row] = convertValueFromFile(New[row0+row][col+1]);//New.get(row0+row, col+1); 
      }

  }

}

/**
 * Given the number of ancestral populations, sets default values for
 * prior allele freq params.
 */
void AlleleFreqs::SetDefaultPriorParams(int i, double defaultpriorparams){
  int NumberOfStates = Loci->GetNumberOfStates(i);
  PriorParams[i] = new double[NumberOfStates* Populations];
  fill(PriorParams[i], PriorParams[i] + NumberOfStates* Populations, defaultpriorparams);
}

void AlleleFreqs::SetDefaultAlleleFreqs(int i){
  int NumberOfStates = Loci->GetNumberOfStates(i);
  // initialize frequencies as equal for all alleles at locus
  for(int j = 0; j < NumberOfStates*Populations; ++j)
    Freqs[i][j] = 1.0/NumberOfStates;
}

// ************************** Sampling and Updating *****************************************
/// samples allele frequencies and prior parameters.
void AlleleFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  // Sample prior parameters
  
  // Sample allele frequencies conditional on Dirichlet priors
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    const unsigned NumberOfStates = Loci->GetNumberOfStates(i);
    //accumulate summary stats for update of priors
    //    double sumlogfreqs1 = 0.0,sumlogfreqs2 = 0.0; 
    //     for(int k = 0; k < Populations; ++k){
    // 	sumlogfreqs1 += mylog(Freqs[i][k*NumberOfStates]);//allele 1
    // 	sumlogfreqs2 += mylog(Freqs[i][k*NumberOfStates + 1]);//allele 2
    //       }

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
void AlleleFreqs::SampleAlleleFreqs(int i, double coolness)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  //if there is, the Dirichlet params are common across populations
  for( int j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    for(unsigned s = 0; s < NumStates; ++s)
      temp[s] = PriorParams[i][j*NumStates + s] + coolness*AlleleCounts[i][s*Populations +j];

    Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
    for(unsigned s = 0; s < NumStates; ++s){
      if(Freqs[i][j*NumStates+s]==0.0) Freqs[i][j*NumStates+s] = 0.000001;
      if(Freqs[i][j*NumStates+s]==1.0) Freqs[i][j*NumStates+s] = 0.999999;
    }

  }
  delete[] freqs;  
  delete[] temp;
}

/** does two things:
    1. allocates AlleleFreqsMAP, 
    2. copies current values from Freqs to AlleleFreqsMAP 

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
const FreqArray& AlleleFreqs::GetAlleleFreqs()const{
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
void AlleleFreqs::OpenOutputFile(const char* filename)
{
  //open file to output freqs as R object
  allelefreqoutput.open(filename, ios::out );
  if( !allelefreqoutput){
    cerr << "Warning: Couldn't open allelefreqoutputfile: " << filename << endl;
    exit( 1 );
  }
  allelefreqoutput << "structure(.Data=c(" << endl;
}

//output of freqs as R object (start and finish of file are written elsewhere)
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

//output of freqs to arbitrary file, tab delimited
void AlleleFreqs::OutputAlleleFreqs(const char* filename, LogWriter& Log)
{
  if(strlen(filename)){
    ofstream outfile(filename);
    if(outfile.is_open()){
      Log << Quiet << "Writing final values of allele freqs to " << filename << "\n";
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
      Log << On << "Error: cannot open " << filename << ", not writing allele freqs.\n";
    }
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

void AlleleFreqs::OutputErgodicAvg( int , std::ofstream *)const
{

}

void AlleleFreqs::OutputAlleleFreqSamplerAcceptanceRates(const char* filename){
  if(FreqSampler.size()){
    ofstream file(filename);
    for(vector<AlleleFreqSampler*>::const_iterator i = FreqSampler.begin(); i !=FreqSampler.end(); ++i)
      file << (*i)->getAcceptanceRate() << endl;
    file.close();
  }
}


#ifdef PARALLEL


#endif
