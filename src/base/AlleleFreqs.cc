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
//#include "bclib/AdaptiveRejection.h"
//#include "bclib/misc.h"
//#include <math.h>
#include <string.h>
#include <numeric>
#include "bclib/LogWriter.h"
#include "bclib/Exceptions.h"

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
  SumKLInfo = 0;
  SumLocusInfo = 0;
}

AlleleFreqs::~AlleleFreqs(){
//  if(AlleleFreqsMAP.array != Freqs.array) {
//    AlleleFreqsMAP.dealloc(NumberOfCompositeLoci);
//  }
//  Freqs.dealloc(NumberOfCompositeLoci);
  AlleleCounts.dealloc(NumberOfCompositeLoci);
  hetCounts.dealloc(NumberOfCompositeLoci);

  if(PriorParams) {
    for(unsigned i = 0; i < NumberOfCompositeLoci; ++i){
      if(PriorParams[i])delete[] PriorParams[i];
    }
    delete[] PriorParams;
  }

  if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER)
    for(vector<AlleleFreqSampler*>::const_iterator i = FreqSampler.begin(); i !=FreqSampler.end(); ++i)
      delete *i;
  FreqSampler.clear();

  //delete[] SumKLInfo;
}

void AlleleFreqs::Initialise(Options* const options, InputData* const data, 
			     Genome *pLoci, bclib::LogWriter &Log, bool MAP){

  Loci = pLoci;
  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  //set model indicators
  RandomAlleleFreqs = !options->getFixedAlleleFreqs();

  if(options->getOutputAlleleFreq()){
//     SumKLInfo = new float[NumberOfCompositeLoci];
//     fill(SumKLInfo, SumKLInfo + NumberOfCompositeLoci, 0.0);

    if(Populations > 1){
      SumLocusInfo = new float*[NumberOfCompositeLoci];
      const int NumPopPairs = (Populations*(Populations-1))/2;
      for(unsigned j = 0; j < NumberOfCompositeLoci; ++j){
	SumLocusInfo[j] = new float[NumPopPairs];
	fill(SumLocusInfo[j], SumLocusInfo[j] + NumPopPairs, 0.0);
      }
    }
  }

  //allocate frequency arrays
  Freqs.array = new double*[NumberOfCompositeLoci];
  AlleleFreqsMAP.array = Freqs.array;

  LoadAlleleFreqs(options, data, Log);
  Log.setDisplayMode(bclib::On);
  //open allelefreqoutputfile
  if(IsRandom() ){
    const char* s = options->getAlleleFreqOutputFilename();
    if(strlen(s))
      OpenOutputFile(s);
  }
  
  if(MAP)
    setAlleleFreqsMAP();    
  for( unsigned i = 0; i < NumberOfCompositeLoci; i++ ){
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
    (*Loci)(i)->InitialiseHapPairProbs(Freqs[i], false);
    //If using Chib algorithm, allocate HapPairProbsMAP and copy values in HapPairProbs
    if(MAP) {
      //set AlleleProbsMAP in Composite Locus to point to AlleleFreqsMAP
      (*Loci)(i)->setAlleleProbsMAP(AlleleFreqsMAP[i]);
      (*Loci)(i)->InitialiseHapPairProbsMAP();
    }
  }//end comp locus loop
  

  AllocateAlleleCountArrays(options->getPopulations());

}

void AlleleFreqs::AllocateAlleleCountArrays(unsigned K){
  const int L = Loci->GetNumberOfCompositeLoci();

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

}

///reads initial values of allele freqs from file
///note: no checking of this file is done; ok if produced by function OutputFinalAlleleFreqs
void AlleleFreqs::LoadInitialAlleleFreqs(const char*filename, bclib::LogWriter &Log){
  ifstream infile(filename);
  if(infile.is_open()){
    Log << bclib::Quiet  << "Reading initial values of allele freqs from " << filename << "\n";
    double phi = 0.0;
    unsigned long count = 0;
    for( unsigned locus = 0; locus < NumberOfCompositeLoci; locus++ ){
      const int NumStates = Loci->GetNumberOfStates(locus);

      Freqs.array[locus] = new double[NumStates* Populations];

      for( unsigned pop = 0; pop < Populations; pop++ ){
        Freqs[locus][NumStates-1 + pop*NumStates] = 1.0;
	
	for( int state = 0; state < NumStates-1; state++ ){
	  if(!(infile >> phi)){
	    stringstream ss;
	    ss << "ERROR: Too few entries in initialallelefreqfile. Expected " 
	       << Populations*NumberOfCompositeLoci 
	       << ", found " << count << " " << phi << "\n"; 
	    throw ss.str() ;
	  }
	  ++count;
	  //check 0 <= phi <= 1
	  if(phi < 0 || phi > 1)
	    throw bclib::DataOutOfRangeException("allele frequency", "between 0 and 1", "initialallelefreqfile");
	  if(phi == 1.0)phi = 0.999;
	  if(phi == 0.0)phi = 0.001;
	  Freqs[locus][state + pop*Loci->GetNumberOfStates(locus)] = phi;
	  Freqs[locus][NumStates-1 + pop*NumStates] -= phi;
 	}
      }
    }

   //see if anything more than whitespace left in file
    string test;
    infile >> test;
    if(test.find_first_not_of(" \t\n\r") != string::npos){
      throw string("ERROR: too many entries in initialallelefreqfile\n");
    }
    infile.close();
  }
  else{
    Log << bclib::On << "Error: cannot open " << filename << "\n";
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

 */
void AlleleFreqs::LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool oldformat){
  // set size of allele freqs array for this locus
  // Freqs array has NumberOfStates - 1 elements for each population
  unsigned NumberOfStates = Loci->GetNumberOfStates(i);
  double d = 0.0;

  if(oldformat){//old format allelefreqfile
    for(unsigned k = 0; k < Populations; ++k){
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
    for( unsigned j = 0; j < Populations; j++ ){
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
      for(unsigned col = 0; col < Populations; ++col){
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
  const unsigned LK = NumberOfStates*Populations;
  for(unsigned j = 0; j < LK; ++j)
    Freqs[i][j] = 1.0/NumberOfStates;
}

/// resets all counts to 0
void AlleleFreqs::ResetAlleleCounts(unsigned K) { 
  const int L = Loci->GetNumberOfCompositeLoci();
  for( int i = 0; i < L; i++ ){
    int NumberOfStates = Loci->GetNumberOfStates(i);

    if(AlleleCounts.array)fill(AlleleCounts[i], AlleleCounts[i] + NumberOfStates*K, 0);
    //TODO: only do next line if thermo = 1 and testoneindiv = 0 and annealing
    if(hetCounts.array && NumberOfStates==2)
      fill(hetCounts[i], hetCounts[i]+K*K, 0);
  }
}

/**
 * Updates the counts of alleles observed in each state of ancestry.
 * Given a haplotype pair, h, and the ordered ancestry states at a locus.
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
    for( unsigned i = 0; i < NumberOfCompositeLoci; ++i ){
      FreqSampler[i]->resetStepSizeApproximator(k);
    }
  }
}

/// samples allele frequencies and prior parameters.
void AlleleFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  // Sample allele frequencies conditional on Dirichlet priors 
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized
  for( unsigned i = 0; i < NumberOfCompositeLoci; i++ ){
   const unsigned NumberOfStates = Loci->GetNumberOfStates(i);
   //accumulate summary stats for update of priors
//    double sumlogfreqs1 = 0.0,sumlogfreqs2 = 0.0; 
//     for(unsigned k = 0; k < Populations; ++k){
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

    //no need to update alleleprobs, they are the same as Freqs
    //set HapPair probs using updated alleleprobs
    (*Loci)(i)->SetHapPairProbs();

  }
}

/** does two things:
    1. allocates AlleleFreqsMAP, 
    2. copies current values from Freqs to AlleleFreqsMAP 

    Used in Chib algorithm which just requires a value near the posterior mode. 
*/
void AlleleFreqs::setAlleleFreqsMAP()
{
  if(!AlleleFreqsMAP.array || (AlleleFreqsMAP.array == Freqs.array) ){
    AlleleFreqsMAP.array = new double*[NumberOfCompositeLoci];
  }
  for(unsigned i = 0; i < NumberOfCompositeLoci;++i) {

    AlleleFreqsMAP.array[i] = new double[Loci->GetNumberOfStates(i) * Populations];

    int numstates = Loci->GetNumberOfStates(i);
    for(unsigned k = 0; k < Populations; ++k) {
      for(int j = 0; j < numstates; ++j) {
	AlleleFreqsMAP[i][k*numstates + j] = Freqs[i][k*numstates + j]; //*Loci->GetNumberOfStates(i)];
      } 
    }
  }
}

void AlleleFreqs::AccumulateKLInfo(){
  if(SumKLInfo){
  for(unsigned j = 0; j < NumberOfCompositeLoci; ++j){
    vector<double> psum(Populations, 0);
    double p_entropy = 0.0;
    const int numstates = Loci->GetNumberOfStates(j);
    for(unsigned k = 0; k < Populations; ++k){
      for(int s = 0; s < numstates; ++s){
	double p = Freqs[j][k*numstates + s];
	psum[k] += p;
	p_entropy -= p*log(p);
      }
    }
    double pbar_entropy = 0.0;
    for(unsigned k = 0; k < Populations; ++k){
      psum[k] /= (double) numstates;
      pbar_entropy -= psum[k]*log(psum[k]);
    }
    SumKLInfo[j] += pbar_entropy - p_entropy / (double)Populations;
  }
  }
}

void AlleleFreqs::AccumulateLocusInfo(){
  if(SumLocusInfo && Populations >1){
    for(unsigned j = 0; j < NumberOfCompositeLoci; ++j){
      const int numstates = Loci->GetNumberOfStates(j);
      int k = 0;
      for(unsigned k1 = 0; k1 < Populations; ++k1)
	for(unsigned k2 = k1+1; k2 < Populations; ++k2){
	  
	  float f = 0.0; 
	  for(int s = 0; s < numstates; ++s){
	    f += Freqs[j][k1*numstates+s]*Freqs[j][k2*numstates + s] / (Freqs[j][k1*numstates+s] + Freqs[j][k2*numstates + s]);
	  }
	  SumLocusInfo[j][k] += 1.0 - 2.0*f;
	  ++k;
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
  allelefreqoutput.open(filename);
  allelefreqoutput.setDecimalPrecision(6);
}

//output of freqs as R object (start and finish of file are written elsewhere)
void AlleleFreqs::OutputAlleleFreqs()
{
  if( IsRandom() ){
    for( unsigned locus = 0; locus < NumberOfCompositeLoci; locus++ ){
      const string LocusLabel = "\"" + (*Loci)(locus)->GetLabel(0) + "\"";
      for( int state = 0; state < Loci->GetNumberOfStates(locus)-1; state++ ){
	allelefreqoutput << LocusLabel;
	for( unsigned pop = 0; pop < Populations; pop++ ){
	  allelefreqoutput <<  Freqs[locus][state + pop*Loci->GetNumberOfStates(locus)] ;
	  // allelefreqoutput << "count " << AlleleCounts[locus][state*Populations + pop] <<"; ";
	}
	allelefreqoutput << bclib::newline;
      }
      allelefreqoutput << bclib::newline;
    }
    allelefreqoutput << bclib::newline;
  }
}

//output of freqs to arbitrary file, tab delimited
void AlleleFreqs::OutputAlleleFreqs(const char* filename, bclib::LogWriter& Log)
{
  if(strlen(filename)){
    ofstream outfile(filename);
    outfile.setf(std::ios::fixed); 
    outfile.precision(12);

    if(outfile.is_open()){
      Log << bclib::Quiet << "Writing final values of allele freqs to " << filename << "\n";
      //if( IsRandom() ){
      for( unsigned locus = 0; locus < NumberOfCompositeLoci; locus++ ){
        for( unsigned pop = 0; pop < Populations; pop++ ){
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
      Log << bclib::On << "Error: cannot open " << filename << ", not writing allele freqs.\n";
    }
  }
}

void AlleleFreqs::CloseOutputFile(int iterations, const Vector_s& PopulationLabels)
{
  int nrows = 0;
  for(unsigned j = 0; j < NumberOfCompositeLoci; ++j)
    nrows += Loci->GetNumberOfStates(j)-1;

  vector<vector<string> > dimnames(1);
  dimnames[0].push_back("Locus");
  for (unsigned i = 0; i < Populations; i++){
    dimnames[0].push_back(PopulationLabels[i]);
  }
  vector<int> dims;
  dims.push_back(Populations+1);
  dims.push_back(nrows);
  dims.push_back(iterations);

  allelefreqoutput.close(dims, dimnames);
} 

#include <iomanip>
#include "Filenames.h"
void AlleleFreqs::OutputAlleleFreqSamplerAcceptanceRates(const string& ResultsDir){
  if(FreqSampler.size()){
    const string filename = ResultsDir + "/" + ALLELE_FREQ_SAMPLER_ACCEPT_RATES;
    ofstream file(filename.c_str());
    file.setf(std::ios::fixed); 
    file.precision(3);

    for(vector<AlleleFreqSampler*>::const_iterator i = FreqSampler.begin(); i !=FreqSampler.end(); ++i)
      file << (*i)->getAcceptanceRate() << endl;
    file.close();
  }
}

void AlleleFreqs::WriteKLInfo(unsigned samples, ostream& os){
  if(SumKLInfo){
    os << "Locus\tKLInfo" << endl;
    for(unsigned j = 0; j < NumberOfCompositeLoci; ++j){
      os << "\""<< (*Loci)(j)->GetLabel(0) << "\"\t" << SumKLInfo[j] / (float)samples << endl;
    }
  }
}

void AlleleFreqs::WriteLocusInfo(unsigned samples, const string& ResultsDir, const vector<string>& PopLabels){
  if(SumLocusInfo && Populations >1){
    ofstream os((ResultsDir + "/" + LOCUS_F_VALUES).c_str());

    os << "Locus";
    for(unsigned k1 = 0; k1 < Populations; ++k1)
      for(unsigned k2 = k1+1; k2 < Populations; ++k2){
	os << "\tf." << PopLabels[k1] << "." << PopLabels[k2];
      }
    os << endl;

    os << std::setfill(' ');
    os.setf(std::ios::fixed); 
    os.precision(3);
    os.width(3);

    const unsigned NumPopPairs = (Populations*(Populations-1))/2;
    for(unsigned j = 0; j < NumberOfCompositeLoci; ++j){
      os << "\""<< (*Loci)(j)->GetLabel(0) << "\"";
      for(unsigned k = 0; k < NumPopPairs; ++k)
	os << "\t" << SumLocusInfo[j][k] / (float)samples;
      os << endl;
    }
    os.close();
  }
}

