/** 
 *   ADMIXMAP
 *   AdmixFreqs.cc
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
#include "AdmixFreqs.h"
#include "AdmixOptions.h"
#include "AdmixFilenames.h"
#include "InputAdmixData.h"
#include "bclib/AdaptiveRejection.h"
#include "bclib/misc.h"
#include "bclib/linalg.h"
#include "bclib/MuSampler.h"
#include <iomanip>
#include <math.h>
#include <numeric>
#include <string.h>

using bclib::Rand;

double AdmixFreqs::convertValueFromFile(const string s){
  double d = atof(s.c_str());
  if(d < 0.000001) d = 0.000001;//values must be strictly positive
  return d;
}

AdmixFreqs::AdmixFreqs(){

}

AdmixFreqs::~AdmixFreqs(){

}

// ************** Initialisation and loading of data  *******************
void AdmixFreqs::Initialise(Options* const a_options, InputData* const a_data, 
			    Genome *pLoci, bclib::LogWriter &Log, bool MAP){

  AdmixOptions const* options = (AdmixOptions*)a_options;

  // set which sampler will be used for allele freqs
  // current version uses conjugate sampler if annealing without thermo integration
  if( (options->getThermoIndicator() && !options->getTestOneIndivIndicator()) ||
      //using default allele freqs or CAF model
      ( !strlen(options->getAlleleFreqFilename()) &&
	!strlen(options->getHistoricalAlleleFreqFilename()) && 
	!strlen(options->getPriorAlleleFreqFilename()) && 
	!options->getCorrelatedAlleleFreqs() ) ) {
    FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
  } else {
    FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
  }
  
  AlleleFreqs::Initialise(a_options, a_data, pLoci, Log, MAP);
}

void AdmixFreqs::PrintPrior(const Vector_s& , bclib::LogWriter& )const{
}

void AdmixFreqs::LoadAlleleFreqs(Options* const a_options, InputData* const a_data, bclib::LogWriter &Log)
{
  AdmixOptions const* options = (AdmixOptions*)a_options;
  InputAdmixData const* data = (InputAdmixData*)a_data;

  int newrow;
  int row = 0;
  const Matrix_s* temporary = 0;

  int offset = 0;
  bool oldformat = false;//flag for old format allelefreqfile
  bool file = false;//flag to indicate if priors have been supplied in a file

  //old format allelefreqfile
  if( strlen( options->getAlleleFreqFilename() ) ){
    temporary = &(data->getAlleleFreqData());

    offset = 1;
    oldformat = true;
    file = true;
  }
  else if( strlen( options->getPriorAlleleFreqFilename() ) ){
    offset = 0;
    oldformat = false;
    file = true;

    //Prior on DispersionFreqs
    temporary = &(data->getPriorAlleleFreqData());
  }

  if(RandomAlleleFreqs){
    PriorParams = new double*[NumberOfCompositeLoci];//2D array    "      "     otherwise
  }

  //set static members of CompositeLocus
  CompositeLocus::SetRandomAlleleFreqs(RandomAlleleFreqs);
  CompositeLocus::SetNumberOfPopulations(Populations);

  for( unsigned i = 0; i < NumberOfCompositeLoci; i++ ){

    Freqs.array[i] = new double[Loci->GetNumberOfStates(i)* Populations];

    if(file){//read allele freqs from file
	  newrow = row + (*Loci)(i)->GetNumberOfStates() - offset;
	  LoadAlleleFreqs( *temporary, i, row+1, oldformat);//row+1 is the first row for this locus (+1 for the header)
	  row = newrow;
    }
    else {  //set default Allele Freqs
      SetDefaultAlleleFreqs(i);
      if(RandomAlleleFreqs){
//prior for hapmix model is set later in Initialise
	// reference prior on allele freqs: all elements of parameter vector set to 0.5
	// this is unrealistic for large haplotypes - should set all elements to sum to 1
	double defaultpriorparams = 0.5;
	SetDefaultPriorParams(i, defaultpriorparams);
      }
    }
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
void AdmixFreqs::LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool ){
  // set size of allele freqs array for this locus
  // Freqs array has NumberOfStates - 1 elements for each population
  unsigned NumberOfStates = Loci->GetNumberOfStates(i);

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
	  /*New.get( row0+k, j +1)*/ ) / sumalpha;
    }

  if(RandomAlleleFreqs){//no need to allocate remaining arrays if fixed allele freqs model
    PriorParams[i] = new double[NumberOfStates* Populations];

    // priorallelefreqs model, with or without correlated allelefreqs
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
void AdmixFreqs::SetDefaultPriorParams(int i, double defaultpriorparams){
  const int NumberOfStates = Loci->GetNumberOfStates(i);
  PriorParams[i] = new double[NumberOfStates* Populations];
  fill(PriorParams[i], PriorParams[i] + NumberOfStates* Populations, defaultpriorparams);
}

// ************************** Sampling and Updating *****************************************
/// samples allele frequencies and prior parameters.
void AdmixFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  
  AlleleFreqs::Update(IC, afterBurnIn, coolness);  

  if(afterBurnIn){
    //AccumulateKLInfo();
    AccumulateLocusInfo();
  }
}

/** samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
 and stores result in array Freqs 
*/
void AdmixFreqs::SampleAlleleFreqs(int i, double coolness)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  //if there is, the Dirichlet params are common across populations
  for( unsigned j = 0; j < Populations; j++ ){

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


