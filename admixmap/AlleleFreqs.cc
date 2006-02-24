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
#include "AdaptiveRejection.h"
#include "functions.h"
#include "MuSampler.h"
#include <math.h>
#include <numeric>

//#define DEBUGETA 1

AlleleFreqs::AlleleFreqs(Genome *pLoci){
  eta = 0;
  psi = 0;
  tau = 0; 
  SumEta = 0;
  psi0 = 0.0;
  Populations = 0;
  NumberOfStates = 0;
  RandomAlleleFreqs = false;
  IsHistoricAlleleFreq = false;
  CorrelatedAlleleFreqs = false;
  AlleleCounts = 0;
  hetCounts = 0;
  Freqs = 0;
  AlleleFreqsMAP = 0;
  HistoricAlleleFreqs = 0;
  HistoricAlleleCounts = 0;
  PriorAlleleFreqs = 0;
  Fst = 0;
  SumFst = 0;
  Loci = pLoci;

  TuneEtaSampler = 0;
  w = 1;
  NumberOfEtaUpdates  = 0;
  etastep = 0;

}

AlleleFreqs::~AlleleFreqs(){
//   for(int i = 0; i < NumberOfCompositeLoci; ++i){
//     //     delete[] Freqs[i];
//     //     delete[] AlleleCounts[i];
//     //     delete[] hetCounts[i];
//      delete[] PriorAlleleFreqs[i];
//      delete[] AlleleFreqsMAP[i];
//      if(IsHistoricAlleleFreq){
//        delete[] HistoricAlleleCounts[i];
//        delete[] HistoricAlleleFreqs[i];
//      }
//    }
  delete[] Freqs;
  delete[] AlleleCounts;
  delete[] hetCounts;
  delete[] PriorAlleleFreqs;
  delete[] HistoricAlleleCounts;
  delete[] AlleleFreqsMAP;
  delete[] HistoricAlleleFreqs;
  delete[] Fst;
  delete[] SumFst;
  delete[] psi;
  delete[] tau;
  delete[] SumEta;
  delete[] NumberOfStates;

#if ETASAMPLER ==1
  delete[] etastep;
  delete[] TuneEtaSampler;
#endif
}

// ************** Initialisation and loading of data  *******************

void AlleleFreqs::Initialise(AdmixOptions* const options, InputData* const data, LogWriter &Log){
  LoadAlleleFreqs(options, data);
  Log.setDisplayMode(On);

  if(IsRandom() &&  options->getOutputAlleleFreq() ){
    OpenOutputFile(options);
  }

  //set up alleleprobs and hap pair probs
  //NB: HaplotypePairProbs in Individual must be set first
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    (*Loci)(i)->Initialise(Freqs[i]);
  }

  // ** settings for sampling of dispersion parameter **
  if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs){
    unsigned dim = 1;
    if(IsHistoricAlleleFreq){
      dim = Populations;

      // ** settings for sampling of PriorAlleleFreqs
      muSampler = new MuSampler[dim*NumberOfCompositeLoci];
      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	for(int k = 0; k < Populations; ++k)
	  muSampler[i*Populations+k].setDimensions(2, NumberOfStates[i], 0.0001, 0.0, 10.0, 0.44);
     }
    else{//correlated allele freq model
      dim = 1;
      muSampler = new MuSampler[NumberOfCompositeLoci];
      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	muSampler[i].setDimensions(Populations, NumberOfStates[i], 0.002, 0.0, 10.0, 0.44);
    }
    
    // ** dispersion parameter(s) and priors **
    eta = new double[ dim ];//dispersion parameters
    psi = new double[ dim ];//gamma prior shape parameter
    tau = new double[ dim ];//gamma prior rate parameter
    SumEta = new double[ dim ];//running sums
    
    //NOTE: if the DispersionSampler is used to sample eta in both historic and correlated allele freq models, we don't need to
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

    if(IsHistoricAlleleFreq){
      Log << "Gamma prior on dispersion parameters with means and variances:\n";
      for( int k = 0; k < Populations; k++ ){
	Log << data->GetPopLabels()[k] << ": "
	    << psi[k]/tau[k] << "  " << psi[k]/(tau[k]*tau[k]) << "\n";
      }
      Log << "\n";
    }
    else{//correlated allele freq model
      Log << "Gamma prior on dispersion parameter with mean and variance:\n"
	  << psi[0]/tau[0] << "  " << psi[0]/(tau[0]*tau[0]) << "\n";
    }
    Log << "\n";
  

    //double maxeta[ dim ];
    for( unsigned k = 0; k < dim; k++ ){
      //       // old method; sets eta to the sum of priorallelefreqs
      //       for( int j = 0; j < NumberOfCompositeLoci; j++ ){
      //        	maxeta[k] =  GetPriorAlleleFreqs(j,k).Sum();
      //        	if( maxeta[k] > eta[k] ){
      //        	  eta[k] = maxeta[k];
      //        	}
      //       }
      
      //Initialise eta at its prior expectation
      eta[k] = psi[k]/tau[k];
      //Rescale priorallelefreqs so the columns sum to eta 
      for(int j = 0; j < NumberOfCompositeLoci; j++ ){
	double sum = 0.0;
	for(int s = 0; s < NumberOfStates[j]; ++s)sum+= PriorAlleleFreqs[j][k*NumberOfStates[j]+s];
	for(int s = 0; s < NumberOfStates[j]; ++s)PriorAlleleFreqs[j][k*NumberOfStates[j]+s]*= eta[k] / sum;
      }
    }

    // ** Settings for random walk sampler
    w = 1;
    etastep0 = 0.1; // sd of proposal distribution for log eta
    etastep = new double[ dim ];
    for(unsigned k = 0; k < dim; ++k) etastep[k] = etastep0;
    NumberOfEtaUpdates = 0;
    TuneEtaSampler = new StepSizeTuner[ dim ];
    for( unsigned k = 0; k < dim; k++ )
      TuneEtaSampler[k].SetParameters( etastep0, 0.01, 10, 0.44 );

//     // ** Settings for Hamiltonian sampler
    //EtaSampler = new DispersionSampler[dim];
    //double initialEtaStepsize = 0.0003;//need a sensible value for this
    //double targetEtaAcceptRate = 0.44;//and this 
    //double min = 0.00;//min and max values for eta stepsize
    //double max = 10.0;

//     if( IsHistoricAlleleFreq)
//       for( unsigned k = 0; k < dim; k++ ){
// 	EtaSampler[k].setDimensions(NumberOfCompositeLoci, 2, NumberOfStates,//<- 2 for admixed & unadmixed pops 
// 				    initialEtaStepsize, min, max, targetEtaAcceptRate);
// 	EtaSampler[k].setEtaPrior(psi[k], tau[k]);
//       }
//       else{//correlated allelefreq model
     //EtaSampler[0].setDimensions(NumberOfCompositeLoci, Populations, NumberOfStates, 
     //			  initialEtaStepsize, min, max, targetEtaAcceptRate);
     //  EtaSampler[0].setEtaPrior(psi[0], tau[0]);
//     }
    //NOTE: only need to set prior if user has not specified

  
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
  }//end if dispersion parameter

}

void AlleleFreqs::OpenOutputFile(const AdmixOptions* const options)
{
  allelefreqoutput.open(options->getAlleleFreqOutputFilename(), ios::out );
  if( !allelefreqoutput && options->getIndAdmixHierIndicator()){
    cerr << "Warning: Couldn't open allelefreqoutputfile: " << options->getAlleleFreqOutputFilename() << endl;
    exit( 1 );
  }
  else{
    allelefreqoutput << "structure(.Data=c(" << endl;
  }
}

void AlleleFreqs::InitializeEtaOutputFile(const AdmixOptions* const options, const std::string* const PopulationLabels, LogWriter &Log)
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

void AlleleFreqs::LoadAlleleFreqs(AdmixOptions* const options, InputData* const data_)
{
  int newrow;
  int row = 0;

  DataMatrix temporary;

  RandomAlleleFreqs = !options->getFixedAlleleFreqs();
  CorrelatedAlleleFreqs = options->getCorrelatedAlleleFreqs();

  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  Freqs = new double*[NumberOfCompositeLoci];
  AlleleFreqsMAP = new double*[NumberOfCompositeLoci];
  HistoricAlleleFreqs = new double*[NumberOfCompositeLoci];
  AlleleCounts = new int*[NumberOfCompositeLoci];
  hetCounts = new int*[NumberOfCompositeLoci];
  HistoricAlleleCounts = new double*[NumberOfCompositeLoci];
  PriorAlleleFreqs = new double*[NumberOfCompositeLoci];
  MuProposal = new std::vector<StepSizeTuner>[NumberOfCompositeLoci];
  NumberOfStates = new int[NumberOfCompositeLoci];

  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    Freqs[i] = 0;
    PriorAlleleFreqs[i] = 0;
    AlleleFreqsMAP[i] = 0;
    HistoricAlleleFreqs[i] = 0;
    AlleleCounts[i] = 0;
    hetCounts[i] = 0;
    HistoricAlleleCounts[i] = 0;
  }
  int offset = 0;
  bool oldformat = false;
  bool file = false;
  //Fixed AlleleFreqs
  if( strlen( options->getAlleleFreqFilename() ) ){
    temporary = data_->getAlleleFreqMatrix();
    temporary = temporary.SubMatrix( 1, temporary.nRows() - 1, 1, Populations );

    offset = 1;
    oldformat = true;
    file = true;
    RandomAlleleFreqs = false;//may be unnecessary; check AdmixOptions
  }
  else if( strlen( options->getHistoricalAlleleFreqFilename() ) || strlen( options->getPriorAlleleFreqFilename() ) ){
    offset = 0;
    oldformat = false;
    file = true;
    //Historic AlleleFreqs
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      temporary = data_->getHistoricalAlleleFreqMatrix();
      IsHistoricAlleleFreq = true;
      
    } else {
      //Prior on AlleleFreqs
      temporary = data_->getPriorAlleleFreqMatrix();
      IsHistoricAlleleFreq = false;
    }
    temporary = temporary.SubMatrix( 1, temporary.nRows() - 1, 1, Populations );
  }

  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    (*Loci)(i)->SetRandomAlleleFreqs(RandomAlleleFreqs);
    (*Loci)(i)->SetNumberOfPopulations(Populations);
    NumberOfStates[i] = (*Loci)(i)->GetNumberOfStates();
    // set size of allele counts array
    // allele counts array has NumberOfStates elements for each population 
    //TODO: may not be necessary to allocate AlleleCounts if fixed allele freqs
    AlleleCounts[i] = new int[NumberOfStates[i] * Populations];
    fill(AlleleCounts[i], AlleleCounts[i]+ NumberOfStates[i]*Populations, 0);
    if(NumberOfStates[i]==2){//fill hetCounts for SNPs
      hetCounts[i] = new int[Populations * Populations];
      fill(hetCounts[i], hetCounts[i]+ Populations*Populations, 0);
    }
    
    if(file){
      newrow = row + (*Loci)(i)->GetNumberOfStates() - offset;
      LoadAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0, Populations - 1 ), i, oldformat);
      row = newrow;
    }
    else{  //Default Allele Freqs
      // reference prior on allele freqs: all elements of parameter vector set to 0.5
      // this is unrealistic for large haplotypes - should set all elements to sum to 1
      double defaultpriorparams = 0.5;
      if(options->getHapMixModelIndicator())defaultpriorparams = 0.1;
      SetDefaultAlleleFreqs(i, defaultpriorparams);
    }
  }

  if(options->getChibIndicator())
    setAlleleFreqsMAP();
}

void AlleleFreqs::LoadAlleleFreqs(const DataMatrix New, int i, bool oldformat){
  /**
   * Initialises the frequencies of each haplotype in the ith
   * composite locus, given Dirichlet priors (oldformat=false) or frequencies (oldformat=true)in matrix New.  
   * If oldformat=false, Allele freqs are set to their prior expectation. 
   * If fixed, allele freqs will be fixed at their prior expectations   
   *
   * New - a matrix containing either frequencies or
   *   parameters for the Dirichlet prior distribution of the allele frequencies. The first dimension is the allele number, 
   *   being in the range of zero to one less than the number of states.
   *   The sum of the prior parameters over all alleles in a population 
   *   (sumalpha) can be interpreted as 
   *   the "prior sample size". The second dimension is the population. Thus, for a 
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


  // set size of allele freqs array for this locus
  // Freqs array has only NumberOfStates - 1 elements for each population
  Freqs[i] = new double[(NumberOfStates[i]-1)* Populations];

  if(oldformat){//old format allelefreqfile
    for(int j = 0; j < NumberOfStates[i]-1; ++j)
      for(int k = 0; k < Populations; ++k)
	Freqs[i][j*Populations+k] = New.get(j,k);
  }
  else{
    double sumalpha;
    // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
    // by dividing each prior parameter by the sum of the parameters.
    for( int j = 0; j < Populations; j++ ){
      vector<double> NewCol = New.getCol(j);
      sumalpha = accumulate( NewCol.begin(), NewCol.end(), 0.0, std::plus<double>() );
      for( int k = 0; k < NumberOfStates[i] - 1; k++ )
	Freqs[i][ k*Populations + j ] = ( New.get( k, (!CorrelatedAlleleFreqs)* j ) ) / sumalpha;
    }
  }
  
  if(RandomAlleleFreqs){
    PriorAlleleFreqs[i] = new double[New.nRows()*New.nCols()];
    if(IsHistoricAlleleFreq){
      HistoricAlleleFreqs[i] = new double[(NumberOfStates[i] - 1)* Populations];
      fill(HistoricAlleleFreqs[i],HistoricAlleleFreqs[i] + (NumberOfStates[i] - 1)* Populations, 0.0);
      HistoricAlleleCounts[i] = new double[New.nRows()*New.nCols()];
      
      for(unsigned row = 0; row < New.nRows(); ++row)
	for(unsigned col = 0; col < New.nCols(); ++col){
	  HistoricAlleleCounts[i][row*New.nCols() +col] = New.get(row, col);
	  PriorAlleleFreqs[i][col*New.nRows() + row] = New.get(row, col) + 0.501; // why add 0.501? 
	}
      
      Fst = alloc2D_d(NumberOfCompositeLoci, Populations);
      SumFst = alloc2D_d(NumberOfCompositeLoci, Populations);
      // set size of vector MuProposal
      if( NumberOfStates[i] > 2 ){
	MuProposal[i].resize( Populations );
	for( int k = 0; k < Populations; k++ ){
	  MuProposal[i][k].SetParameters( 0.01, 0.001, 10.0, 0.23 );
	}
      }
    }
    else{ // priorallelefreqs model, with or without correlated allelefreqs
      RandomAlleleFreqs = true;
      for(unsigned row = 0; row < New.nRows(); ++row)
	for(unsigned col = 0; col < New.nCols(); ++col){
	  PriorAlleleFreqs[i][col*New.nRows() + row] = New.get(row, col); 
	}
    }
  }

}

void AlleleFreqs::SetDefaultAlleleFreqs(int i, double defaultpriorparams){
  /**
   * Given the number of ancestral populations, sets default values for
   * allele frequencies (in Freqs) and prior allele frequencies (in PriorAlleleFreqs).
   */

  if(CorrelatedAlleleFreqs){
    PriorAlleleFreqs[i] = new double[NumberOfStates[i]];
    // reference prior on allele freqs: all elements of parameter vector set to 0.5
    // this is unrealistic for large haplotypes - should set all elements to sum to 1
    fill(PriorAlleleFreqs[i], PriorAlleleFreqs[i] + NumberOfStates[i], 0.5);
  }
  else{
    PriorAlleleFreqs[i] = new double[NumberOfStates[i]* Populations];
    fill(PriorAlleleFreqs[i], PriorAlleleFreqs[i] + NumberOfStates[i]* Populations, defaultpriorparams);
  }
  
  // initialize frequencies as equal for all alleles at locus
  Freqs[i] = new double[(NumberOfStates[i]-1)*Populations];
  for(int j = 0; j < (NumberOfStates[i]-1)*Populations; ++j)
    Freqs[i][j] = 1.0/NumberOfStates[i];
  

}

// ************************** Sampling and Updating *****************************************

// samples allele frequency and prior allele frequency parameters.
void AlleleFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness, bool annealUpdate){
  // Sample for prior frequency parameters mu, using eta, the sum of the frequency parameters for each locus.
  if(IsHistoricAlleleFreq ){
    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      if( NumberOfStates[i] == 2 )
	SampleDirichletParams1D( i);
      else
	SampleDirichletParamsMultiDim( i);
    }
  }
  else if(CorrelatedAlleleFreqs){
    SampleDirichletParams();
  }
  
  // Sample allele frequencies conditional on Dirichlet priors 
  // use these frequencies to set AlleleProbs in CompositeLocus
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    if(annealUpdate){//use long method when computing marginal likelihood
      if(NumberOfStates[i]==2) //shortcut for SNPs
	FreqSampler.SampleSNPFreqs(Freqs[i], PriorAlleleFreqs[i], AlleleCounts[i], hetCounts[i], i, Populations, coolness);
      else FreqSampler.SampleAlleleFreqs(Freqs[i], PriorAlleleFreqs[i], IC, i, NumberOfStates[i], Populations, coolness);
    }
    else //use standard conjugate update
      SampleAlleleFreqs(i, coolness);
    (*Loci)(i)->SetAlleleProbs(Freqs[i], afterBurnIn);
    (*Loci)(i)->SetHapPairProbs();
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
    SampleEtaWithRandomWalk(0, afterBurnIn);
    //if(!(NumberOfEtaUpdates % 100))TuneEtaSampler[0].Reset();
  }
  
  if( afterBurnIn && IsHistoricAlleleFreq ){
    UpdateFst();
  }

}

/*
  Given a haplotype pair, h, and the ordered ancestry states at a locus,
  updates the counts of alleles observed in each state of ancestry.
  * should use hap pairs stored in Individual object
  */
void AlleleFreqs::UpdateAlleleCounts(int locus, const int h[2], const int ancestry[2], bool diploid )
{
  AlleleCounts[locus][ h[0]*Populations + ancestry[0] ]++;
  if(diploid)AlleleCounts[locus][ h[1]*Populations + ancestry[1] ]++;
  //if haploid(ie diploid = false), h[0]==h[1]==genotypes[locus] and ancestry[0]==ancestry[1]
  //and we only count once
}
void AlleleFreqs::UpdateAlleleCounts(int locus, std::vector<unsigned short> genotype, const int ancestry[2], bool diploid )
{//case of SNP when annealing to compute marginal likelihood by thermo method
  if(NumberOfStates[locus]>2)return; //incase called when not a SNP
  if( (genotype[0] != genotype[1]) && (ancestry[0] !=ancestry[1]))//heterozygous with distinct ancestry states
    ++hetCounts[locus][ancestry[0]*Populations + ancestry[1]];
  else{
    ++AlleleCounts[locus][genotype[0]*Populations + ancestry[0]];
    if(diploid)++AlleleCounts[locus][genotype[1]*Populations + ancestry[1]];
  }
  //TODO: check haploid case
}

void AlleleFreqs::SampleAlleleFreqs(int i, double coolness)
{
  // samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
  // and stores result in array Freqs 
  unsigned NumStates = NumberOfStates[i];
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  int c = CorrelatedAlleleFreqs? 0 : 1; //indicates whether correlated allelefreqmodel
  //if there is, the Dirichlet params are common across populations
  for( int j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    for(unsigned s = 0; s < NumStates; ++s)
      temp[s] = PriorAlleleFreqs[i][c*j*NumberOfStates[i] + s] + coolness*AlleleCounts[i][s*Populations +j];
    gendirichlet(NumStates, temp, freqs);
    for(unsigned s = 0; s < NumStates-1; ++s)Freqs[i][s*Populations+j] = freqs[s];

    // sample HistoricAlleleFreqs as a conjugate Dirichlet update with prior specified by PriorAlleleFreqs
    if( IsHistoricAlleleFreq ){
      for(unsigned s = 0; s < NumStates; ++s)
	temp[s] = PriorAlleleFreqs[i][c*j*NumberOfStates[i] + s] + HistoricAlleleCounts[i][s*Populations+j];
      gendirichlet(NumStates, temp, freqs);
      for(unsigned s = 0; s < NumStates-1; ++s)HistoricAlleleFreqs[i][s*Populations+j] = freqs[s];
    }
  }
  delete[] freqs;  
  delete[] temp;
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
  double *stats = new double[ NumberOfStates[locus] ];
  fill(stats, stats+NumberOfStates[locus], 0.0);
  if(IsHistoricAlleleFreq){
    // calculates sufficient stats for update of dispersion parameter
    double sumHistoric = 0.0, sum = 0.0;
    for( int i = 0; i < NumberOfStates[locus] - 1; i++ ){
      stats[ i ] = log( Freqs[locus][ i*Populations+ population ] ) + log( HistoricAlleleFreqs[locus][ i*Populations + population ] );
      sum +=Freqs[locus][ i*Populations + population ];
      sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + population ];
    }
    stats[ NumberOfStates[locus] - 1 ] = log( 1 - sum ) + log( 1 - sumHistoric );
  }
  else{//correlated allelefreqs model, sum over populations
    for(int k = 0; k < Populations; ++k){
      double sum = 0.0;
      for( int i = 0; i < NumberOfStates[locus] - 1; i++ ){
	stats[ i ] += log( Freqs[locus][ i*Populations+ k ] );
	sum +=Freqs[locus][ i*Populations + k ];
      }
      stats[ NumberOfStates[locus] - 1 ] += log( 1 - sum );
    }
  }

  return stats;
}

void AlleleFreqs::UpdatePriorAlleleFreqs(int j, const vector<vector<double> >& mu)
//sets PriorAlleleFreqs to sum to eta after sampling of eta
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      for(int h = 0; h < NumberOfStates[i]; ++h)
	PriorAlleleFreqs[i][j*NumberOfStates[i]+h] = mu[i][h];
  }
}

void AlleleFreqs::SampleDirichletParamsMultiDim( int locus)
  // problem here is to sample the Dirichlet parameters of a multinomial-Dirichlet distribution
  // for a multi-allelic locus, we sample the Dirichlet proportion parameters conditional on the  
  // dispersion parameter and the allele counts in admixed and historic populations
  // by a Metropolis random walk
{

  double LogAccProb;
  double* mu1 = new double[NumberOfStates[locus]];
  double* mu2 = new double[NumberOfStates[locus]];// mu1 is current vector of proportion parameters, mu2 is proposal

  for( int j = 0; j < Populations; j++ ){
    double Proposal1=0, Proposal2=0, f1=0, f2=0;

    for(int s = 0; s< NumberOfStates[locus]; ++s)
      mu1[s] = PriorAlleleFreqs[locus][j*NumberOfStates[locus] +s] / (eta[j]*MuProposal[locus][j].getStepSize());
    // propose mu2 from Dirichlet distribution with vector of expectations given by mu1
    // step size parameter controls variance: small step size gives small variance

    gendirichlet(NumberOfStates[locus], mu1, mu2);
 
        
    for( int i = 0; i < NumberOfStates[locus]; i++ ){
      mu1[i] *= MuProposal[locus][j].getStepSize();
      // priors on proportion parameters are apparently Dirichlet(1.1, ,,, 1,1) 
      // f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) ); 
      //f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );

      //using Di(1,...,1) prior
      Proposal1 += (eta[j] * mu2[i] - 1) * log( mu1[i] ) - gsl_sf_lngamma( eta[j] * mu2[i] );
      Proposal2 += (eta[j] * mu1[i] - 1) * log( mu2[i] ) - gsl_sf_lngamma( eta[j] * mu1[i] );
    }
      
    for( int k = 0; k < NumberOfStates[locus]; k++ ){
      f1 -= 2.0 * gsl_sf_lngamma( mu1[ k ]* eta[j] );
      f2 -= 2.0 * gsl_sf_lngamma( mu2[ k ]* eta[j] );
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
    }
    LogAccProb = f2-f1-Proposal2 + Proposal1;
    if(LogAccProb > 0.0) LogAccProb = 0.0;

    if( log(myrand()) < LogAccProb ){
      for(int s = 0; s < NumberOfStates[locus]; ++s)PriorAlleleFreqs[locus][j*NumberOfStates[locus] + s] = mu2[s]*eta[j];
    }
    MuProposal[locus][j].UpdateStepSize(exp(LogAccProb));
  }
  delete[] mu1;
  delete[] mu2;
}

void AlleleFreqs::SampleDirichletParams1D( int locus)
  // with a dispersion model, we sample PriorAlleleFreqs conditional on the observed counts 
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

    PriorAlleleFreqs[locus][ j*2 ] = SampleMu.Sample(&MuParameters, ddfMu); //first state/allele
    //?? wrong should be eta * SampleMu.Sample(...
    // Last (second) prior frequency parameter is determined by sum of mu's = eta.
    PriorAlleleFreqs[locus][ j*2 +1 ] = eta[j] - PriorAlleleFreqs[locus][ j*2 ];
  }
  delete[] counts0;
  delete[] counts1;
}

//sampling of Dirichlet parameters for allelefreqs in historic or correlated allele freq model
void AlleleFreqs::SampleDirichletParams(){
  if(IsHistoricAlleleFreq){
    for(int k = 0; k < Populations; ++k){
      for(int i = 0; i < NumberOfCompositeLoci; ++i){
	int* counts = new int[2*NumberOfStates[i]];
	for(int j = 0; j < NumberOfStates[i]; ++j)
	  {
	    counts[j*NumberOfStates[i]] = AlleleCounts[i][j*Populations +k];
	    counts[j*NumberOfStates[i]+1] = (int)HistoricAlleleCounts[i][j*Populations+k];
	  }
	
	muSampler[i*Populations + k].Sample(PriorAlleleFreqs[i] + (k*NumberOfStates[i]), eta[k], counts);

	//EtaSampler[k].addAlphas(i, PriorAlleleFreqs[i] + (k*NumberOfStates[i]));
	//EtaSampler[k].addCounts(i, counts);
	delete[] counts;
      }

      //sample eta
      //eta[k] = EtaSampler[k].Sample();
      //SumEta[k]+=eta[k];
      //rescale PriorAlleleFreqs so they sum to eta
      //for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      //double sum = accumulate(PriorAlleleFreqs[i][k*NumberOfStates[i]], PriorAlleleFreqs[i][k*NumberOfStates[i]]+NumberOfStates[i], 
      //			0.0, std::plus<double>());
      //sum = eta[k] / sum;
      //for(int t = 0; t < NumberOfStates[i]; ++t){
      //  PriorAlleleFreqs[i][t] *= sum;
      //}
      //}

    }
  }
  
  else if(CorrelatedAlleleFreqs){
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      muSampler[i].Sample(PriorAlleleFreqs[i], eta[0], AlleleCounts[i]);
      //EtaSampler[0].addAlphas(i, PriorAlleleFreqs[i]);
      //EtaSampler[0].addCounts(i, AlleleCounts[i]);
	}

    //sample eta
    //eta[0] = EtaSampler[0].Sample();
    //SumEta[0] += eta[0];

    //rescale PriorAlleleFreqs so they sum to eta
    //for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    //double sum = accumulate(PriorAlleleFreqs[i], PriorAlleleFreqs[i]+NumberOfStates[i], 0.0, std::plus<double>());
    //sum = eta[0] / sum;
    //for(int t = 0; t < NumberOfStates[i]; ++t){
    //PriorAlleleFreqs[i][t] *= sum;
    //}
    //}
  }
}

void AlleleFreqs::SampleEtaWithRandomWalk(int k, bool updateSumEta){
  double etanew, LogPostRatio = 0.0, LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0, AccProb = 0.0;
  double Denom = 0.0;
  double mineta = 0;
  vector< vector<double> > munew;
  // propose etanew from truncated log-normal distribution.
  do{
    etanew = exp( gennor( log( eta[k] ), etastep[k] ) );
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
  if( log( myrand() ) < LogPostRatio && mineta < etanew ){
    eta[k] = etanew;
    UpdatePriorAlleleFreqs( k, munew );
  }
  
  if( !( NumberOfEtaUpdates % w ) ){
    etastep[k] = TuneEtaSampler[k].UpdateStepSize( AccProb );
  }
  
  if( updateSumEta )
    SumEta[k]+=eta[k];
}

void AlleleFreqs::ResetAlleleCounts() { // resets all counts to 0
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    fill(AlleleCounts[i], AlleleCounts[i]+NumberOfStates[i]*Populations, 0);
    //TODO: only do next line if thermo = 1 and testoneindiv = 0 and annealing
    if(NumberOfStates[i]==2)
      fill(hetCounts[i], hetCounts[i]+Populations*Populations, 0);
  }

}

// sets posterior mode of allelefreqs as current value
// probably doesn't matter where there is strong prior on allele freqs, and used in Chib algorithm 
// which just requires a value near the posterior mode
// set AlleleFreqsMAP and getAlleleFreqs are called by Individual object
void AlleleFreqs::setAlleleFreqsMAP()
{
  for(int i = 0; i < NumberOfCompositeLoci;++i){
    if(AlleleFreqsMAP[i]==0)
      AlleleFreqsMAP[i] = new double[(NumberOfStates[i]-1)*Populations];
    for(int j = 0; j < NumberOfStates[i] - 1; ++j)for(int k = 0; k < Populations; ++k)
      AlleleFreqsMAP[i][j*Populations+k] = Freqs[i][j*Populations+k];
  }
}

// ************ Accessors **********************************************************

// Indicates whether allele frequencies are fixed or random.
// returns a boolean true if allele frequencies are random, false otherwise.
bool AlleleFreqs::IsRandom()const
{
  return( RandomAlleleFreqs );
}

/**PriorAlleleFreqs
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
  std::vector<double> counts(NumberOfStates[locus]);
  for(int s = 0; s < NumberOfStates[locus]; ++s)
    counts[s] = PriorAlleleFreqs[locus][population*NumberOfStates[locus] + s];
  return counts;
}

// AlleleCounts is a 2D array: 1st dimension is locus, 2nd dimension is state*population
// this function returns counts for all populations 
const int *AlleleFreqs::GetAlleleCounts(int locus)const
{
  return( AlleleCounts[locus] );
}
// this function returns counts for the population specified
std::vector<int> AlleleFreqs::GetAlleleCounts( int locus, int population )const
{
  std::vector<int> counts(NumberOfStates[locus]);
  for(int s = 0; s < NumberOfStates[locus]; ++s)
    counts[s] = AlleleCounts[locus][s*Populations + population];
  return counts;
}
std::vector<double> AlleleFreqs::getAlleleFreqsMAP( int locus, int population )const
{
  std::vector<double> A(NumberOfStates[locus]-1);
  for(int i = 0; i < NumberOfStates[locus]-1; ++i)
    A[i] = AlleleFreqsMAP[locus][i*Populations+population];

  return A;
}
/**AlleleFreqs
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
  vector<double> A(NumberOfStates[locus]-1);
  for(int i = 0; i < NumberOfStates[locus]-1; ++i)
    A[i] = Freqs[locus][i*Populations+population];
  return A;
}
const double* const* AlleleFreqs::GetAlleleFreqs()const{
  return Freqs;
}
// get posterior mode of frequency of allele x, given locus and subpopulation
double AlleleFreqs::GetAlleleProbsMAP( int x, int ancestry , int locus)const
{
  double P;
  if( x < (*Loci)(locus)->GetNumberOfAllelesOfLocus(0) - 1 )
    // calculate posterior mode
    P = AlleleFreqsMAP[locus][ x*Populations+ ancestry ];
  else // frequency of last allele is set by subtracting sum of posterior modes of other alleles from 1 
    {
      P = 1;
      for( int j = 0; j < (*Loci)(locus)->GetNumberOfAllelesOfLocus(0) - 1; j++ )
	P -= AlleleFreqsMAP[locus][ j *Populations+ ancestry ];
    }
  return P;
}

// ******************** Output **************************
void AlleleFreqs::OutputAlleleFreqs()
{
  if( IsRandom() ){
    for( int locus = 0; locus < NumberOfCompositeLoci; locus++ ){
      for( int state = 0; state < NumberOfStates[locus]-1; state++ ){
	allelefreqoutput << "\"" << (*Loci)(locus)->GetLabel(0) << "\",";
	for( int pop = 0; pop < Populations; pop++ ){
	  allelefreqoutput <<  Freqs[locus][state*Populations + pop] << ",";
	  // allelefreqoutput << "count " << AlleleCounts[locus][state*Populations + pop] <<"; ";
	}
	allelefreqoutput << endl;
      }
      allelefreqoutput << endl;
    }
  }
  allelefreqoutput << endl;
}

void AlleleFreqs::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  if( IsHistoricAlleleFreq ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta[j] / samples << " ";
    }
  }
  else if(CorrelatedAlleleFreqs){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumEta[0] / samples << " ";
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

void AlleleFreqs::CloseOutputFile(int iterations, const string* const PopulationLabels)
{
  int nrows = 0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    nrows += NumberOfStates[j]-1;
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
    double q_admix,q_parental,f,H_admix, H_parental, H_combined, pbar;
    for( int k = 0; k < Populations; k++ ){
      H_admix = 0;
      H_parental = 0;
      H_combined = 0;
      q_admix=1.0;
      
      for( int i = 0; i < NumberOfStates[locus] - 1; i++ ){
	H_admix += Freqs[locus][ i *Populations+ k ] * Freqs[locus][ i*Populations+ k ];
	H_parental += HistoricAlleleFreqs[locus][ i*Populations+ k ] * HistoricAlleleFreqs[locus][ i*Populations+ k ];
	pbar = 0.5 * ( Freqs[locus][ i *Populations+ k ] + HistoricAlleleFreqs[locus][ i*Populations+ k ] );
	H_combined += pbar * pbar;
	q_admix -= Freqs[locus][i*Populations +k];
      }
      
      H_admix += q_admix * q_admix;
      double sumHistoric = 0.0;
      for( int i = 0; i < NumberOfStates[locus] - 1; i++ ){
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

// ********** log density and derivatives for Adaptive Rejection sampling of PriorAlleleFreqs  ***********

double fMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  //counts 0 are counts for admixed pop, counts1 for historic/unadmixed pop
  double mu = alpha / eta;
  double logprior = 0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1,1) prior
  double f = logprior - 2 * gsl_sf_lngamma( alpha ) - 2 * gsl_sf_lngamma( eta - alpha );

  f += gsl_sf_lngamma( alpha+counts0[pop*2] ) + gsl_sf_lngamma( eta-alpha+counts0[1+pop*2] );//state 1 + state 2, admixed pop
  f += gsl_sf_lngamma( alpha+counts1[pop*2] ) + gsl_sf_lngamma( eta-alpha+counts1[1+pop*2] );//historic pop

  return f;
}

double dfMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x, y1, y2;
  double logprior = 0.1 / alpha - 0.1 / ( eta - alpha );//Beta(1,1) prior
  double f = logprior;
  x = eta - alpha;
  if(alpha < 0)cout<<"\nError in dfMu in allelefreqs.cc - arg mu to ddigam is negative\n"; 
  ddigam( &alpha, &y1 );
  if(x < 0)cout<<"\nError in dfMu in allelefreqs.cc - arg x to ddigam is negative\n"; 
  ddigam( &x, &y2 );
  f += 2 * ( y2 - y1 );

  //admixed pop
  //first state/allele
  x = alpha + counts0[pop*2];
  ddigam( &x, &y2 );
  f += y2;
  //second state/allele
  x = eta - alpha + counts0[1+pop*2];
  ddigam( &x, &y2 );
  f -= y2;

  //unadmixed pop
  x = alpha + counts1[pop*2];
  ddigam( &x, &y2 );
  f += y2;
  x = eta - alpha + counts1[1+pop*2];
  ddigam( &x, &y2 );
  f -= y2;

  return f;
}

double ddfMu( double alpha, const void* const args ) {
  // 
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x, y1, y2;
  double prior = -0.1 / (alpha*alpha) - 0.1 / (( eta - alpha ) * ( eta - alpha ) );
  double f = prior;
  x = eta - alpha;
  trigam( &alpha, &y1 );
  trigam( &x, &y2 );
  f -= 2 * ( y2 + y1 );

  x = alpha + counts0[pop*2];
  trigam( &x, &y2 );
  f += y2;
  x = eta - alpha + counts0[1+pop*2];
  trigam( &x, &y2 );
  f += y2;

  x = alpha + counts1[pop*2];
  trigam( &x, &y2 );
  f += y2;
  x = eta - alpha + counts1[1+pop*2];
  trigam( &x, &y2 );
  f += y2;

  return f;
}

float AlleleFreqs::getEtaRWSamplerAcceptanceRate(int k)const{
  return TuneEtaSampler[k].getExpectedAcceptanceRate();
}
float AlleleFreqs::getEtaRWSamplerStepsize(int k)const{
  return etastep[k];
}

float AlleleFreqs::getAlphaSamplerAcceptanceRate(int i)const{
  return muSampler[i].getAcceptanceRate();
}
float AlleleFreqs::getAlphaSamplerStepsize(int i)const{
  return muSampler[i].getStepsize();
}
float AlleleFreqs::getEtaSamplerAcceptanceRate(int i)const{
  return EtaSampler[i].getAcceptanceRate();
}
float AlleleFreqs::getEtaSamplerStepsize(int i)const{
  return EtaSampler[i].getStepsize();
}
