/** 
 *   ADMIXMAP
 *   AlleleFreqs.cc 
 *   Class to hold and update allele frequencies, their prior parameters, allele counts and sums. Also holds and updates dispersion
 *   parameter eta and its prior parameters, for a dispersion model. Also computes Fst if required.
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include "AlleleFreqs.h"
#include "DARS.h"
#include "functions.h"

AlleleFreqs::AlleleFreqs(Genome *pLoci){
  eta = 0;
  psi = 0;
  tau = 0; 
  SumEta = 0;
  psi0 = 0.0;
  Populations = 0;
  RandomAlleleFreqs = false;
  IsHistoricAlleleFreq = false;
  AlleleCounts = 0;
  Freqs = 0;
  AlleleFreqsMAP = 0;
  HistoricAlleleFreqs = 0;
  HistoricAlleleCounts = 0;
  PriorAlleleFreqs = 0;
  //SumAlleleFreqs = 0;
  Fst = 0;
  SumFst = 0;
  Loci = pLoci;

#if ETASAMPLER ==1
  TuneEtaSampler = 0;
  w = 0;
  NumberAccepted = 0;
  Number  = 0;
  etastep = 0;
  etastep0 = 0.0;
  SumAcceptanceProb = 0; 
#elif ETASAMPLER ==2
  initialEtaStepsize = 0.03;//need a sensible value for this
  targetEtaAcceptRate = 0.5;//and this 
#endif 
}

AlleleFreqs::~AlleleFreqs(){
  delete[] Freqs;
  delete[] AlleleCounts;//do properly
  delete[] PriorAlleleFreqs;
  delete[] HistoricAlleleCounts;
  //delete[] SumAlleleFreqs;
  delete[] AlleleFreqsMAP;
  delete[] HistoricAlleleFreqs;
  delete[] Fst;
  delete[] SumFst;
  delete[] psi;
  delete[] tau;
  delete[] SumEta;

#if ETASAMPLER ==1
  delete[] NumberAccepted;
  delete[] SumAcceptanceProb;
  delete[] etastep;
  delete[] TuneEtaSampler;
#elif ETASAMPLER ==2
  delete[] logeta;
  delete[] EtaArgs;
#endif
}

void AlleleFreqs::Initialise(AdmixOptions *options, InputData *data, LogWriter *Log){
  LoadAlleleFreqs(options, data);

  Populations = options->getPopulations();
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ) IsHistoricAlleleFreq = true;
  else IsHistoricAlleleFreq = false;

  if(IsRandom() &&  options->getOutputAlleleFreq() ){
    OpenOutputFile(options);
  }

  //set up alleleprobs and hap pair probs
  //NB: HaplotypePairProbs in Individual must be set first
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    (*Loci)(i)->Initialise(Freqs[i]);
  }

  // ** settings for sampling of dispersion parameter **
  if( IsHistoricAlleleFreq ){
    eta = new double[ Populations ];//dispersion parameters
    psi = new double[ Populations ];//gamma prior shape parameter
    tau = new double[ Populations ];//gamma prior rate parameter
    SumEta = new double[ Populations ];//running sums
    
#if ETASAMPLER == 1
    // ** Settings for random walk sampler
    SumAcceptanceProb = new double[ Populations ];
    w = 10;
    etastep0 = 2.0;
    etastep =new double[ Populations ];
    for(int k = 0; k < Populations; ++k)etastep[k] = etastep0;

    Number = 0;
    NumberAccepted =new int[ Populations ];
    TuneEtaSampler = new AdaptiveRandomWalkMH[ Populations ];
    for( int k = 0; k < Populations; k++ )
      TuneEtaSampler[k].SetParameters( w, etastep0, 0.1, 100, 0.44 );
#elif ETASAMPLER == 2
    // ** Settings for Hamiltonian sampler
    logeta = new double[Populations];
    transform(eta, eta+Populations, logeta, xlog);//logeta = log(eta) 
    EtaArgs = new double*[3];// ** need to finish this

    EtaSampler.SetDimensions(Populations, initialEtaStepsize, 20, targetEtaAcceptRate, etaEnergyFunction, etaGradient);
#endif

    // ** set eta priors **
    if( strlen(options->getEtaPriorFilename()) ){
      Log->logmsg(true,"Loading gamma prior parameters for allele frequency dispersion from ");
      Log->logmsg(true,options->getEtaPriorFilename());
      Log->logmsg(true,".\n");
      const DataMatrix& etaprior = data->getEtaPriorMatrix();

      for( int k = 0; k < Populations; k++ ){
	psi[k] = etaprior.get( k, 0 );
	tau[k] = etaprior.get( k, 1 );
      }
    }
    else{//default priors on eta
      // default gamma prior with mean 400, variance 80 000 
      fill(psi, psi + Populations, 2.0);
      fill(tau, tau+Populations, 0.005);
    }

    Log->logmsg(false, "Gamma prior on dispersion parameters with means and variances:\n");
    for( int k = 0; k < Populations; k++ ){
      Log->logmsg(false, data->GetPopLabels()[k]);Log->logmsg(false, ": ");
      Log->logmsg(false, psi[k]/tau[k]);Log->logmsg(false, "  ");Log->logmsg(false, psi[k]/(tau[k]*tau[k]));Log->logmsg(false, "\n");
    }
    Log->logmsg(false, "\n");

    //double maxeta[ Populations ];
    for( int k = 0; k < Populations; k++ ){
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
	PriorAlleleFreqs[j].SetColumn(k, PriorAlleleFreqs[j].GetColumn(k) * eta[k] / PriorAlleleFreqs[j].GetColumn(k).Sum());
	// double sum = 0.0;
	// 	for(int s = 0; s < (*Loci)(j)->GetNumberOfStates(); ++s)sum+= PriorAlleleFreqs[j](s,k);
	// 	for(int s = 0; s < (*Loci)(j)->GetNumberOfStates(); ++s)PriorAlleleFreqs[j](s,k)*= eta[k] / sum;
	//cout<<PriorAlleleFreqs[j].GetColumn(k).Sum()<<eta[k]<<endl;
	//system("pause");
      }
    }
  
    // ** Open output file for eta **
    if ( options->getIndAdmixHierIndicator()){
      if (strlen( options->getEtaOutputFilename() ) ){
	InitializeEtaOutputFile(options, data->GetPopLabels(), Log); 
      }
      else{
	Log->logmsg(true,"No dispparamfile given\n");
	//exit(1);
      }
    }
    // ** open fst output file if specified **
    if( options->getOutputFST() ){
      OpenFSTFile(options,Log);
    }
  }//end if historicallelefreqs

}

void AlleleFreqs::LoadAlleleFreqs(AdmixOptions *options, InputData *data_)
{
  data_->CheckAlleleFreqs(options, Loci->GetNumberOfCompositeLoci(), Loci->GetNumberOfStates());
  int newrow;
  int row = 0;

  Matrix temp2, elements;
  Matrix_d temporary;

  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  Freqs = new double*[NumberOfCompositeLoci];
  AlleleFreqsMAP = new double*[NumberOfCompositeLoci];
  HistoricAlleleFreqs = new double*[NumberOfCompositeLoci];
  AlleleCounts = new int*[NumberOfCompositeLoci];
  HistoricAlleleCounts = new double*[NumberOfCompositeLoci];
  PriorAlleleFreqs = new Matrix_d[NumberOfCompositeLoci];
  //SumAlleleFreqs = new double*[NumberOfCompositeLoci];
  MuProposal = new std::vector<AdaptiveRandomWalkMH>[NumberOfCompositeLoci];

  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    Freqs[i] = 0;
    AlleleFreqsMAP[i] = 0;
    HistoricAlleleFreqs[i] = 0;
    AlleleCounts[i] = 0;
    HistoricAlleleCounts[i] = 0;
  }
  //Fixed AlleleFreqs
  if( strlen( options->getAlleleFreqFilename() ) ){
    temporary = data_->getAlleleFreqMatrix();
    if( options->getTextIndicator() ){
      temporary = temporary.SubMatrix( 1, temporary.GetNumberOfRows() - 1, 1, Populations );
    }

    for( int i = 0; i < NumberOfCompositeLoci; i++ )
      {
	newrow = row + (*Loci)(i)->GetNumberOfStates() - 1;
	InitialiseAlleleFreqs( temporary.Double().SubMatrix( row, newrow - 1, 0, Populations - 1 ), i, Populations);
	row = newrow;
      }
  }
  else if( strlen( options->getHistoricalAlleleFreqFilename() ) || strlen( options->getPriorAlleleFreqFilename() ) ){
    bool Historic;
    //Historic AlleleFreqs
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      temporary = data_->getHistoricalAlleleFreqMatrix();
      Historic = true; 
 
    } else {
      //Prior on AlleleFreqs
      temporary = data_->getPriorAlleleFreqMatrix();
      Historic = false;
    }
    temporary = temporary.SubMatrix( 1, temporary.GetNumberOfRows() - 1, 1, Populations );

    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      newrow = row + (*Loci)(i)->GetNumberOfStates();
      InitialisePriorAlleleFreqs( temporary.SubMatrix( row, newrow - 1, 0, Populations - 1 ), i,options->getFixedAlleleFreqs(), Historic);
      row = newrow;
    }
  }
  //Default Allele Freqs
  else{
    SetDefaultAlleleFreqs( Populations );
  }

}

void AlleleFreqs::InitialiseAlleleFreqs(Matrix_d NewAlleleFreqs, int i, int Pops){
  /**
   * Sets the frequencies of each haplotype at the ith composite locus.
   * obsolescent - maintained for compatibility with old format 
   * in which fixed allele freqs were specified in input data with last allele omitted 
   * 
   * NewAlleleFreqs - a matrix containing allele frequencies, read from input data
   * Rows index alleles or haplotypes, omitting the last allele
   * as its frequency is 1 minus sum of all other frequencies. 
   * Cols index populations.  Thus, for a composite locus with four states
   *   and European and African populations, the matrix might be:
   *
   *             Population
   *
   *            | EUR | AFR |
   *         ---|-----|-----|
   *          0 | 0.5 | 0.2 |
   *   State  1 | 0.1 | 0.2 |
   *          2 | 0.2 | 0.5 |
   */

  int NumberOfStates = (*Loci)(i)->GetNumberOfStates();
  
  // initialize Freqs
  Freqs[i] = new double[(NumberOfStates-1)*Pops];
  for(int j = 0; j < NumberOfStates-1; ++j)
    for(int k = 0; k < Pops; ++k)
      Freqs[i][j*Populations+k] = NewAlleleFreqs(j,k);
  (*Loci)(i)->SetNumberOfPopulations(Pops);
  (*Loci)(i)->SetRandomAlleleFreqs(RandomAlleleFreqs);
  // set size of allele counts matrix at this locus
  AlleleCounts[i] = new int[NumberOfStates*Pops];
  fill(AlleleCounts[i], AlleleCounts[i]+NumberOfStates*Pops, 0);
}

void AlleleFreqs::InitialisePriorAlleleFreqs(Matrix_d New, int i, bool fixed, bool Historic){
  /**
   * Initialises the frequencies of each allele in the ith
   * composite locus, given Dirichlet priors in matrix New.  Allele freqs
   * are set to their prior expectation 
   * If fixed, allele freqs will be fixed at their prior expectations   
   *
   * New - a matrix containing 
   *   parameters for the Dirichlet prior distribution of the allele frequencies. The first dimension is the allele number, 
   *   being in the range of zero to one ??? two less than the number of states
   *   [see GetNumberOfStates()]. The sum of the prior parameters over all alleles in a population 
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
   * If Historic is true, sets "historical allele frequencies", where the model has been specified to allow the 
   * allele freqs in the admixed population 
   * to vary from the historical allele frequencies in the unadmixed ancestral populations that have 
   * been sampled. 
   */
  double sumalpha;
  //int NumberOfStates = (*Loci)(i)->GetNumberOfStates();

  int Pops = New.GetNumberOfCols();
  (*Loci)(i)->SetNumberOfPopulations(Pops);
  // set size of allele freqs array for this locus
  // Freqs array has only NumberOfStates - 1 elements for each population
  Freqs[i] = new double[((*Loci)(i)->GetNumberOfStates()-1)* Pops];
  //SumAlleleFreqs[i] = new double[((*Loci)(i)->GetNumberOfStates() - 1)* Pops];
  // set size of allele counts array
  // allele counts array has NumberOfStates elements for each population 
  AlleleCounts[i] = new int[(*Loci)(i)->GetNumberOfStates() * Pops];
  fill(AlleleCounts[i], AlleleCounts[i]+(*Loci)(i)->GetNumberOfStates()*Pops, 0);

  // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
  // by dividing each prior parameter by the sum of the parameters.     
  for( int j = 0; j < Pops; j++ ){
    sumalpha = ( New.GetColumn(j) ).Sum();
    for( int k = 0; k < (*Loci)(i)->GetNumberOfStates() - 1; k++ )
      Freqs[i][ k*Populations + j ] = ( New( k, j ) ) / sumalpha;
  }
  
  if(Historic){
    HistoricAlleleFreqs[i] = new double[((*Loci)(i)->GetNumberOfStates() - 1)* Pops];
    fill(HistoricAlleleFreqs[i],HistoricAlleleFreqs[i]+ ((*Loci)(i)->GetNumberOfStates() - 1)* Pops, 0.0);
    HistoricAlleleCounts[i] = new double[New.GetNumberOfRows()*New.GetNumberOfCols()];
    for(int row = 0; row < New.GetNumberOfRows(); ++row)
      for(int col = 0; col < New.GetNumberOfCols(); ++col)
	HistoricAlleleCounts[i][row*New.GetNumberOfCols() +col] = New(row, col);
    PriorAlleleFreqs[i] = New + 0.501; // why add 0.501? 
    RandomAlleleFreqs = true;
    Fst = alloc2D_d(NumberOfCompositeLoci, Pops);
    SumFst = alloc2D_d(NumberOfCompositeLoci, Pops);
    // set size of vector MuProposal
    if( (*Loci)(i)->GetNumberOfStates() > 2 ){
      MuProposal[i].resize( Populations );
      for( int k = 0; k < Populations; k++ ){
	MuProposal[i][k].SetParameters( 10, 0.01, 0.001, 0.1, 0.23 );
      }
    }
  }
  else{ // priorallelefreqs model
    if(!fixed){
      PriorAlleleFreqs[i] = New;
      RandomAlleleFreqs = true;
    }
  }
  (*Loci)(i)->SetRandomAlleleFreqs(RandomAlleleFreqs);
}

void AlleleFreqs::SetDefaultAlleleFreqs(int Pops){
  /**
   * Given the number of ancestral populations, sets default values for
   * allele frequencies (in Freqs) and prior allele frequencies (in PriorAlleleFreqs).
   * also creates array AlleleCounts
   * populations - the number of ancestral populations
   */

  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    int NumberOfStates = (*Loci)(i)->GetNumberOfStates();
    // more duplicated code - should do this within method InitializePriorAlleleFreqs
    (*Loci)(i)->SetNumberOfPopulations(Pops);
    PriorAlleleFreqs[i].SetNumberOfElements((*Loci)(i)->GetNumberOfStates(), Pops);
    // reference prior on allele freqs: all elements of parameter vector set to 0.5
    // this is unrealistic for large haplotypes - should set all elements to sum to 1
    PriorAlleleFreqs[i].SetElements(0.5);
    // initialize frequencies as equal for all alleles at locus
    Freqs[i] = new double[(NumberOfStates-1)*Pops];
    for(int j = 0; j < (NumberOfStates-1)*Pops; ++j)
      Freqs[i][j] = 1.0/(*Loci)(i)->GetNumberOfStates();

    AlleleCounts[i] = new int[(*Loci)(i)->GetNumberOfStates() * Pops];
    fill(AlleleCounts[i], AlleleCounts[i]+(*Loci)(i)->GetNumberOfStates()*Pops,0);
    //SumAlleleFreqs[i] = new double[((*Loci)(i)->GetNumberOfStates() - 1)* Pops];
    RandomAlleleFreqs = true;
    (*Loci)(i)->SetRandomAlleleFreqs(RandomAlleleFreqs);
  }
}

// Method samples allele frequency and prior allele frequency
// parameters.
void AlleleFreqs::Update(int iteration,int BurnIn){
  //Reset();

  if( IsRandom() ){
    // Sample for prior frequency parameters mu, using eta, the sum of the frequency parameters for each locus.
    if(IsHistoricAlleleFreq ){
      for( int i = 0; i < NumberOfCompositeLoci; i++ ){
	if( (*Loci)(i)->GetNumberOfStates() == 2 )
	  SamplePriorAlleleFreqs1D( i);
	else
	  SamplePriorAlleleFreqsMultiDim( i);
      }
    }
    
    // Sample allele frequencies conditional on Dirichlet priors 
    // use these frequencies to set AlleleProbs in CompositeLocus
    // then use AlleleProbs to set HapPairProbs in CompositeLocus
    // this is the only point at which SetHapPairProbs is called, apart from when 
    // the composite loci are initialized
    for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      SampleAlleleFreqs(i, 1 );
      (*Loci)(i)->SetAlleleProbs(Freqs[i]);
      //      if( (*Loci)(i)->GetNumberOfLoci() > 1 ) // ************* this condition must be wrong - remarked out
      (*Loci)(i)->SetHapPairProbs();
    }
    
    // Sample for allele frequency dispersion parameters, eta, using
    // Metropolis random-walk.
    if(  IsHistoricAlleleFreq ){
#if ETASAMPLER == 1
      double etanew, LogPostRatio;
      Number++;
      for( int k = 0; k < Populations; k++ ){
	double mineta = 0;
	vector< Vector_d > munew;
	// Sample eta from truncated log-normal distribution.
	do{
	  etanew = exp( gennor( log( eta[k] ), etastep[k] ) );
	}while( etanew > 5000.0 );
	// Prior log-odds ratio (proposal ratio cancels with a part of the prior ratio)   
	LogPostRatio = ( psi[k] - 1 ) * (log(etanew) - log(eta[k]))
	  - tau[k] * ( etanew - eta[k] );
	// Log-likelihood ratio; numerator of integrating constant
	LogPostRatio += 2 * NumberOfCompositeLoci
	  * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta[k] ) );
	for(int j = 0; j < NumberOfCompositeLoci; j++ ){
	  Vector_d mu = GetPriorAlleleFreqs(j,k);
	  //mineta is a lower bound for proposal etanew
	  if( mineta < 0.1 * eta[k] / mu.MinimumElement() )
	    mineta = 0.1 * eta[k] / mu.MinimumElement();
	  
	  //rescale munew so munew sums to etanew
	  munew.push_back( mu * etanew / eta[k] );
	  double *SumLogFreqs = GetStatsForEta(j,k);
	  for( int l = 0; l < (*Loci)(j)->GetNumberOfStates(); l++ ){
	    // Denominator of integrating constant
	    LogPostRatio += 2*(gsl_sf_lngamma( mu(l) ) - gsl_sf_lngamma( munew[j](l) ));
	    // SumLogFreqs = log phi_1 + log phi_2
	    LogPostRatio += (munew[j](l) - mu(l))*SumLogFreqs[l];
	  }
	  delete[] SumLogFreqs;
	}
	
	// Log acceptance probability = Log posterior ratio since the
	// proposal ratio (log-normal) cancels with prior.
	
	// Acceptance test.
	if( log( myrand() ) < LogPostRatio && mineta < etanew ){
	  eta[k] = etanew;
	  //UpdatePriorAlleleFreqs( k, munew );
	  SumAcceptanceProb[k]++;
	  NumberAccepted[k]++;
	}
	
	if( !( Number % w ) ){
	  etastep[k] = TuneEtaSampler[k].UpdateSigma( NumberAccepted[k] );
	  NumberAccepted[k] = 0;
	}
      }
      
      if( !( Number % w ) ){
	Number = 0;
      }
#elif ETASAMPLER == 2
	EtaSampler.Sample(logeta, EtaArgs);//sample new values for eta
	transform(logeta, logeta+Populations, eta, xexp);//eta = exp(logeta)
	if(!((iteration+1) % 10)){
	  EtaSampler.Tune();//tune Hamiltonian sampler every 10 iterations
    }
#endif
      
	    if( iteration > BurnIn )
	      for(int i=0;i<Populations;++i)SumEta[i]+=eta[i];
    }//end isHistoric
    
	if( iteration > BurnIn && IsHistoricAlleleFreq ){
	  UpdateFst();
	}
  }//end isRandom
    //ResetSumAlleleFreqs();
}

/*
  Given a haplotype pair, h, and the ordered ancestry states at a locus,
  updates the counts of alleles observed in each state of ancestry.
  * should use hap pairs stored in Individual object
  */
void AlleleFreqs::UpdateAlleleCounts(int locus, int h[2], int ancestry[2], bool diploid )
{
  AlleleCounts[locus][ h[0]*Populations + ancestry[0] ]++;
  if(diploid)AlleleCounts[locus][ h[1]*Populations + ancestry[1] ]++;
  //if haploid(ie diploid = false), h[0]==h[1]==genotypes[locus] and ancestry[0]==ancestry[1]
  //and we only count once
}
// get posterior mode of frequency of allele x, given locus and subpopulation
double AlleleFreqs::GetAlleleProbsMAP( int x, int ancestry , int locus)
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


/* 
   gets probability of genotypes given ancestry states. 
   Probs is a KxK array in which rows and cols index 
   paternal and maternal locus ancestry states for each population.  
   this function will be redundant if GetGenotypeProbs is fixed to work with a haploid locus. can then call 
   GetGenotypeProbs directly
   Called by UpdateParameters method in Chromosome class and OnePopulationUpdate in Individual
   chibindicator is only to facilitate the Chib algorithm in Individual; instructs CompositeLocus to use HapPairProbsMAP
   instead of HapPairProbs when allelefreqs are not fixed.

*/
//TODO: haploid case
//void AlleleFreqs::GetGenotypeProbsHaploid( double **Probs, int locus, unsigned short **genotype)
//{

//     // lines below should be replaced by a call to GetGenotypeProbs, which should be extended to 
//     // work with a haploid locus
//     //here Probs has a single column 
//     for(int i=0;i<Populations;++i)Probs[i][0]=0.0;
//      if( (*Loci)(locus)->GetNumberOfLoci() == 1 ){
//         for( int pop = 0; pop < Populations; pop++ ){
// 	  Probs[pop][ 0 ] = GetAlleleProbs( genotype[0] - 1, pop , locus);
//         }
//      }
//      else{
//        Vector_i x = (*Loci)(locus)->decodeGenotype(genotype);
//        int xx = (*Loci)(locus)->HapLoopGetDecimal( x );
//         for( int pop = 0; pop < Populations; pop++ ){
// 	  Probs[pop][ 0 ] = GetAlleleProbs( xx - 1, pop , locus);
//         }
//      }
//}

/**
 * Whether the object should remember the results of sampling.
 *
 * flag - integer representing a boolean. Set true (one) to remember
 *   sampled data. Set false (zero) during burn in.
 * i - locus at which to update
 */
void AlleleFreqs::SampleAlleleFreqs(int i, int flag )
{
  // samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
  // and stores result in array Freqs 
  unsigned NumStates = (*Loci)(i)->GetNumberOfStates();
  double temp[NumStates];
  double *freqs = new double[NumStates];
  
  for( int j = 0; j < Populations; j++ ){
    for(unsigned s = 0; s < NumStates; ++s)temp[s] = 
      PriorAlleleFreqs[i](s,j)+AlleleCounts[i][s*Populations +j];
    gendirichlet(NumStates, temp, freqs);
    for(unsigned s = 0; s < NumStates-1; ++s)Freqs[i][s*Populations+j] = freqs[s];

    // sample HistoricAlleleFreqs as a conjugate Dirichlet update with prior specified by PriorAlleleFreqs
    if( IsHistoricAlleleFreq ){
      for(unsigned s = 0; s < NumStates; ++s)temp[s] = PriorAlleleFreqs[i](s,j)+HistoricAlleleCounts[i][s*Populations+j];
      gendirichlet(NumStates, temp, freqs);
      for(unsigned s = 0; s < NumStates-1; ++s)HistoricAlleleFreqs[i][s*Populations+j] = freqs[s];
    }
  }
  
  //if( flag > 0 ){
  //for(int i = 0; i < NumStates-1; ++i)
  //for(int j = 0; j < Populations; ++j)
  //SumAlleleFreqs[i][i*Populations +j] += Freqs[i][i*Populations+j];
  //}
}


/**
 * Sets the sum of allele frequencies over all sampling iterations to
 * zero.
 */
// void AlleleFreqs::ResetSumAlleleFreqs()
// {
// for( int i = 0; i < NumberOfCompositeLoci; i++ )      
//   for(int j = 0; j < ((*Loci)(i)->GetNumberOfStates()-1)*Populations; ++j){
//     SumAlleleFreqs[i][j] = 0.0;
//  }
// }

/**
 * Used when model is specified to allow allele freqs in admixed
 * population to vary from the "historical" allele freqs in the
 * unadmixed population.  Dirichlet distribution for allele freqs at
 * locus with k alleles is specified with (k - 1) frequency parameters
 * (mu) and with a single dispersion parameter (eta) this method
 * samples mu and eta, and updates Dirichlet parameters for the allele
 * frequencies
 */
double *AlleleFreqs::GetStatsForEta( int locus, int population)
{
  // calculates sufficient stats for update of dispersion parameter
  double *stats = new double[ (*Loci)(locus)->GetNumberOfStates() ];
  double sumHistoric = 0.0, sum = 0.0;
  for( int i = 0; i < (*Loci)(locus)->GetNumberOfStates() - 1; i++ ){
    stats[ i ] = log( Freqs[locus][ i*Populations+ population ] ) + log( HistoricAlleleFreqs[locus][ i*Populations + population ] );
    sum +=Freqs[locus][ i*Populations + population ];
    sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + population ];
  }
  stats[ (*Loci)(locus)->GetNumberOfStates() - 1 ] = log( 1 - sum ) + log( 1 - sumHistoric );
  return stats;
}

// void AlleleFreqs::UpdatePriorAlleleFreqs(int j, const vector<Vector_d>& mu)
// {
//   for( int i = 0; i < NumberOfCompositeLoci; i++ ){
//     PriorAlleleFreqs[i].SetColumn( j, mu[i] );
//     //    double sum;
//     //    Vector_d freqs;
//     //    for( int j = 0; j < Populations; j++ ){
//     //       freqs = PriorAlleleFreqs[i].GetColumn(j);
//     //       sum = freqs.Sum();
//     //       PriorAlleleFreqs[i].SetColumn( j, freqs * eta[j] / sum );
//     //    }
//   }
// }

void AlleleFreqs::SamplePriorAlleleFreqsMultiDim( int locus)
  // problem here is to sample the Dirichlet parameters of a multinomial-Dirichlet distribution
  // for a multi-allelic locus, we sample the Dirichlet proportion parameters conditional on the  
  // dispersion parameter and the allele counts in admixed and historic populations
  // by a Metropolis random walk
{
  vector<int> accept(Populations,0);
  Vector_d mu1, mu2; // mu1 is current vector of proportion parameters, mu2 is proposal
  for( int j = 0; j < Populations; j++ ){
    double Proposal1=0, Proposal2=0, f1=0, f2=0;
    mu1 = PriorAlleleFreqs[locus].GetColumn(j) / eta[j];
    mu2 = gendirichlet( mu1 / MuProposal[locus][j].GetSigma() );
        
    for( int i = 0; i < (*Loci)(locus)->GetNumberOfStates(); i++ ){
      // priors on proportion parameters are apparently Dirichlet(0.1, ,,, 0,1) 
      f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) ); 
      f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );
      Proposal1 += (eta[j] * mu2(i) - 1) * log( mu1(i) ) - gsl_sf_lngamma( eta[j] * mu2(i) );
      Proposal2 += (eta[j] * mu1(i) - 1) * log( mu2(i) ) - gsl_sf_lngamma( eta[j] * mu1(i) );
    }
      
    int numberofstates = mu1.GetNumberOfElements();
    for( int k = 0; k < numberofstates; k++ ){
      f1 -= Populations * gsl_sf_lngamma( mu1( k )* eta[j] );
      f2 -= Populations * gsl_sf_lngamma( mu2( k )* eta[j] );
      f1 += gsl_sf_lngamma( mu1( k )* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f2 += gsl_sf_lngamma( mu2( k )* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f1 += gsl_sf_lngamma( mu1( k )* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
      f2 += gsl_sf_lngamma( mu2( k )* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
    }
    if( log(myrand()) < f2 - f1 - Proposal2 + Proposal1 ){
      PriorAlleleFreqs[locus].SetColumn( j, mu2 * eta[j] );
      accept[j] = 1;
      MuProposal[locus][j].Event(true);
    }
    else
      MuProposal[locus][j].Event(false);
  }
}

void AlleleFreqs::SamplePriorAlleleFreqs1D( int locus)
  // with a dispersion model, we sample PriorAlleleFreqs conditional on the observed counts 
  // (with the realized allele freqs integrated out) from a distribution that is 
  // proportional to the product of two binomial-beta likelihoods (for a diallelic locus)
  // we sample the proportion parameter of the beta distribution
  // conditional on the realized allele counts in the admixed population, 
  // the allele counts in the historic population, and the dispersion parameter
  // using an adaptive rejection sampler
  //Note: here NumberOfStates == 2  
{
  double lefttruncation = 0.1;
  double MuParameters[2];
  int counts0[2 * Populations];
  double counts1[2 * Populations];

  // Construct adaptive rejection sampler for mu.
  for(int i = 0; i < 2; ++i)
    for(int j = 0; j < Populations; ++j) {
      counts0[i + j*2] = AlleleCounts[locus][i*Populations +j];
      counts1[i+ j*2] = HistoricAlleleCounts[locus][i*Populations +j];
    }
  //warning message - move to wherever loci data are read
  //    for(int row=0;row<counts0(0).GetNumberOfRows();++row)for(int col=0;col<counts0(0).GetNumberOfCols();++col)
  //      if(counts0(0)(row,col)==0 && counts1(0)(row,col)==0){
  //        cout<<"Warning: zero copies of allele ("<<row<<","<<col<<") in both admixed and unadmixed samples"<<endl;
  //    //poss return name of comp locus
  //    }
  DARS SampleMu( 0, 0, 0, MuParameters, fMu, dfMu, ddfMu, counts0, counts1 );

  SampleMu.SetLeftTruncation( lefttruncation );
  for( int j = 0; j < Populations; j++ ){
    //      sum = 0.0;
    MuParameters[0] = eta[j];
    MuParameters[1] = j;
    //       MuParameters(1) = (float)AlleleCounts( NumberOfStates - 1, j );
    //       MuParameters(2) = (float)HistoricAlleleCounts[ (NumberOfStates - 1)*Populations + j ];
    //       MuParameters(3) = (float)AlleleCounts( k, j );
    //       MuParameters(4) = (float)HistoricAlleleCounts[ k*Populations + j ];
    SampleMu.SetRightTruncation( eta[j] - lefttruncation );
    SampleMu.UpdateParameters( MuParameters);
    PriorAlleleFreqs[locus]( 0, j ) = SampleMu.Sample();
    // Last prior frequency parameter is determined by; sum of mu's = eta.
    PriorAlleleFreqs[locus]( 1, j ) = eta[j] - PriorAlleleFreqs[locus]( 0, j );
  }
}


void AlleleFreqs::InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log)
{
  outputstream.open( options->getEtaOutputFilename(), ios::out );
  if( !outputstream )
    {
      Log->logmsg(true,"ERROR: Couldn't open dispparamfile\n");
      exit( 1 );
    }
  else{
    Log->logmsg(true,"Writing dispersion parameters to ");
    Log->logmsg(true,options->getEtaOutputFilename());
    Log->logmsg(true,"\n");
    if( options->getTextIndicator()  && options->getAnalysisTypeIndicator() >= 0)
      {
	//Dispersion parameters (eta)
	if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	  for( int k = 0; k < Populations; k++ ){
	    outputstream << "\"eta." << PopulationLabels[k].substr(1);
	  }
	}
	outputstream << endl;
      }
  }
}

void AlleleFreqs::OutputErgodicAvg( int samples,AdmixOptions *options, std::ofstream *avgstream)
{
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta[j] / samples << " ";
    }
  }
}

void AlleleFreqs::OutputEta(int iteration, AdmixOptions *options, LogWriter *Log){
  if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
    //output to logfile
    if( !options->useCOUT() || iteration == 0 )
      {
	for( int j = 0; j < Populations; j++ ){
	  Log->width(9);
	  Log->write(eta[j],6);
	}
      }
    //output to screen
    if( options->useCOUT() )
      {
	for( int j = 0; j < Populations; j++ ){
	  cout.width(9);
	  cout << setprecision(6) << eta[j] << " ";
	}
      }
    //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      for( int j = 0; j < Populations; j++ ){
	outputstream.width(9);
	outputstream << setprecision(6) << eta[j];
      }
      outputstream << endl;
    }
  }
}

Genome *AlleleFreqs::getLoci(){
  return Loci;
}
CompositeLocus *AlleleFreqs::getLocus(int i){
  return (CompositeLocus *)((*Loci)(i));
}
//any call to this function should be replaced by a call to the same function in Genome
int AlleleFreqs::GetNumberOfCompositeLoci(){
  return NumberOfCompositeLoci;
}

void AlleleFreqs::OutputFST(bool IsPedFile){ // should remove all refs to IsPedFile
  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    if(IsPedFile)
      fstoutputstream << "\"" << (*Loci)(j)->GetLabel(0) << "\"";
    else
      fstoutputstream << (*Loci)(j)->GetLabel(0);
    for(int k=0; k<Populations; ++k)fstoutputstream << " " << Fst[j][k];
    fstoutputstream << endl;
  }
}

void AlleleFreqs::ResetAlleleCounts(){
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    fill(AlleleCounts[i], AlleleCounts[i]+(*Loci)(i)->GetNumberOfStates()*Populations, 0);
  }
}
void AlleleFreqs::OpenFSTFile(AdmixOptions *options,LogWriter *Log){
  Log->logmsg(true, "Writing ergodic averages of FSTs to: ");
  Log->logmsg(true,options->getFSTOutputFilename());
  Log->logmsg(true,"\n");
  fstoutputstream.open( options->getFSTOutputFilename(), ios::out );
  if( !fstoutputstream ){
    Log->logmsg(true,"ERROR: Couldn't open fstoutputfile\n");
    exit( 1 );
  }
}

// should set allele freqs to the posterior mode
// but apparently just sets to current value
// probably doesn't matter where there is strong prior on allele freqs, and used in Chib algorithm 
// which just requires a value near the posterior mode
// set AlleleFreqsMAP and getAlleleFreqs are called by Individual object
void AlleleFreqs::setAlleleFreqsMAP()
{
  int NumberOfStates;
  for(int i = 0; i < NumberOfCompositeLoci;++i){
    NumberOfStates = (*Loci)(i)->GetNumberOfStates();
    if(AlleleFreqsMAP[i]==0)
      AlleleFreqsMAP[i] = new double[(NumberOfStates-1)*Populations];
    for(int j = 0; j < NumberOfStates - 1; ++j)for(int k = 0; k < Populations; ++k)
      AlleleFreqsMAP[i][j*Populations+k] = Freqs[i][j*Populations+k];
  }
}

// Indicates whether allele frequencies are fixed or random.
// returns a boolean true if allele frequencies are random, false otherwise.
bool AlleleFreqs::IsRandom()
{
  return( RandomAlleleFreqs );
}

/**PriorAlleleFreqs
 * Returns Dirichlet parameters for allele frequencies for a particular population and locus.
 * 
 * population - the number of the population (zero based)
 * locus - the number of the locus
 * returns:
 * a vector_d containing Dirichlet parameters for frequencies of each allele at the locus. 
 * Expected frequencies are calculated by dividing each parameter by the sum of parameters
 */
Vector_d AlleleFreqs::GetPriorAlleleFreqs( int locus, int population )
{
  return( PriorAlleleFreqs[locus].GetColumn( population ) );
}

// AlleleCounts is a 2D array: 1st dimension is locus, 2nd dimension is state*population
// this function returns counts for all populations 
int *AlleleFreqs::GetAlleleCounts(int locus)
{
  return( AlleleCounts[locus] );
}
// this function returns counts for the population specified
Vector_i AlleleFreqs::GetAlleleCounts( int locus, int population )
{
  Vector_i counts((*Loci)(locus)->GetNumberOfStates());
  for(int s = 0; s < (*Loci)(locus)->GetNumberOfStates(); ++s)
    counts(s) = AlleleCounts[locus][s*Populations + population];
  return counts;
}
Vector_d AlleleFreqs::getAlleleFreqsMAP( int locus, int population )
{
  Vector_d A((*Loci)(locus)->GetNumberOfStates()-1);
  for(int i = 0; i < (*Loci)(locus)->GetNumberOfStates()-1; ++i)
    A(i) = AlleleFreqsMAP[locus][i*Populations+population];
  //return( AlleleFreqsMAP[locus].GetColumn( population ) );
  return A;
}
/**AlleleFreqs
 * Gets the frequencies of each haplotype in a composite locus.
 *
 * returns:
 * array containing the frequencies of each allele
 * one row per allele state, one column per population
 */
double *AlleleFreqs::GetAlleleFreqs(int locus)
{
  return( Freqs[locus] );
}
Vector_d AlleleFreqs::GetAlleleFreqs( int locus, int population )
{
  Vector_d A((*Loci)(locus)->GetNumberOfStates()-1);
  for(int i = 0; i < (*Loci)(locus)->GetNumberOfStates()-1; ++i)
    A(i) = Freqs[locus][i*Populations+population];
  return A;
}
double **AlleleFreqs::GetAlleleFreqs(){
  return Freqs;
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
      
      for( int i = 0; i < (*Loci)(locus)->GetNumberOfStates() - 1; i++ ){
	H_admix += Freqs[locus][ i *Populations+ k ] * Freqs[locus][ i*Populations+ k ];
	H_parental += HistoricAlleleFreqs[locus][ i*Populations+ k ] * HistoricAlleleFreqs[locus][ i*Populations+ k ];
	pbar = 0.5 * ( Freqs[locus][ i *Populations+ k ] + HistoricAlleleFreqs[locus][ i*Populations+ k ] );
	H_combined += pbar * pbar;
	q_admix -= Freqs[locus][i*Populations +k];
      }
      
      H_admix += q_admix * q_admix;
      double sumHistoric = 0.0;
      for( int i = 0; i < (*Loci)(locus)->GetNumberOfStates() - 1; i++ ){
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

void AlleleFreqs::OpenOutputFile(AdmixOptions *options)
{
  allelefreqoutput.open(options->getAlleleFreqOutputFilename(), ios::out );
  if( !allelefreqoutput && options->getAnalysisTypeIndicator() >= 0){
    cerr << "Warning: Couldn't open allelefreqsoutputfile: " << options->getAlleleFreqOutputFilename() << endl;
    //exit( 1 );
  }
  else
    allelefreqoutput << "structure(.Data=c(" << endl;
}

void AlleleFreqs::OutputAlleleFreqs()
{
  if( IsRandom() ){
    for( int locus = 0; locus < NumberOfCompositeLoci; locus++ ){
      for( int state = 0; state < (*Loci)(locus)->GetNumberOfStates()-1; state++ ){
	allelefreqoutput << getLocus(locus)->GetLabel(0) << ",";
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

void AlleleFreqs::CloseOutputFile(int iterations, string* PopulationLabels)
{
  int nrows = 0;
  for(int j = 0; j < NumberOfCompositeLoci; ++j)
    nrows += (*Loci)(j)->GetNumberOfStates()-1;
  allelefreqoutput << ")," << endl;
  allelefreqoutput << ".Dim = c(";
  allelefreqoutput << Populations+1 << ",";
  allelefreqoutput << nrows << ",";
  allelefreqoutput << iterations;
  allelefreqoutput << ")," << endl;
  allelefreqoutput << ".Dimnames=list(c(\"Locus\",";
  for (int i = 0; i < Populations; i++){
    allelefreqoutput << PopulationLabels[i];
    if(i < Populations-1){
      allelefreqoutput << ",";
    }
  }
  allelefreqoutput << "), character(0), character(0)))" << endl;
  allelefreqoutput.close();
}

double fMu( const double* parameters, const int *counts0,  const double *counts1, double mu )
{
  //    Vector_d rn(2), ri(2);
  //    rn(0) = parameters(1);
  //    rn(1) = parameters(2);
  //    ri(0) = parameters(3);
  //    ri(1) = parameters(4);
  int pop = (int)parameters[1];
  double eta = parameters[0];
  double prior = 0.1 * log( mu / eta ) + 0.1 * log( 1 - mu / eta );
  double f = prior - 2 * gsl_sf_lngamma( mu ) - 2 * gsl_sf_lngamma( eta - mu );
  //   for( int i = 0; i < 2; i++ )
  //      f += gsl_sf_lngamma( mu + ri(i) ) + gsl_sf_lngamma( eta - mu + rn(i) );
  f += gsl_sf_lngamma( mu+counts0[pop*2] ) + gsl_sf_lngamma( eta-mu+counts0[1+pop*2] );
  f += gsl_sf_lngamma( mu+counts1[pop*2] ) + gsl_sf_lngamma( eta-mu+counts1[1+pop*2] );

  return f;
}

double dfMu( const double* parameters, const int *counts0, const double *counts1, double mu )
{
  //    Vector_d rn(2), ri(2);
  //    rn(0) = parameters(1);
  //    rn(1) = parameters(2);
  //    ri(0) = parameters(3);
  //    ri(1) = parameters(4);
  int pop = (int)parameters[1];
  double eta = parameters[0], x, y1, y2;
  double prior = 0.1 / mu - 0.1 / ( eta - mu );
  double f = prior;
  x = parameters[0] - mu;
  if(mu < 0)cout<<"\nError in dfMu in compositelocus.cc - arg mu to ddigam is negative\n"; 
  ddigam( &mu, &y1 );
  if(x < 0)cout<<"\nError in dfMu in compositelocus.cc - arg x to ddigam is negative\n"; 
  ddigam( &x, &y2 );
  f += 2 * ( y2 - y1 );

  //    for( int i = 0; i < 2; i++ ){
  //       x = mu + ri(i);
  //       ddigam( &x, &y2 );
  //       f += y2;
  //       x = eta - mu + rn(i);
  //       ddigam( &x, &y2 );
  //       f -= y2;
  //    }

  x = mu + counts0[pop*2];
  ddigam( &x, &y2 );
  f += y2;
  x = eta - mu + counts0[1+pop*2];
  ddigam( &x, &y2 );
  f -= y2;

  x = mu + counts1[pop*2];
  ddigam( &x, &y2 );
  f += y2;
  x = eta - mu + counts1[1+pop*2];
  ddigam( &x, &y2 );
  f -= y2;

  return f;
}

double ddfMu( const double* parameters, const int *counts0, const double *counts1, double mu )
{
  //    Vector_d rn(2), ri(2);
  //    rn(0) = parameters(1);
  //    rn(1) = parameters(2);
  //    ri(0) = parameters(3);
  //    ri(1) = parameters(4);
  int pop = (int)parameters[1];
  double eta = parameters[0], x, y1, y2;
  double prior = -0.1 / (mu*mu) - 0.1 / (( eta - mu ) * ( eta - mu ) );
  double f = prior;
  x = parameters[0] - mu;
  trigam( &mu, &y1 );
  trigam( &x, &y2 );
  f -= 2 * ( y2 + y1 );

  //    for( int i = 0; i < 2; i++ ){
  //       x = mu + ri(i);
  //       trigam( &x, &y2 );
  //       f += y2;
  //       x = eta - mu + rn(i);
  //       trigam( &x, &y2 );
  //       f += y2;
  //    }

  x = mu + counts0[pop*2];
  trigam( &x, &y2 );
  f += y2;
  x = eta - mu + counts0[1+pop*2];
  trigam( &x, &y2 );
  f += y2;

  x = mu + counts1[pop*2];
  trigam( &x, &y2 );
  f += y2;
  x = eta - mu + counts1[1+pop*2];
  trigam( &x, &y2 );
  f += y2;

  return f;
}
#if ETASAMPLER == 2
/*
  Energy function (-log density) and gradient function for dispersion parameters eta, used in Hamiltonian sampler
  args[0] = psi, shape parameters of gamma prior
  args[1] = tau, rate parameters  "   "     "

*/
double AlleleFreqs::etaEnergyFunction(unsigned dim, const double* const theta, const double* const*args){
  //TODO
}
void AlleleFreqs::etaGradient(unsigned dim,const double* const theta, const double* const* args, double *g){
  //TODO
}
#endif
