//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file Individual.cc
/// Implementation of the Individual class.
//=============================================================================

#include "Individual.h"
#include "config.h" // AGGRESSIVE_RANGE_CHECK, USE_GENOTYPE_PARSER
#include "bclib/rand.h"
#include <cmath>


using genepi::RhoType;


#define PR(x) cout << #x << " = " << x << endl;


// ******* Static Member Declarations
bool Individual::Xdata;
unsigned Individual::X_posn;
unsigned int Individual::numChromosomes;
Genome * Individual::Loci;
int Individual::NumHiddenStates;



//--------------------------------------------------------------------------
// Public rho-proposal and psi-proposal methods, overridden from PedBase,
// ignored for individuals.
void Individual::setRho( double /*nv*/ ) {}
void Individual::startRhoProposal () {}
void Individual::acceptRhoProposal() {}
void Individual::rejectRhoProposal() {}

void Individual::startPsiProposal () {}
void Individual::acceptPsiProposal() {}
void Individual::rejectPsiProposal() {}
//--------------------------------------------------------------------------



//******** Constructors **********

Individual::Individual(unsigned number) :
	myNumber	 ( number ) ,
	PossibleHapPairs ( 0      ) ,
	GenotypesMissing ( 0      ) ,
	missingGenotypes ( 0      ) , // allocated later, if needed
	LocusAncestry    ( 0      )
    {
    }



// Individual::Individual(int number, const Options* const options, const InputData* const Data) {
//   missingGenotypes = 0;//allocated later, if needed
//   Initialise(number, options, Data);
// }

void Individual::Initialise(const Options* const options, const InputData* const Data){
  if( options->isRandomMatingModel() && !isHaploid )
    NumGametes = 2;
  else
    NumGametes = 1;

  // Read sex value if present.
  #if USE_GENOTYPE_PARSER
    SexIsFemale = Data->isFemale( myNumber - 1 );
  #else
    SexIsFemale = Data->isFemale( myNumber );
  #endif

  double L = Loci->GetLengthOfGenome();
  double LX = 0.0;
  if(Xdata) LX = Loci->GetLengthOfXchrm();
  // effective length of genome is L + 0.5*LX if there is an X chrm: i.e. if g=1 or sex is female
  if(SexIsFemale) {
    EffectiveL[0] = L + 0.5*LX;
  } else {
    EffectiveL[0] = L;
  }
  EffectiveL[1] = L + 0.5*LX;

  int numCompositeLoci = Loci->GetNumberOfCompositeLoci();

  // vector of possible haplotype pairs - 2 integers per locus if diploid, 1 if haploid
  PossibleHapPairs = new vector<hapPair>[numCompositeLoci];

  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers

  //initialise genotype probs array and array of indicators for genotypes missing at locus

  size_t AncestrySize; // set size of locus ancestry array
  //gametes holds the number of gametes for each chromosome, either 1 or 2
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if(isHaploid || (!SexIsFemale && Loci->isXChromosome(j))){//haploid on this chromosome
      AncestrySize = Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(1);
    }
    else{
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
    }

    LocusAncestry[j] = new int[ AncestrySize ];
    for ( unsigned i = 0; i < AncestrySize; ++i )
	LocusAncestry[j][i] = 0;
  }

  logLikelihood.value = 0.0;
  logLikelihood.ready = false;
  logLikelihood.HMMisOK = false;

}

//********** Destructor **********
Individual::~Individual() {
  if (PossibleHapPairs) delete[] PossibleHapPairs;

  if(GenotypesMissing){
    for(unsigned j = 0; j < numChromosomes; ++j) {
      delete[] GenotypesMissing[j];
    }
    delete[] GenotypesMissing;
  }
  if(LocusAncestry){
    for(unsigned j = 0; j < numChromosomes; ++j) {
      delete[] LocusAncestry[j];
    }
    delete[] LocusAncestry;
  }

  delete[] missingGenotypes;
}

void Individual::SetUniformAdmixtureProps() {
  Theta.setTo(1.0 / NumHiddenStates);
}

void Individual::setOutcome(double* Y){
  Outcome = Y;
}
void Individual::setCovariates(double* X){
  Covariates = X;
}

void Individual::setGenotypesToMissing(){
  for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c)
    for(unsigned j = 0; j < Loci->GetSizeOfChromosome(c); ++j)
      GenotypesMissing[c][j] = true;
}

void Individual::DeleteGenotypes(){
  unsigned noCompositeLoci = Loci->GetNumberOfCompositeLoci();
  for(unsigned j = 0; j < noCompositeLoci; ++j){

    #if 0 // DDF: this should not be necessary as std::vector::clear() will
	  // execute the destructors of each of its elements
	for ( int k = Loci->getNumberOfLoci(j) ; k-- != 0 ; )
	  genotypes[j][k].clear();
    #endif

    genotypes[j].clear();
  }
  genotypes.clear();
}


/// sets static members, including allocation and deletion of static objects for score tests
void Individual::SetStaticMembers( Genome & pLoci, const Options & options ) {
  Loci = &pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  NumHiddenStates = options.getPopulations();
  Xdata = Loci->isX_data();
  X_posn = 9999; //position of the X chromosome in the sequence of chromosomes in the input data
  if(Xdata) {
    X_posn = Loci->GetChrNumOfLocus(Loci->getFirstXLocus());//too clunky, should simplify
  }
}


const RhoType & Individual::getRho() const {
   return _rho;
}

const PsiType & Individual::getPsi() const {
   return psi;
}

/// Sets genome (Loci)
/// Needed for tests which need to substitute only this one static
/// variable.
// void Individual::setGenome(Genome *pLoci)
// {
//   Loci = pLoci;
// }

// void Individual::setPopulations(int p)
// {
//   Populations = p;
// }

void Individual::HMMIsBad(bool loglikisbad) {
  logLikelihood.HMMisOK = false;
  if(loglikisbad)logLikelihood.ready = false;
}

//******************** Accessors ***********************************************************
const double* Individual::getAdmixtureProps(bool isXChrom) const {
  if (isXChrom)
    return Theta.flatXChromosome(psi);
  else
    return Theta.flat();
}

const std::vector<hapPair > &Individual::getPossibleHapPairs( unsigned int locus ) const {
  #if AGGRESSIVE_RANGE_CHECK
    if ( locus >= Loci->GetNumberOfCompositeLoci() )
	throw std::invalid_argument( "Failed assertion: locus out-of-range" );
  #endif
  return PossibleHapPairs[locus];
}

const int* Individual::getSampledHapPair(int locus)const{
  return sampledHapPairs[locus].haps;
}

/**
 *  Get the locus ancestry by the absolute locus index, ignoring the chromosome.
 */
void Individual::GetLocusAncestry(int locus, int Ancestry[2])const{
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  GetLocusAncestry(c, l, Ancestry);
}

void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const{
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if(isHaploid || ((unsigned)chrm == X_posn && !SexIsFemale))
    Ancestry[1] = Ancestry[0];
  else
    Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm] + locus];
}

///returns value of LocusAncestry at a locus for a particular gamete
int Individual::GetLocusAncestry(int chrm, int gamete, int locus)const{
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

///Indicates whether genotype is missing at all simple loci within a composite locus
bool Individual::GenotypeIsMissing(unsigned int locus)const {
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  return GenotypesMissing[c][l];
}
///Indicates whether genotype is missing at a simple locus
//used by HW score test
bool Individual::simpleGenotypeIsMissing(unsigned locus)const{
  if(!missingGenotypes)throw string("missingGenotypes not allocated");
  return missingGenotypes[locus];
}

bool Individual::isHaploidatLocus(unsigned j)const{
  return (bool)(isHaploid || (!SexIsFemale && Loci->isXLocus(j)));
}
bool Individual::isHaploidIndividual()const{
  return isHaploid;
}
//****************** Log-Likelihoods **********************
// public function:
// calls private function to get log-likelihood at current parameter values, and stores it either as loglikelihood.value or as loglikelihood.tempvalue
// store should be false when calculating energy for an annealed run, or when evaluating proposal for global sum-intensities
double Individual::getLogLikelihood( const Options& options, const bool forceUpdate, const bool store) {

  if(!forceUpdate && logLikelihood.ready)
    return logLikelihood.value;
  else{
    logLikelihood.tempvalue = getLogLikelihood(options, Theta, _rho, (forceUpdate || !logLikelihood.HMMisOK));
    if(store) {
      logLikelihood.value = logLikelihood.tempvalue;
      logLikelihood.ready = true;
      logLikelihood.HMMisOK = true; //because forward probs now correspond to current parameter values
    }				    //and call to UpdateHMM has set this to false
    return logLikelihood.tempvalue;
  }
}

// private function: gets log-likelihood at parameter values specified as arguments, but does not update loglikelihoodstruct
double Individual::getLogLikelihood(const Options& options,
                                    const AdmixtureProportions& theta,
				    const RhoType & rho, bool updateHMM) {
  double LogLikelihood = 0.0;
  for( unsigned int j = 0; j < numChromosomes; j++ ) {
    //cout << Loci->isXChromosome(j) << " ";
    if(updateHMM){// force update of forward probs
      UpdateHMMInputs(j, options, theta, rho);
    }
    LogLikelihood += Loci->getChromosome(j)->HMM->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
    if(isnan(LogLikelihood)) {
      throw string("HMM returns log-likelihood as nan (not a number)\n");
    }
  }
  return LogLikelihood;
}

// public function:
// calls private function to get log-likelihood at current parameter values
// only for the X chromosome
double Individual::getLogLikelihoodXChr(const Options& options,
                                        const bool forceUpdate,
                                        const bool /* store */) {

  return getLogLikelihoodXChr(options, Theta, _rho,
                              forceUpdate || !logLikelihood.ready);
}

// private function: gets log-likelihood for the X chromosome only
// at parameter values specified as arguments
double Individual::getLogLikelihoodXChr(const Options& options,
                                        const AdmixtureProportions& theta,
                                        const RhoType& rho, bool updateHMM) {
  gp_assert(Xdata);

  if (updateHMM) // force update of forward probs
    UpdateHMMInputs(X_posn, options, theta, rho);

  Chromosome *XChr = Loci->getChromosome(X_posn);
  double LogLikelihood = XChr->HMM->getLogLikelihood(!isHaploid && SexIsFemale);
  if (isnan(LogLikelihood))
    throw string("HMM returns log-likelihoodXChr as nan (not a number)\n");

  return LogLikelihood;
}

void Individual::storeLogLikelihood(const bool setHMMAsOK) { // to call if a Metropolis proposal is accepted
    logLikelihood.value = logLikelihood.tempvalue;
    logLikelihood.ready = true;
    if(setHMMAsOK) logLikelihood.HMMisOK = true;
}

double Individual::getLogLikelihoodAtPosteriorMeans(const Options& options) {
  // should set allele freqs also to posterior means, and recalculate prob genotypes at these freqs before calling getloglikelihood
  double LogLikelihood = 0.0;
  for( unsigned int j = 0; j < numChromosomes; j++ ) {
    UpdateHMMInputs(j, options, Theta, _rho);
    LogLikelihood += Loci->getChromosome(j)->HMM->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
  }
  return LogLikelihood;
}

//************** Updating (Public) **********************************************************
void Individual::SampleHiddenStates(const Options& options){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome *C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted
    //Update Forward/Backward probs in HMM
    if( !logLikelihood.HMMisOK ) {
      UpdateHMMInputs(j, options, Theta, _rho);
    }
    // sampling locus ancestry can use current values of forward probability vectors alpha in HMM
    C->HMM->SampleHiddenStates(LocusAncestry[j], (!isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)));
  } //end chromosome loop
  logLikelihood.HMMisOK = true;
}

void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts){
  if( !skipMissingGenotypes || !GenotypesMissing[j][jj]) {
    int anc[2];//to store ancestry states
    GetLocusAncestry(j,jj,anc);
    if(PossibleHapPairs[locus].size() > 1 && !annealthermo) {
      // no need to sample if only one possible hap pair or if annealing for thermo integration
      (*Loci)(locus)->SampleHapPair(&(sampledHapPairs[locus]), PossibleHapPairs[locus], anc);
    }
    // now update allelecounts in AlleleFreqs using sampled hap pair
    // UpdateAlleleCounts does nothing if annealthermo and > 2 alleles
    if(UpdateCounts && !GenotypesMissing[j][jj])
      A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, (gametes[j]==2), annealthermo);
  }
}

void Individual::UpdateAlleleCounts(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool annealthermo)const{
  if(!GenotypesMissing[j][jj]){
    int anc[2];//to store ancestry states
    GetLocusAncestry(j,jj,anc);
    A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, (gametes[j]==2), annealthermo);
  }
}

void Individual::SampleMissingOutcomes(bclib::DataMatrix *Outcome, const vector<bclib::Regression*>& R){
  int NumOutcomes = Outcome->nCols();
  // sample missing values of outcome variable
  for( int k = 0; k < NumOutcomes; k++ ){
    if( Outcome->isMissing( getIndex(), k ) ){
      if( R[k]->getRegressionType() == Linear)
	Outcome->set( getIndex(), k, bclib::Rand::gennor( R[k]->getExpectedOutcome(getIndex()), 1 / sqrt( R[k]->getlambda() ) ));
      else{
	if( bclib::Rand::myrand() * R[k]->getExpectedOutcome(getIndex()) < 1 )
	  Outcome->set( getIndex(), k, 1);
	else
	  Outcome->set( getIndex(), k, 0);
      }
    }
  }
}

void Individual::setPsi(const PsiType& _psi) {
  psi = _psi;
}
