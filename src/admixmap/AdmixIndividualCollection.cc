//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AdmixIndividualCollection.cc
/// Implementation of the AdmixIndividualCollection class.
//=============================================================================

#include "AdmixIndividualCollection.h"
#include "InputAdmixData.h"
#include "bclib/LogWriter.h"
#include "bclib/Regression.h"
#include <iomanip>
#include <iostream>


#define DEBUG_DTYPE	0 // DDF: debug which objects in the collection are of which class

#if DEBUG_DTYPE
    #include <cxxabi.h>
#endif

using namespace std;

using genepi::RhoType;
using genepi::Pedigree;



// **** CONSTRUCTORS  ****
// AdmixIndividualCollection::AdmixIndividualCollection() {
//     SetNullValues();
// }
void AdmixIndividualCollection::SetNullValues(){
  SUPER::SetNullValues();
  TestInd = 0;
  sizeTestInd = 0;
  indadmixoutput = 0;
  SumLogTheta = 0;
}

AdmixIndividualCollection::AdmixIndividualCollection( const AdmixOptions &   options ,
						      const InputAdmixData & data    ,
						      Genome &		     loci    ) :
    SUPER( data.isPedFile() ? data.getPeds().size() : data.getNumberOfIndividuals() ,
	   options.getPopulations() ,
	   loci.GetNumberOfCompositeLoci() )
  {
  SetNullValues();

  size = NumInd;

  AdmixedIndividual::SetStaticMembers( loci, options );

  unsigned i0 = 0; // used to offset numbering of other individuals (not the one under test)

  // Fill separate individuals.
  if ( options.getTestOneIndivIndicator() ) {
    sizeTestInd = options.getNumAnnealedRuns()+1;
    // Create array of copies of TestInd
    TestInd = new PedBase*[ sizeTestInd ];
    SumEnergy = new double[ sizeTestInd ];
    fill( SumEnergy, SumEnergy+sizeTestInd, 0.0 );
    SumEnergySq = new double[ sizeTestInd ];
    fill( SumEnergySq, SumEnergySq+sizeTestInd, 0.0 );

    for ( int i = 0; i < sizeTestInd; ++i ) {
      TestInd[i] = new AdmixedIndividual(1, "", &options, &data, true);
    }
    ++i0;
    --size;
    if ( ! TestInd[0]->isHaploidIndividual() )
	++NumDiploidIndividuals;
  }

  _child = new PedBase*[ size ];

  if ( data.isPedFile() || options.getUsePedForInd() )
      {
      bool hasPedigrees  = false;
      bool hasUnrelateds = false;

      for ( unsigned int i = 0; i < size; i++ )
	{
	const Pedigree & ped = data.getPed( i + i0 );
	if ( options.getUsePedForInd() || (ped.getNMembers() > 1) )
          {
            hasPedigrees = true;
	    _child[i] = const_cast<genepi::Pedigree*>( &ped );
          }
	else
          {
            hasUnrelateds = true;
	    _child[i] = new AdmixedIndividual( i+i0+1, data.getFamId(i + i0),
					       &options, &data, false );
          }
	++NumDiploidIndividuals;
	}

      if ( hasPedigrees && hasUnrelateds && !options.getUsePedForInd() )
	throw("Cannot mix pedigrees and unrelated individuals with "
	      "use-pedigree-for-individual=0.");
      }
  else
      {
      for (unsigned int i = 0; i < size; i++ )
	{
	_child[i] = new AdmixedIndividual( i+i0+1, data.getFamId(i + i0),
					   &options, &data, false);
	if ( ! _child[i]->isHaploidIndividual() )
	  ++NumDiploidIndividuals;
	}
      }
}


// ************** DESTRUCTOR **************
AdmixIndividualCollection::~AdmixIndividualCollection() {
  //cout << "\n Deleting individual objects\n" << flush;
  //  AdmixedIndividual::DeleteStaticMembers();

  if ( TestInd )
    {
    for ( int i = 0 ; i < sizeTestInd ; ++i )
	delete TestInd[i];
    }

  delete[] TestInd;
  delete indadmixoutput;
  delete[] SumLogTheta;

  #if DEBUG_DTYPE
    for ( size_t i = 0 ; i < size ; i++ )
      {
      int status;
      char * const realname = abi::__cxa_demangle( typeid(*(_child[i])).name(), 0, 0, &status );
      const char * cl = (status == 0) ? realname : typeid(*(_child[i])).name();
      fprintf( stderr, "Class of %zu: %s\n", i, cl );
      if ( realname != 0 )
	  free( realname );
      }
  #endif

  for(unsigned int i = 0; i < size; i++)
      {
      PedBase * const pb = _child[i];
      if ( dynamic_cast<AdmixedIndividual *>( pb ) != 0 )
	delete pb;
      }

  delete[] _child;
}


// ************** INITIALISATION AND LOADING OF DATA **************

void AdmixIndividualCollection::Initialise(const AdmixOptions& options, const Genome & Loci,
				      const Vector_s& PopulationLabels, bclib::LogWriter &Log){
  Log.setDisplayMode(bclib::Quiet);
  //Open indadmixture file
  if ( strlen( options.getIndAdmixtureFilename() ) ){
    Log << "Writing individual-level parameters to "
        << options.getIndAdmixtureFilename() << "\n";
    indadmixoutput = new IndAdmixOutputter(options, Loci, PopulationLabels);
  }
  else {
    Log << "No indadmixturefile given\n";
  }

  //Set locusfortest if specified
  if( options.getLocusForTestIndicator() )
    _locusfortest = Loci.GetChrmAndLocus( options.getLocusForTest() );

  // allocate array of sufficient statistics for update of population admixture parameters
  SumLogTheta = new double[ options.getPopulations() ];

  //allocate array of sufficient statistics for update of locus-specific sumintensities
}

void AdmixIndividualCollection::DrawInitialAdmixture( const PedBase::AlphaType & alpha ){
    //draw initial values for individual admixture proportions
  for(unsigned int i = 0; i < size; i++)
    getElement(i).drawInitialAdmixtureProps(alpha);
    if(TestInd)
      for(int i = 0; i < sizeTestInd; i++)
	TestInd[i]->drawInitialAdmixtureProps(alpha);

}

void AdmixIndividualCollection::resetStepSizeApproximators(int k) {
  for(unsigned i = 0; i < size; i++)
    getElement(i).resetStepSizeApproximator(k);
}

void AdmixIndividualCollection::LoadData( const AdmixOptions & options, const InputAdmixData & data_ ) {

  SUPER::LoadData( options, data_,
			(!options.getTestForAdmixtureAssociation() && options.getPopulations() > 1));
  if ( strlen( options.getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry( &data_ );
  }
}

void AdmixIndividualCollection::LoadRepAncestry(const InputAdmixData* const data_){
  ReportedAncestry = new bclib::DataMatrix[NumInd];
  bclib::DataMatrix& temporary = (bclib::DataMatrix&)data_->getReportedAncestryMatrix();
  for( unsigned i = 0; i < temporary.nRows() / 2; i++ )
    ReportedAncestry[i] = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.nCols() - 1 );

}

// ************** UPDATING **************
//TODO: ?? have next two functions use those in base class and have ones here only operate on test individuals
void AdmixIndividualCollection::setGenotypeProbs(const Genome* const Loci){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  const unsigned *sizeOfChromosome = Loci->GetSizesOfChromosomes();

  for (unsigned j = 0; j < nchr; ++j) {

    for (unsigned int i = 0; i < size; ++i)
      getElement(i).SetGenotypeProbs(j, sizeOfChromosome[j], false);

    if(TestInd)
      for(int i = 0; i < sizeTestInd; ++i)
        TestInd[i]->SetGenotypeProbs(j, sizeOfChromosome[j], false);
  }
}

void AdmixIndividualCollection::annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses){
  for(unsigned j = 0; j < nchr; ++j){
    if(TestInd) { // anneal test individual only
      for(int i = 0; i < sizeTestInd; ++i)
	TestInd[i]->AnnealGenotypeProbs(j, Coolnesses[i]);

    } else { // anneal all individuals
      for(unsigned int i = 0; i < size; i++) {
	getElement(i).AnnealGenotypeProbs(j, coolness);
      }
    }
  }
}

/// HMMUpdates() performs the following actions:
/// -# Samples admixture with random walk, on even-numbered iterations
/// -# Samples Locus Ancestry (after updating HMM)
/// -# Samples Jump Indicators and accumulates sums of (numbers of arrivals)
///    and (ancestry states where there is an arrival)
/// -# updates score, info and score squared for ancestry score tests
///
/// Coolness is not passed as argument to this function because annealing has
/// already been implemented by calling annealGenotypeProbs().
void AdmixIndividualCollection::HMMUpdates(int iteration, const AdmixOptions & options,
						const vector<bclib::Regression*> &R,
						const PopAdmix::PopThetaType & poptheta,
						const PedBase::AlphaType & alpha,
						AffectedsOnlyTest& affectedsOnlyTest,
						CopyNumberAssocTest& ancestryAssocTest, bool anneal){
  vector<double> lambda;
  vector<const double*> beta;
  double dispersion = 0.0;

  const bool even_numbered_iteration = ((iteration & 1) == 0);
  const bool _anneal = (anneal && !options.getTestOneIndivIndicator());
  const bool updateScoreTests = iteration > options.getBurnIn()
    && (options.getTestForAffectedsOnly() || options.getTestForLinkageWithAncestry());
  const bool updateSumLogTheta = !_anneal && iteration > options.getBurnIn();

  if(even_numbered_iteration){//lambda, beta and dispersion are only required for random-walk update of admixture
    for(int i = 0; i < options.getNumberOfOutcomes(); ++i){
      lambda.push_back( R[i]->getlambda());
      beta.push_back( R[i]->getbeta());
    }
    if(R.size()>0) dispersion = R[0]->getDispersion();
  }

  // skip test individuals when obtaining derivative inverse-link
  // (test indivs are not included in regression)
  int i0 = 0;
  if (options.getTestOneIndivIndicator())
    i0 = 1;

  //reset sufficient stats for update of pop admixture params to 0
  //  if((iteration %2))
  fill(SumLogTheta, SumLogTheta+options.getPopulations(), 0.0);

  // reset arrays used in score test to 0. This must be done here as the
  // B matrix is updated after sampling admixture
  if(iteration >= options.getBurnIn()){
    ancestryAssocTest.Reset();
    affectedsOnlyTest.Reset();
  }

  // Now loop over individuals
  const int ssize = int(size);
  #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
    #pragma omp parallel for default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
  #endif
  for ( int i = 0 ; i < ssize ; i++ ){

    PedBase & el = getElement(i);

    // ** set SumLocusAncestry and SumNumArrivals to 0
    el.ResetSufficientStats();

    if(Populations >1){//no updates required for single-population model
      // ** update theta with random-walk proposal on even-numbered iterations
      if(even_numbered_iteration){
        double DinvLink = 1.0;
	// Accessing global data here, must protect -- let's find a better way?
	#if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP && SAMPLE_THETA_CALL_CRITICAL
	    #pragma omp critical
	#endif
	    { // BEGIN SCOPE BLOCK
	    if(R.size())DinvLink = R[0]->DerivativeInverseLinkFunction(i+i0);
	    el.SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates,
                                     &Covariates, beta, poptheta, options, alpha, DinvLink,
                                     dispersion, ancestryAssocTest, true, updateSumLogTheta);
	    } // END SCOPE BLOCK
      }


      // Since these methods are not implemented for pedigrees, we explicitly
      // exclude them.
      if ( ! el.isPedigree() )
	{

	// ** Run HMM forward recursions, if required, and sample locus ancestry
	el.SampleHiddenStates( options );

	// ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
	el.SampleJumpIndicators( ! options.isGlobalRho() );

	}

      // ** Update score, info and varscore for ancestry score tests.
      // NOTE: The access to the quasi-global test-objects (Outcome, Covariates,
      //    affectedsOnlyTest, etc.) here raises concurrency issues, currently
      //    addressed by marking certain sections as critical _inside_
      //    Pedigree::UpdateScores().
      if ( updateScoreTests )
        el.UpdateScores(options, &Outcome, &Covariates, R, affectedsOnlyTest, ancestryAssocTest);
    }

  }
}

///samples individual-level sumintensities and admixture
void AdmixIndividualCollection::SampleParameters(int iteration, const AdmixOptions& options,
						const vector<bclib::Regression*> &R, const PopAdmix::PopThetaType & poptheta,
						const PedBase::AlphaType &alpha, double rhoalpha, double rhobeta,
						CopyNumberAssocTest& ancestryAssocTest, bool anneal=false){
  //sufficient statistics have been stored in Individuals
  // coolness is not passed as argument to this function because annealing has already been implemented by
  // calling annealGenotypeProbs
  vector<double> lambda;
  vector<const double*> beta;
  double dispersion = 0.0;
  const bool updateSumLogs = !anneal && iteration > options.getBurnIn();
  const bool isGlobalRho   = options.isGlobalRho();

  for(unsigned i = 0; i < R.size(); ++i){
    lambda.push_back( R[i]->getlambda());
    beta.push_back( R[i]->getbeta());
  }
  if (R.size() > 0)
    dispersion = R[0]->getDispersion();

  // Test Individuals: all vars updated here as they contribute nothing to
  // score tests or update of allele freqs
  int i0 = 0;
  if(options.getTestOneIndivIndicator()) {// anneal likelihood for test individual only
    i0 = 1;
    for(int i = 0; i < sizeTestInd; ++i){
      // ** set SumLocusAncestry and SumNumArrivals to 0
      TestInd[i]->ResetSufficientStats();

      if (Populations > 1) {

        // ** update theta with random-walk proposal on even-numbered iterations
        if ( !(iteration % 2) )
          TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType,
                                  lambda, NumCovariates, &Covariates,
                                  beta, poptheta, options, alpha, 0.0,
                                  dispersion, ancestryAssocTest, true, updateSumLogs);

        // ** Run HMM forward recursions and Sample Locus Ancestry
        TestInd[i]->SampleHiddenStates(options);

        // ** Sample JumpIndicators, update SumLocusAncestry and SumNumArrivals
        TestInd[i]->SampleJumpIndicators(!isGlobalRho);

        // ** Sample individual- or gamete-specific sumintensities
        if (!isGlobalRho)
          TestInd[i]->SampleRho(options, rhoalpha, rhobeta, updateSumLogs);

        // ** update admixture props with conjugate proposal on odd-numbered iterations
        if (iteration % 2)
          TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType,
                                  lambda, NumCovariates, &Covariates,
                                  beta, poptheta, options, alpha, 0.0,
                                  dispersion, ancestryAssocTest, false, updateSumLogs);
      }

      // ** Sample missing values of outcome variable
      //TestInd[i]->SampleMissingOutcomes(&Outcome, R);

    }//end loop over test individuals
  }

  // ** Non-test individuals - conjugate updates only
  #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
    #pragma omp parallel for default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
  #endif
  for (unsigned int i = 0; i < size; ++i) {

    PedBase & el = getElement(i);

    if (Populations > 1) {

      // ** Sample individual- or gamete-specific sumintensities
      if (!isGlobalRho)
        el.SampleRho(options, rhoalpha, rhobeta, updateSumLogs);

      // ** update admixture props with conjugate proposal on odd-numbered iterations
      if (iteration % 2) {
        double DinvLink = 1.0;

	// Accessing global data here, must protect -- let's find a better way?
	#if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP && SAMPLE_THETA_CALL_CRITICAL
	    #pragma omp critical
	#endif
	    { // BEGIN SCOPE BLOCK
	    if(R.size())DinvLink = R[0]->DerivativeInverseLinkFunction(i+i0);
            el.SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType,
                           lambda, NumCovariates, &Covariates,
                           beta, poptheta, options, alpha, DinvLink,
                           dispersion, ancestryAssocTest,
                           el.isPedigree(), updateSumLogs);
	    } // END SCOPE BLOCK
      }
    }

    // ** Sample missing values of outcome variable
    el.SampleMissingOutcomes(&Outcome, R);
  }
}

void AdmixIndividualCollection::setChibNumerator(const AdmixOptions& options, const PedBase::AlphaType & alpha,
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
   getElement(0).setChibNumerator(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat*,*/ &MargLikelihood, A);
}

void AdmixIndividualCollection::updateChib(const AdmixOptions& options, const PedBase::AlphaType & alpha,
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
   getElement(0).updateChib(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat,*/ &MargLikelihood, A);
}

void AdmixIndividualCollection::FindPosteriorModes(const AdmixOptions& options,
						   const vector<bclib::Regression*> &R,
					      const PedBase::AlphaType &alpha, double rhoalpha, double rhobeta,
					      AlleleFreqs* A,
					      const Vector_s& PopulationLabels){
  if(options.getDisplayLevel()>1)
    cout<< endl << "Finding posterior mode of individual parameters ..." << endl;
  //open output file and write header
  ofstream modefile(options.getIndAdmixModeFilename());
  modefile << "Individual \t";
  if(!options.isGlobalRho()){
    if(options.isRandomMatingModel())modefile << "rho0 \t rho1 \t";
    else modefile << "rho \t";
  }
  if(options.isRandomMatingModel()){
    for(int g = 0; g < 2; ++g)
      for(int k = 0; k < options.getPopulations(); ++k) modefile << "mu" <<g<<"."<<PopulationLabels[k]<<"\t";
  }
  else{
    for(int k = 0; k < options.getPopulations(); ++k)modefile << "mu"<<PopulationLabels[k]<<"\t";
  }
  modefile <<endl;

  fill(SumLogTheta, SumLogTheta+options.getPopulations(), 0.0);//reset to 0 as mode-finding function changes it
  //may be unecessary if SumLogTheta zeroed after call to this function(FindPosteriorModes)

  vector<double> lambda;
  vector<const double*> beta;
  for(unsigned i = 0; i < R.size(); ++i){
    lambda.push_back( R[i]->getlambda());
    beta.push_back( R[i]->getbeta());
  }
  if(options.getTestOneIndivIndicator()) {// find posterior mode for test individual only
    TestInd[sizeTestInd-1]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta, A,
					       modefile/*, thetahat, rhohat*/);
  }
  for(unsigned int i = 0; i < size; i++ ){
     getElement(i).FindPosteriorModes(options, alpha, rhoalpha, rhobeta,A,
				  modefile/*, thetahat, rhohat*/);
    modefile << endl;
  }
  modefile.close();
}


// ************** ACCESSORS **************

double AdmixIndividualCollection::GetSumrho()const
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < size; i++ )
      Sumrho += getElement(i).getSumrho();
   return Sumrho;
}


// ************** OUTPUT **************

void AdmixIndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  if(TestInd)
    indadmixoutput->visitIndividual(*(TestInd[sizeTestInd-1]), _locusfortest);
  for(unsigned int i = 0; i < size; i++){
    indadmixoutput->visitIndividual( getElement(i), _locusfortest );
  }
}

void AdmixIndividualCollection::OutputChibResults(bclib::LogWriter& Log) const {
  MargLikelihood.outputResults(Log);
//   Log.setDisplayMode(On);
//   Log << "\nCalculation of Chib algorithm using posterior mode of admixture and sum-intensities, prior mean of allele freqs"
//       << "\nDeviance\t" << -2.0*MargLikelihood.getLogLikelihood()
//       << "\nLogLikelihood\t" << MargLikelihood.getLogLikelihood()
//       << "\nLogPrior\t" << MargLikelihood.getLogPrior()
//       << "\nLogPosterior\t" << MargLikelihood.getLogPosterior()
//       << "\n\nLogMarginalLikelihoodFromChibAlgorithm\t" << MargLikelihood.getLogMarginalLikelihood()
//       << "\n";
}

void AdmixIndividualCollection::getOnePopOneIndLogLikelihood(bclib::LogWriter &Log, const Vector_s& PopulationLabels) {
  Log.setDisplayMode(bclib::On);
  Log << "Log-likelihood for unadmixed " << PopulationLabels[0] << ": "
      << getElement(0).getLogLikelihoodOnePop() << "\n";
}

void AdmixIndividualCollection::accumulateEnergyArrays(const Options& options) {
  #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
    #pragma omp parallel for default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
  #endif
  for(int i = 0; i < sizeTestInd; ++i){ // loop over coolnesses - one copy of test individual at each coolness
    const double Energy = -TestInd[i]->getLogLikelihood(options, true, false); // force HMM update, do not store result
    SumEnergy[i] += Energy;
    SumEnergySq[i] += Energy*Energy;
    TestInd[i]->HMMIsBad(true); // HMM is bad, stored loglikelihood bad
  }
}

void AdmixIndividualCollection::HMMIsBad(bool b){
  if(TestInd)
    for(int i = 0; i < sizeTestInd; ++i)
      TestInd[i]->HMMIsBad(b);
  for(unsigned i = 0; i < size; i++)
    _child[i]->HMMIsBad(b);
}


void AdmixIndividualCollection::ResetChib(){
  MargLikelihood.Reset();
}


void AdmixIndividualCollection::OutputErgodicChib(std::ofstream *avgstream, bool fixedallelefreqs) {
  *avgstream << MargLikelihood.getLogPrior()<< " " << MargLikelihood.getLogPosterior() << " "
	     <<  getElement(0).getLogPosteriorTheta() << " " <<   getElement(0).getLogPosteriorRho()<< " ";
  if(!fixedallelefreqs)
    *avgstream  <<  getElement(0).getLogPosteriorAlleleFreqs() << " "; // not if fixedallelefreqs
  *avgstream << MargLikelihood.getLogMarginalLikelihood();
}


double* AdmixIndividualCollection::getSumEnergy()const{
    return SumEnergy;
}


double* AdmixIndividualCollection::getSumEnergySq()const{
    return SumEnergySq;
}


double AdmixIndividualCollection::getDevianceAtPosteriorMean(
						const Options &			options		,
						vector<bclib::Regression *> &	R		,
						Genome *			Loci		,
						bclib::LogWriter &		Log		,
						const RhoType &			/*SumLogRho*/	,
						unsigned			/*something*/	,
						AlleleFreqs *			/*A*/		) {

  //SumRho = ergodic sum of global sumintensities
  int iterations = options.getTotalSamples()-options.getBurnIn();
  const unsigned numberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();


  // *********************************************************************
  // *********************************************************************
  // *****
  // ***** DDF 2009-06-27:
  // *****
  // ***** This is problematic (i.e. buggy): when this method is called from
  // ***** AdmixMapModel::getDevianceAtPosteriorMean(), which is called from
  // ***** AdmixMapModel::Finalize(), it passes L->getSumLogRho()
  // ***** (i.e. PopAdmix::sumLogRho) in as SumLogRho.  As nearly as I can tell
  // ***** from examining the PopAdmix code as well as observing it while
  // ***** running, PopAdmix::SumLogRho will never have more than 1 element, yet
  // ***** here we access it past the end of its elements provided the genome
  // ***** has more than 1 composite locus.
  // *****
  // ***** This is not overly problematic because the computed RhoBar is passed
  // ***** into Genome::SetLocusCorrelation(), which, it appears to me, promptly
  // ***** ignores it since it ignores any rho vector of length other than 1 (I
  // ***** have added some comments to that method indicating this).  I have
  // ***** therefore #if'd out the buggy code here, although perhaps someday
  // ***** someone should try to figure out what it is trying to do.
  // *****
  // *********************************************************************
  // ***** VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
  #if 0
      // Update chromosomes using globalrho, for globalrho model
      if ( (options.getPopulations() > 1) && options.isGlobalRho() ) {

	RhoType RhoBar(numberOfCompositeLoci);
	for (unsigned i = 0; i < numberOfCompositeLoci; ++i)
	    RhoBar[i] = exp( SumLogRho[i] / double(iterations) );
	//set locus correlation
	Loci->SetLocusCorrelation( RhoBar );
      }
  #endif
  // ***** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // ***** End of buggy+ineffectual code
  // *********************************************************************
  // *********************************************************************



  //set haplotype pair probs to posterior means
  for (unsigned int j = 0; j < numberOfCompositeLoci; ++j)
    (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);

  //set genotype probs using happair probs calculated at posterior means of allele freqs
  setGenotypeProbs(Loci);

  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  const int signed_size = size;
  #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
    #pragma omp parallel for reduction(+:Lhat) default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
  #endif
  for( int i = 0; i < signed_size; i++ ){
    Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }

  Log << bclib::Quiet << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
  for(unsigned c = 0; c < R.size(); ++c){
    double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(iterations, getOutcome(c));
    Lhat += RegressionLogL;
    Log << "DevianceAtPosteriorMean(Regression " << c+1 << ")"
	<< -2.0*RegressionLogL << "\n";
  }

  return(-2.0*Lhat);
}


/// Set the number of decimal places to @a prec
static void ofstreamSettings(ofstream& ofstr, int prec) {
  ofstr << std::setfill(' ');
  ofstr.setf(std::ios::fixed);
  ofstr.precision(prec);
  ofstr.width(prec);
}

#include "AdmixFilenames.h"
///write posterior means of individual admixture params to file
void AdmixIndividualCollection::WritePosteriorMeans(const AdmixOptions& options,
                                                    const vector<string>& PopLabels,
                                                    const Genome& Loci) const {

  const PopIdx K = options.getPopulations();
  const bool isRandomMating = options.isRandomMatingModel();

  ofstream meanfile((options.getResultsDir() + "/" + IND_ADMIXTURE_POSTERIOR_MEANS).c_str());
  ofstreamSettings(meanfile, 3);

  // write header
  meanfile << "ID\t";
  for (PopIdx k = 0; k < K; ++k) {
    meanfile << PopLabels[k];
    if (isRandomMating)
      meanfile << "1";
    meanfile << "\t";
  }
  if (isRandomMating)
    for (PopIdx k = 0; k < K; ++k)
      meanfile << PopLabels[k] << "2\t";
  if (!options.isGlobalRho()) {
    meanfile << "rho";
    if (isRandomMating)
      meanfile << "1\trho2";
  }
  meanfile << endl;

  const unsigned samples = options.getTotalSamples() - options.getBurnIn();
  for (unsigned int i = 0; i < size; ++i) {
    meanfile << getElement(i).getId() << "\t";
    getElement(i).WritePosteriorMeans(meanfile, samples, options.isGlobalRho());
    meanfile << endl;
  }
  meanfile.close();

  // report the admixtures for the X chromosome
  if (Loci.isX_data()) {

    ofstream xmeanfile((options.getResultsDir() + "/" +
                        IND_ADMIXTURE_POSTERIOR_MEANS_XCHR).c_str());
    ofstreamSettings(xmeanfile, 3);

    // write header
    xmeanfile << "ID\t";
    for (PopIdx k = 0; k < K; ++k) {
      xmeanfile << PopLabels[k];
      if (isRandomMating)
        xmeanfile << "1";
      xmeanfile << "\t";
    }
    if (isRandomMating)
      for (PopIdx k = 0; k < K; ++k)
        xmeanfile << PopLabels[k] << "2\t";
    xmeanfile << endl;

    // write the values
    for (unsigned int i = 0; i < size; ++i) {
      xmeanfile << getElement(i).getId() << "\t";
      getElement(i).WritePosteriorMeansXChr(xmeanfile, samples);
      xmeanfile << endl;
    }
    xmeanfile.close();
  }

  if (options.getLocusAncestryProbsIndicator()) {

    // write posterior probs locus ancestry to file
    ofstream locifile((options.getResultsDir() + "/" +
                       LOCUS_ANCESTRY_POSTERIOR_PROBS).c_str());
    ofstreamSettings(locifile, 3);

    // write the values as an R object
    locifile << "structure(c(";
    for (unsigned int i = 0; i < size; ++i) {
      getElement(i).WritePosteriorMeansLoci(locifile);

      // do not append a comma after the values for last individual
      locifile << (i != size - 1 ? ",\n" : "");
    }

    // write the dimensions and their names
    locifile << "), .Dim=c(3, " << K << ", "
             << Loci.GetNumberOfCompositeLoci() << ", " << size << "),"
             << " .Dimnames=list(c(\"0\",\"1\",\"2\"), c(\"";
    for (PopIdx k = 0; k < K; ++k)
      locifile << PopLabels[k] << (k != K - 1 ? "\",\"" : "\"");
    locifile << "), character(0), character(0)))\n";

    locifile.close();
  }
}
