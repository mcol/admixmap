#include "HapMixIndividualCollection.h"
#include "HapMixOptions.h"
#include "HapMixIndividual.h"
#include "HapMixFreqs.h"
#include "regression/Regression.h"
#include "Comms.h"

HapMixIndividualCollection
::HapMixIndividualCollection(const HapMixOptions* const options, 
			     const InputData* const Data, Genome* Loci, const HapMixFreqs* A, const double* theta){
  SetNullValues();
  NumCompLoci = Loci->GetNumberOfCompositeLoci();
  GlobalConcordanceCounts = 0;
  GlobalSumArrivalCounts = 0;
  ConcordanceCounts = new int[NumCompLoci*2*(options->getNumberOfBlockStates())];
  SumArrivalCounts = new int[NumCompLoci*options->getNumberOfBlockStates()];
#ifdef PARALLEL
  if(Comms::isMaster()){
    GlobalConcordanceCounts = new int[NumCompLoci*2*(options->getNumberOfBlockStates())];
    GlobalSumArrivalCounts = new int[NumCompLoci*options->getNumberOfBlockStates()];
  }
#endif
  Populations = options->getPopulations();
  NumInd = Data->getNumberOfIndividuals();//number of individuals, including case-controls
  size = NumInd;
  NumCaseControls = Data->getNumberOfCaseControlIndividuals();
  worker_rank = 0;
  NumWorkers = 1;
#ifdef PARALLEL
  int global_rank = MPI::COMM_WORLD.Get_rank();
  //create communicator for workers and find size of and rank within this group
  workers = MPI::COMM_WORLD.Split( Comms::isWorker(), global_rank);
  NumWorkers = workers.Get_size();
  if(global_rank >1)
    worker_rank = workers.Get_rank();
  else worker_rank = size;//so that non-workers will not loop over Individuals
#endif
  
  //  Individual::SetStaticMembers(Loci, options);
  HapMixIndividual::SetStaticMembers(Loci, options, A->getHaploidGenotypeProbs(), A->getDiploidGenotypeProbs());

  if(worker_rank < (int)size){
    _child = new Individual*[size];
    for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
      // _child[i] = new Individual(i+1, options, Data);//NB: first arg sets Individual's number
      _child[i] = new HapMixIndividual(i+1, options, Data, theta);//NB: first arg sets Individual's number
    }
  }
  if(options->OutputCGProbs())GPO.Initialise(options->GetNumMaskedIndividuals(), options->GetNumMaskedLoci());
}

HapMixIndividualCollection::~HapMixIndividualCollection(){
  delete[] ConcordanceCounts;
#ifdef PARALLEL
  delete[] GlobalConcordanceCounts;
#endif
}
// Individual* HapMixIndividualCollection::getIndividual(int num)const
// {
//   if (num < (int)size){
//     return _child[num];
//   } else {
//     return 0;
//   }
// }
void HapMixIndividualCollection::SampleHiddenStates(const HapMixOptions* const options, unsigned iteration){

#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(15, 0, "Sampleancestry");
#endif

  fill(ConcordanceCounts, ConcordanceCounts + NumCompLoci*(options->getPopulations()+1), 0);
  fill(SumArrivalCounts, SumArrivalCounts + NumCompLoci*options->getPopulations(), 0);

  const vector<unsigned>& maskedIndividuals = options->getMaskedIndividuals();
  const vector<unsigned>::const_iterator mi_begin = maskedIndividuals.begin();
  const vector<unsigned>::const_iterator mi_end = maskedIndividuals.end();

  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** Run HMM forward recursions and sample hidden states
    _child[i]->SampleLocusAncestry(options);

    if (
	(//If it's after the burnin 
	 (int)iteration > options->getBurnIn()
	 // and it's a case or control individual
	 && isCaseControl(i)
	 // and if the score tests are switched on
	 && options->getTestForAllelicAssociation()
	 // FIXME: The next condition shouldn't be necessary, but it is.
	 && (not _child[i]->isHaploidIndividual())
	 )
	//or it's a masked individual
	// maskedIndividuals indices are 1-based, offset of 1 is needed
	|| (find(mi_begin, mi_end, i+1) != mi_end)
	)
    {
      _child[i]->calculateUnorderedGenotypeProbs();
    } 

    //accumulate sufficient statistics for update of arrival rates and mixture proportions
    _child[i]->AccumulateConcordanceCounts(ConcordanceCounts);
    _child[i]->SampleJumpIndicators(SumArrivalCounts);
  }
#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(16, 0, "Sampledancestry");
  Comms::ReduceAncestryCounts(ConcordanceCounts, GlobalConcordanceCounts, NumCompLoci*2*(options->getPopulations()));
  //TODO: reduce SumArrivalCounts
#endif
}

const int* HapMixIndividualCollection::getConcordanceCounts()const{
#ifdef PARALLEL
  return GlobalConcordanceCounts;
#else
  return ConcordanceCounts;
#endif
}

const int* HapMixIndividualCollection::getSumArrivalCounts()const{
#ifdef PARALLEL
  return GlobalSumArrivalCounts;
#else
  return SumArrivalCounts;
#endif
}

int HapMixIndividualCollection::getNumberOfIndividualsForScoreTests()const{
  if(NumCaseControls > 0)return NumCaseControls;
  else return size;
}

unsigned int HapMixIndividualCollection::getFirstScoreTestIndividualNumber()const{
  if(NumCaseControls > 0)return size - NumCaseControls;
  else return 0;
}

///determines if individual i is a case/control ie its genotype came from ccgenotypesfile
bool HapMixIndividualCollection::isCaseControl(unsigned i)const{
  return ( i > (size - NumCaseControls) );
}

//TODO: alternative for parallel version
void HapMixIndividualCollection::AccumulateConditionalGenotypeProbs(const HapMixOptions* const options, const Genome& Loci){
  
  const std::vector<unsigned>& MaskedLoci = options->getMaskedLoci();
  const std::vector<unsigned>& MaskedIndividuals = options->getMaskedIndividuals();
  // vu_ci stands for vector<unsigned>::const_iterator
  typedef std::vector<unsigned>::const_iterator vu_ci;
  unsigned j = 0;
  //NB: indices count from 1 so must be offset by -1
  for(vu_ci locus_i = MaskedLoci.begin(); locus_i!= MaskedLoci.end(); ++j, ++locus_i) {
    if (*locus_i <= Loci.GetNumberOfCompositeLoci()){
      unsigned i = 0;
      for(vu_ci indiv_i = MaskedIndividuals.begin(); indiv_i!= MaskedIndividuals.end(); ++i, ++indiv_i) {
        if(*indiv_i <= size) {
          GPO.Update(i, j, _child[(*indiv_i) - 1]->getUnorderedProbs((*locus_i) - 1));
        }
      }
    }
  }
}
void HapMixIndividualCollection::OutputCGProbs(const char* filename){
  GPO.Output(filename);
}

double HapMixIndividualCollection
::getDevianceAtPosteriorMean(const Options* const options, vector<Regression *> &R, Genome* Loci,
                             LogWriter &Log, const vector<double>& SumLogRho, unsigned numChromosomes
                             , AlleleFreqs* 
#ifdef PARALLEL
//A is not required in serial version
    A
#endif
    ){

  //TODO: broadcast SumLogRho to workers
  //SumRho = ergodic sum of global sumintensities
  int iterations = options->getTotalSamples()-options->getBurnIn();
  
  //update chromosomes using globalrho, for globalrho model
  if(options->getPopulations() >1 && (options->isGlobalRho() || options->getHapMixModelIndicator()) ){
    vector<double> RhoBar(Loci->GetNumberOfCompositeLoci());
    if(Comms::isMaster())//master only
      for(unsigned i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)RhoBar[i] = (exp(SumLogRho[i] / (double)iterations));
#ifdef PARALLEL
    if(!Comms::isFreqSampler()) 
      Comms::BroadcastVector(RhoBar);
#endif
    //set locus correlation
    if(Comms::isWorker()){//workers only
      Loci->SetLocusCorrelation(RhoBar);
      if(options->getHapMixModelIndicator())
	for( unsigned int j = 0; j < numChromosomes; j++ )
	  //set global state arrival probs in hapmixmodel
	  //TODO: can skip this if xonly analysis with no females
	  //NB: assumes always diploid in hapmixmodel
	  //KLUDGE: should use global theta as first arg here; Theta in Individual should be the same
	  Loci->getChromosome(j)->SetStateArrivalProbs(options->isRandomMatingModel(), true);
    }
  }
  
  //set haplotype pair probs to posterior means (in parallel version, sets AlleleProbs(Freqs) to posterior means
  if(Comms::isFreqSampler())
    for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
      (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);
  
#ifdef PARALLEL
  //broadcast allele freqs
  if(!Comms::isMaster())A->BroadcastAlleleFreqs();
#endif

  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  //if(Comms::isWorker())setGenotypeProbs(Loci, A);

  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    if(options->getHapMixModelIndicator())Lhat += _child[i]->getLogLikelihood(options, false, false);
    else Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }
#ifdef PARALLEL
  if(!Comms::isFreqSampler()){
    Comms::Reduce(&Lhat);
  }
#endif

  if(Comms::isMaster()){
    Log << Quiet << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
    for(unsigned c = 0; c < R.size(); ++c){
      double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(iterations, getOutcome(c));
      Lhat += RegressionLogL;
      Log << "DevianceAtPosteriorMean(Regression " << c+1 << ")"
	  << -2.0*RegressionLogL << "\n";
    }
  }
  return(-2.0*Lhat);
}
