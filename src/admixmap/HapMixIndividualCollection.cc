#include "HapMixIndividualCollection.h"
#include "HapMixOptions.h"
#include "HapMixIndividual.h"
#include "HapMixFreqs.h"
#include "regression/Regression.h"
#include "Comms.h"

HapMixIndividualCollection::HapMixIndividualCollection(const HapMixOptions* const options, const InputData* const Data, IGenome* Loci, const HapMixFreqs* A){
  SetNullValues();
  GlobalSumAncestry = 0;
  SumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#ifdef PARALLEL
  if(Comms::isMaster())GlobalSumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#endif
  Populations = options->getPopulations();
  NumInd = Data->getNumberOfIndividuals();//number of individuals, including case-controls
  size = NumInd;
  NumCaseControls = Data->getNumberOfCaseControlIndividuals();
  NumCompLoci = Loci->GetNumberOfCompositeLoci();
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
      _child[i] = new HapMixIndividual(i+1, options, Data);//NB: first arg sets Individual's number
    }
  }
  if(options->OutputCGProbs())GPO.Initialise(options->GetNumMaskedIndividuals(), options->GetNumMaskedLoci());
}

HapMixIndividualCollection::~HapMixIndividualCollection(){
  delete[] SumAncestry;
#ifdef PARALLEL
  delete[] GlobalSumAncestry;
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
void HapMixIndividualCollection::SampleLocusAncestry(const Options* const options){

  fill(SumAncestry, SumAncestry + 2*NumCompLoci, 0);
#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(15, 0, "Sampleancestry");
#endif

  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** Run HMM forward recursions and sample locus ancestry
    _child[i]->SampleLocusAncestry(options);

    if (
      // If it's a control individual
      i >= getFirstScoreTestIndividualNumber()
      // And if the score tests are switched on
      && strlen(options->getAllelicAssociationScoreFilename()) > 0
      // FIXME: The next condition shouldn't be necessary.
      && (not _child[i]->isHaploidIndividual()))
    {
      _child[i]->calculateUnorderedGenotypeProbs();
    }

    _child[i]->AccumulateAncestry(SumAncestry);
  }
#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(16, 0, "Sampledancestry");
  Comms::ReduceAncestryCounts(SumAncestry, GlobalSumAncestry, 2*NumCompLoci);
#endif
}

const int* HapMixIndividualCollection::getSumAncestry()const{
#ifdef PARALLEL
  return GlobalSumAncestry;
#else
  return SumAncestry;
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
//TODO: alternative for parallel version
void HapMixIndividualCollection::AccumulateConditionalGenotypeProbs(const HapMixOptions* const options, const IGenome& Loci){
  const std::vector<unsigned>& MaskedLoci = options->getMaskedLoci();
  const std::vector<unsigned>& MaskedIndividuals = options->getMaskedIndividuals();
  int anc[2];
  unsigned j = 0;
  //NB: indices count from 1 so must be offset by -1
  for(std::vector<unsigned>::const_iterator locus = MaskedLoci.begin(); locus!= MaskedLoci.end(); ++j, ++locus)
    if (*locus <= Loci.GetNumberOfCompositeLoci()){
      unsigned i = 0;
      for(std::vector<unsigned>::const_iterator indiv = MaskedIndividuals.begin(); indiv!= MaskedIndividuals.end(); ++i, ++indiv)
	if(*indiv <= size){
	  _child[*indiv-1]->GetLocusAncestry(*locus-1, anc);
	  GPO.Update(i, j, Loci(*locus-1), _child[*indiv-1]->getPossibleHapPairs(*locus-1), anc);
    }
  }
}
void HapMixIndividualCollection::OutputCGProbs(const char* filename){
  GPO.Output(filename);
}

double HapMixIndividualCollection::getDevianceAtPosteriorMean(
    const Options* const options, vector<Regression *> &R, IGenome* Loci,
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
