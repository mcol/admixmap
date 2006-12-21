#include "HapMixIndividualCollection.h"
#include "Comms.h"

HapMixIndividualCollection::HapMixIndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci){
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

  Individual::SetStaticMembers(Loci, options);
  if(worker_rank < (int)size){
    _child = new Individual*[size];
    for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
      _child[i] = new Individual(i+1, options, Data);//NB: first arg sets Individual's number
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
void HapMixIndividualCollection::SampleLocusAncestry(const AdmixOptions* const options){

  fill(SumAncestry, SumAncestry + 2*NumCompLoci, 0);
#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(15, 0, "Sampleancestry");
#endif

  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** Run HMM forward recursions and sample locus ancestry
    _child[i]->SampleLocusAncestry(options);
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

int HapMixIndividualCollection::getFirstScoreTestIndividualNumber()const{
  if(NumCaseControls > 0)return size - NumCaseControls;
  else return 0;
}
//TODO: alternative for parallel version
void HapMixIndividualCollection::AccumulateConditionalGenotypeProbs(const AdmixOptions* const options, const Genome& Loci){
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
