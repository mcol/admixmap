#include "HapMixIndividualCollection.h"


HapMixIndividualCollection::HapMixIndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci) :
  IndividualCollection(options, Data, Loci){
  GlobalSumAncestry = 0;
  SumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#ifdef PARALLEL
  if(Comms::isMaster())GlobalSumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#endif
}

HapMixIndividualCollection::~HapMixIndividualCollection(){
  delete[] SumAncestry;
#ifdef PARALLEL
  delete[] GlobalSumAncestry;
#endif
}
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

