// This file contains the communication code for the parallel version of ADMIXMAP
// 
#include "Comms.h"
#include <algorithm>

int Comms::global_rank = 0;
int Comms::NumProcesses = 1;

bool Comms::isMaster(){
  return((bool)(global_rank==0));
}
int Comms::getNumWorkers(){
  if(NumProcesses > 2)
    return NumProcesses-2;
  else if(NumProcesses == 1)//serial version
    return 1;
  else return 0;//trying to run on two processors
}
int Comms::getWorkerRank(){
  if(NumProcesses == 1)
    return global_rank;
  else return global_rank-2;
}
#ifdef PARALLEL
bool Comms::isFreqSampler(){
  return((bool)(global_rank==1));
}
bool Comms::isWorker(){
  return((bool)(global_rank>1));
}

MPI::Intracomm Comms::workers_and_master;
MPI::Intracomm  Comms::workers_and_freqs;
unsigned Comms::max_doubles = 0;
double* Comms::double_send = 0;
double* Comms::double_recv = 0;
unsigned Comms::max_ints = 0;
int* Comms::int_send = 0;
int* Comms::int_recv = 0;

void Comms::Initialise(){
  global_rank = MPI::COMM_WORLD.Get_rank();
  NumProcesses = MPI::COMM_WORLD.Get_size();
  if(NumProcesses <= 2)
    throw("ERROR: too few processes\n");

  //set up communicators
  int ranks[1];
  MPI::Group world_group = MPI::COMM_WORLD.Get_group();
  ranks[0] = 1;
  MPI::Group master_group = world_group.Excl(1, ranks);//exclude process (1) that will update freqs
  ranks[0] = 0;
  MPI::Group freqs_group = world_group.Excl(1, ranks);//exclude process 0 (master)
  workers_and_master = MPI::COMM_WORLD.Create(master_group);
  workers_and_freqs = MPI::COMM_WORLD.Create(freqs_group);
  world_group.Free();
  master_group.Free();
  freqs_group.Free();
}

void Comms::SetDoubleWorkspace(unsigned size, bool isMaster){
  //expand double workspace to size, if necessary. isMaster indicates if this process is a receive (no need to allocate receive space if not) 
  if(size > max_doubles){
    max_doubles = size;
    delete[] double_send;
    double_send = new double[max_doubles];
    if(isMaster){
      delete[] double_recv;
      double_recv = new double[max_doubles];
    }
  }
}
void Comms::SetIntegerWorkspace(unsigned size, bool isMaster){
  //expand int workspace to size, if necessary. isMaster indicates if this process is a receive (no need to allocate receive space if not) 
  if(size > max_ints){
    max_ints = size;
    delete[] int_send;
    int_send = new int[max_ints];
    if(isMaster){
      delete[] int_recv;
      int_recv = new int[max_ints];
    }
  }
}

void Comms::Finalise(){
  //free workspace
  delete[] double_send;
  delete[] double_recv;
  delete[] int_send;
  delete[] int_recv;

  //free communicators
  if(global_rank!=1)workers_and_master.Free();
  if(global_rank!=0)workers_and_freqs.Free();
}

// *** IndividualCollection: reduction of LocusAncestry
void Comms::ReduceAncestryCounts(const int* const SumAncestry, int* GlobalSumAncestry, const unsigned size){
  //if(workers_and_master.Get_rank()==0)fill(GlobalSumAncestry, GlobalSumAncestry + size, 0);
  MPE_Log_event(1, 0, "BarrierStart");
  workers_and_master.Barrier();
  MPE_Log_event(2, 0, "Barrierend");
  MPE_Log_event(3, 0, "RedAncStart");    
  workers_and_master.Reduce(SumAncestry, GlobalSumAncestry, size, MPI::INT, MPI::SUM, 0); 
  MPE_Log_event(4, 0, "RedAncEnd");
}

// *** AlleleFreqs: reduction of allele counts
void Comms::reduceAlleleCounts(int* counts, int* globalcounts, const unsigned size){
//synchronise processes
    MPE_Log_event(1, 0, "CountsBarrier");
    workers_and_freqs.Barrier();
    MPE_Log_event(2, 0, "CountsBarrierEnd");
    MPE_Log_event(5, 0, "RedCountstart");
    workers_and_freqs.Reduce(counts, globalcounts, size, MPI::INT, MPI::SUM, 0);
    MPE_Log_event(6, 0, "RedCountend");

    //put totals back into AlleleCounts on top process, by swapping addresses
    if(isFreqSampler()){
      int* temp = counts;
      counts = globalcounts;
      globalcounts = temp;
    }
}

// ** AlleleFreqs : broadcast frequencies
void Comms::BroadcastAlleleFreqs(double* const Freqsptr, const int size){
  MPE_Log_event(1, 0, "FreqsBarrier");
  workers_and_freqs.Barrier();
  MPE_Log_event(2, 0, "FreqsBarrierEnd");
  MPE_Log_event(19, 0, "BcastFreqs"); 
  workers_and_freqs.Bcast(Freqsptr, size, MPI::DOUBLE, 0 );
  MPE_Log_event(20, 0, "FreqsBcasted"); 
}

// ** Latent: broadcast sumintensities
void Comms::BroadcastRho(std::vector<double>& rho){
  MPE_Log_event(17, 0, "Bcastrho");
  workers_and_master.Barrier();
  workers_and_master.Bcast(&(*(rho.begin())), rho.size(), MPI::DOUBLE, 0);
  MPE_Log_event(18, 0, "Bcasted");
}

//Regression:: broadcsat regression parameters
void Comms::BroadcastRegressionParameters(const double* beta, const int NumCovariates){
  workers_and_master.Barrier();
  workers_and_master.Bcast(const_cast<double*>(beta), NumCovariates, MPI::DOUBLE, 0);
}
void Comms::ReduceLogLikelihood(double* LogLik){
  double globalLogLik = 0.0;
  workers_and_master.Barrier();
  workers_and_master.Reduce(LogLik, &globalLogLik, 1, MPI::DOUBLE, MPI::SUM, 0);
  *LogLik = globalLogLik;
}

void Comms::ReduceResidualLDScores(const std::vector<std::vector<std::vector<double> > >& Score, 
				   const std::vector<std::vector<std::vector<double> > >& Info, 
				   std::vector<std::vector<std::vector<double> > >& SumScore, 
				   std::vector<std::vector<std::vector<double> > >& SumScore2, 
				   std::vector<std::vector<std::vector<double> > >& SumInfo){
  //TODO: ?define define MPI type for score and info to save having to pack and unpack


  //pack score into array ready to send
  int count = 0;
  for(unsigned c = 0; c < Score.size(); ++c){//number of chromosomes
    for(unsigned k = 0; k < Score[c].size(); ++k){//number of comp loci on this chromosome
      copy(Score[c][k].begin(), Score[c][k].end(), double_send+count);
      count += Score[c][k].size();
    }
  }
  //reduce into receive arrays on master process
  workers_and_master.Barrier();
  workers_and_master.Reduce(double_send, double_recv, count, MPI::DOUBLE, MPI::SUM, 0);

  if(workers_and_master.Get_rank()==0){
    //accumulate score, square of score on master process
    int scoreindex = 0; 
    for(unsigned c = 0; c < Score.size(); ++c)
      for(unsigned k = 0; k < Score[c].size(); ++k){
	unsigned dim = Score[c][k].size();

	for(unsigned j = 0; j < dim; ++j){
	  SumScore[c][k][j] += double_recv[scoreindex + j];
	  for(unsigned jj = 0; jj < dim; ++jj){
	    SumScore2[c][k][j*dim +jj] += double_recv[scoreindex + j]*double_recv[scoreindex + jj];
	  }
	}
	scoreindex += dim;
      }

  }
  //pack info into array ready to send
  count = 0;
  for(unsigned c = 0; c < Info.size(); ++c){//number of chromosomes
    for(unsigned k = 0; k < Info[c].size(); ++k){//number of comp loci on this chromosome
      unsigned dim = Info[c][k].size();// (#states in this comp locus -1 )*(#states in next comp locus -1)
      copy(Info[c][k].begin(), Info[c][k].end(), double_send+count);
      count += dim*dim;
    }
  }

  workers_and_master.Barrier();
  workers_and_master.Reduce(double_send, double_recv, count, MPI::DOUBLE, MPI::SUM, 0);

  if(workers_and_master.Get_rank()==0){
    //accumulate info on master process
    int infoindex = 0;
    for(unsigned c = 0; c < Info.size(); ++c)
      for(unsigned k = 0; k < Info[c].size(); ++k){
	unsigned dim = Info[c][k].size();

	for(unsigned j = 0; j < dim; ++j){
	  for(unsigned jj = 0; jj < dim; ++jj){
	    SumInfo[c][k][j*dim +jj] += double_recv[infoindex + j*dim + jj];
	  }
	}
	infoindex+= dim*dim;
      }
  }
}

void Comms::ReduceAdmixtureAssocScores(double* Score, double* Info, int size){
  workers_and_master.Reduce(Score, double_send, size, MPI::DOUBLE, MPI::SUM, 0);
  if(workers_and_master.Get_size()==0)std::copy(double_send, double_send+size, Score);
  workers_and_master.Reduce(Info, double_send, size, MPI::DOUBLE, MPI::SUM, 0);
  if(workers_and_master.Get_size()==0)std::copy(double_send, double_send+size, Info);
}

void Comms::ReduceAllelicAssocScores(double** Score, double** Info, unsigned NumLoci, unsigned* sizes, int NumCovars){
  //pack score into send array
  int count = 0;
  for(unsigned int j = 0; j < NumLoci; j++ ){
    std::copy(Score[j], Score[j] + sizes[j]+NumCovars, double_send+count);
    count += sizes[j]+NumCovars;
  }

  //sum over individuals on master process
  workers_and_master.Barrier();
  workers_and_master.Reduce(double_send, double_recv, count, MPI::DOUBLE, MPI::SUM, 0);
  
  //unpack
  if(workers_and_master.Get_rank()==0){
    count = 0;
    for(unsigned int j = 0; j < NumLoci; j++ ){
      std::copy(double_recv+count, double_recv+count+ sizes[j]+NumCovars, Score[j]);
      count += sizes[j]+NumCovars;
    }
  }


  //pack info into send array
  count = 0;
  for(unsigned int j = 0; j < NumLoci; j++ ){
    std::copy(Info[j], Info[j] + (sizes[j]+NumCovars)*(sizes[j]+NumCovars), double_send+count);
    count += (sizes[j]+NumCovars)*(sizes[j]+NumCovars);
  }
  //sum over individuals on master process
  workers_and_master.Barrier();
  workers_and_master.Reduce(double_send, double_recv, count, MPI::DOUBLE, MPI::SUM, 0);
  //unpack
  if(workers_and_master.Get_rank()==0){
    count = 0;
    for(unsigned int j = 0; j < NumLoci; j++ ){
      std::copy(double_recv+count, double_recv+count+ (sizes[j]+NumCovars)*(sizes[j]+NumCovars), Info[j] );
      count += (sizes[j]+NumCovars)*(sizes[j]+NumCovars);
    }
  }
}

void Comms::AllReduce_int(int* x, int size){
  workers_and_master.Barrier();
  workers_and_master.Allreduce(x, int_send, size, MPI::INT, MPI::SUM);
  std::copy(int_send, int_send+size, x);
}

#else //serial version
bool Comms::isFreqSampler(){
  return true;
}
bool Comms::isWorker(){
  return true;
}
#endif
