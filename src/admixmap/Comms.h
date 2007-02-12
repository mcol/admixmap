// *-*-C++-*-*
// header file for Comms class, whcih handles communication between processes in parallel version.
#ifndef COMMS_H
#define COMMS_H 1
#include "config.h"

#ifdef PARALLEL
#include <mpi2c++/mpi++.h>
#include <mpe.h>
#include <vector>
#endif
class Comms{

public:
  static void Initialise();
  static void Finalise();
  static bool isMaster();
  static bool isWorker();
  static bool isFreqSampler();
  static int getNumProcesses(){return NumProcesses;}
  static int getNumWorkers();
  static int getRank(){return global_rank;}
  static int getWorkerRank();
#ifdef PARALLEL
  static void SetDoubleWorkspace(unsigned size, bool isMaster);
  static void SetIntegerWorkspace(unsigned size, bool isMaster);

  static void Reduce(double* x, int size);
  static void Reduce(double* x);
  static void Reduce(int* x);
  static void AllReduce_int(int* x, int size);
  static void BroadcastVector(std::vector<double>& x);
  static void Broadcast(double* x);
  static void Broadcast(double* x, int size);
  static void Broadcast(int* x);

  static void ReduceAncestryCounts(const int* const SumAncestry, int* GlobalSumAncestry, const unsigned size);
  static void reduceAlleleCounts(int* counts, int* globalcounts, const unsigned size);
  static void BroadcastAlleleFreqs(double* const Freqsptr, const int size);
  static void BroadcastRegressionParameters(const double* beta, const int NumCovariates);

  static void ReduceResidualLDScores(const std::vector<std::vector<std::vector<double> > >& Score, 
				     const std::vector<std::vector<std::vector<double> > >& Info, 
				     std::vector<std::vector<std::vector<double> > >& SumScore, 
				     std::vector<std::vector<std::vector<double> > >& SumScore2, 
				     std::vector<std::vector<std::vector<double> > >& SumInfo);
  
  static void ReduceAllelicAssocScores(double** Score, double** Info, unsigned NumLoci, unsigned* sizes, int NumCovars);
  static void ReduceAdmixtureAssocScores(double* Score, double* Info, int size);
#endif

private:
  Comms();//not allowed to instantiate this class
#ifdef PARALLEL
  //communicators
  static MPI::Intracomm workers_and_master;
  static MPI::Intracomm  workers_and_freqs;

  //workspace
  static unsigned max_doubles;
  static double* double_send;
  static double* double_recv;
  static unsigned max_ints;
  static int* int_send;
  static int* int_recv;
#endif
  static int global_rank;
  static int NumProcesses;
};

#endif
