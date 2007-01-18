// *-*-C++-*-*
/** 
 *   HAPMIXMAP
 *   PopHapMix.h 
 *   header file for PopHapMix class
 *   Copyright (c) 2006, 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef POPHAPMIX_H
#define POPHAPMIX_H 1

#include "AdmixOptions.h"
#include "Genome.h"
#include "samplers/HamiltonianMonteCarlo.h"
#include "samplers/StepSizeTuner.h"//for sampling globalrho and globaltheta
#include "samplers/AdaptiveRejection.h"

class InputData;
class IndividualCollection;

///Struct to hold arguments for sampling sumintensities in hapmixmodel
typedef struct {
  unsigned NumPops;
  unsigned NumIntervals; // ? necessary
  int NumConcordant;
  int NumDiscordant;
  double Distance;
  //const double* theta;
  double h;
  double beta;//rate parameter of prior on rho
  double beta_shape;
  double beta_rate;
}LambdaArguments;

///Struct to hold arguments for sampling hyperparameters of sumintensities in hapmixmodel
typedef struct {
    unsigned NumIntervals;
    // const std::vector<double>* rho;
    //double priormeans[3];
    //double priorvars[3];
  double sumlambda;
}LambdaPriorArguments;

///struct to hold args for adaptive rejection sampling of h, the shape multiplier of the prior on expected number of arrivals
typedef struct{
  double sum_lngamma_hd;
  double sum_dloglambda;
  double Dlogbeta;
  double shape;///h has gamma prior with this shape
  double rate;///and this rate
  unsigned NumIntervals;
  double* distances;
}h_args;

///Class to hold and update population admixture and sumintensities parameters and their priors
class PopHapMix
{
public:
  PopHapMix(AdmixOptions* op, Genome* loci);
  
  ~PopHapMix();
  
  void Initialise(const string& distanceUnit, LogWriter& Log);  
  
  void SampleHapMixLambda(const int* SumAncestry, bool sumlogrho) ;

  void UpdateGlobalTheta(int iteration, IndividualCollection* individuals);
  
  void OutputParams(int iteration, LogWriter &Log);
  void OutputParams(ostream* out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const double* getGlobalTheta()const;
  
  void printAcceptanceRates(LogWriter &Log);
  
  //functions required by base class, not implemented
//   void UpdateGlobalSumIntensities(const IndividualCollection* const , bool ){};
//   void UpdatePopAdmixParams(int , const IndividualCollection* const, LogWriter& ){};
//   const double *getpoptheta()const{return 0;};//TODO: check compatibility with regression models

//   const std::vector<std::vector<double> > dummy_alpha;
//   const std::vector<double> &getalpha0()const{return dummy_alpha[0];};
//   const std::vector<std::vector<double> > &getalpha()const{return dummy_alpha;};
//   double getrhoalpha()const{return 0.0;};
//   double getrhobeta()const{return 0.0;};
//   double getglobalrho()const{return 0.0;};
   const vector<double>& getlambda()const{return lambda;};
   const vector<double>& getSumLogRho()const{return SumLogLambda;};
  void OutputLambda(const char* filename)const;
  void OutputLambdaPosteriorMeans(const char* filename, int samples)const;
  
private:
  int K;///< number of subpopulations / block states
  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 

  std::vector<double> lambda;
  std::vector<double> SumLogLambda; //ergodic sum of log(rho)

  LambdaArguments LambdaArgs;
  LambdaPriorArguments LambdaPriorArgs;
  HamiltonianMonteCarlo* HapMixLambdaSampler;
  h_args hargs;
  StepSizeTuner hTuner;  
  AdaptiveRejection hARS;

  double* globaltheta;//global admixture proportions in a hapmixmodel
  //double* globalthetaproposal;//for random walk update
  //StepSizeTuner ThetaTuner;
  //double thetastep;

  void InitializeOutputFile(const string& distanceUnit);
  void ConjugateUpdateGlobalTheta(const vector<int> sumLocusAncestry);
  void UpdateGlobalThetaWithRandomWalk(IndividualCollection* IC);
  void Accept_Reject_Theta( double logpratio, int Populations);
  void SampleRateParameter();
  void Sampleh_RandomWalk();
  void Sampleh_ARS();

  static double LambdaEnergy(const double* const x, const void* const vargs);
  static void LambdaGradient( const double* const x, const void* const vargs, double* g );

  static double hlogf(double h, const void* const vargs);
  static double hdlogf(double h, const void* const vargs);
  static double hd2logf(double h, const void* const vargs);

  // UNIMPLEMENTED
  // to avoid use
  PopHapMix();
  PopHapMix(const PopHapMix&);
  PopHapMix& operator=(const PopHapMix&);

};

#endif /* !defined POPHAPMIX_H */

