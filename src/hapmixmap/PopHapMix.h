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

#include "HapMixOptions.h"
#include "HapMixGenome.hh"
#include "bclib/HamiltonianMonteCarlo.h"
//#include "bclib/StepSizeTuner.h"//for sampling globalrho and globaltheta
#include "bclib/AdaptiveRejection.h"
#include "bclib/DirichletParamSampler.h"
#include "bclib/DelimitedFileWriter.h"

class InputData;
class IndividualCollection;

///Struct to hold arguments for sampling sumintensities in hapmixmodel
typedef struct {
  unsigned NumBlockStates;
  unsigned NumIntervals; // ? necessary
  //int NumConcordant;
  //int NumDiscordant;
  const int* ConcordanceCounts;
  double Distance;
  const double* theta;
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
  PopHapMix(const HapMixOptions& op, HapMixGenome& loci);
  
  ~PopHapMix();
  
  void Initialise(const string& distanceUnit, bclib::LogWriter& Log);
  void ReadInitialMixturePropsFromFile(const char* initfilename, bclib::LogWriter& Log);
  void ReadInitialArrivalRatesFromFile(const char* initfilename, bclib::LogWriter& Log);  
  
  void SampleArrivalRate(const int* SumAncestry, bool sumlogrho) ;

  //void UpdateGlobalTheta(int iteration, IndividualCollection* individuals);
  void SampleMixtureProportions(const int* SumArrivalCounts);
  
  void SetHMMStateArrivalProbs(bool diploid);

  void OutputParams();
  void OutputParams(bclib::Delimitedostream& out); 
  ///output average mixture props to an ostream
  void OutputAverageMixtureProps(ostream& out)const;
  
  void OutputErgodicAvg( int, std::ofstream& avgstream);

  ///indicates whether the rate parameter is fixed
  bool fixedRateParameter()const;
  const double* getGlobalMixtureProps()const;
  
  void printAcceptanceRates(bclib::LogWriter &Log);
  
  const vector<double>& getlambda()const{return lambda;};
  const vector<double>& getSumLogRho()const{return SumLogLambda;};
  void OutputArrivalRates(const char* filename)const;
  void OutputMixtureProps(const char* filename)const;
  void OutputArrivalRatePosteriorMeans(const char* filename, int samples, const string& distanceUnit)const;
  
private:
  const unsigned K;///< number of subpopulations / block states
  bclib::DelimitedFileWriter outputstream;//output to paramfile

  const HapMixOptions& options;
  HapMixGenome& Loci; 

  std::vector<double> lambda;///< arrival rates
  std::vector<double> SumLogLambda; //ergodic sum of log(rho)

  LambdaArguments LambdaArgs;
  LambdaPriorArguments LambdaPriorArgs;
  bclib::HamiltonianMonteCarlo* ArrivalRateSampler;
  h_args hargs;
  bclib::StepSizeTuner hTuner;  
  bclib::AdaptiveRejection hARS;
  bool fixRateParameter;///< indicates whether to fix the lambda rate parameter

  double* MixtureProps;///<global admixture proportions
  double* MixturePropsPrior;///<parameters of Dirichlet prior on MixureProps
  double* dirparams;//</parameters of Dirichlet distribution from which MixtureProps are sampled
  double* SumTheta;///< sum of mixture props over loci, for calculation of variance
  double* SumThetaSq;///< sum of square of mixture props over loci, for calculation of variance
  double eta;///< observed mixture props precision

  bclib::DirichletParamSampler MixturePropsPrecisionSampler;

  void InitialiseMixtureProportions(bclib::LogWriter& Log);
  void InitialiseArrivalRates(bclib::LogWriter& Log);
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

