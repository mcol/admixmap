// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Latent.h 
 *   header file for Latent class
 *   Copyright (c) 2002-2006 LSHTM
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef LATENT_H
#define LATENT_H 1

#include "common.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>

#include "samplers/rand.h"

#include "AdmixOptions.h"
#include "Individual.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "utils/LogWriter.h"

#include "samplers/StepSizeTuner.h"//for sampling globalrho and globaltheta
#include "samplers/DirichletParamSampler.h"//for sampling pop admix

class InputData;
class IndividualCollection;

///Struct to hold arguments for sampling sumintensities in hapmixmodel
typedef struct {
  unsigned NumPops;
  unsigned NumIntervals; // ? necessary
  const int* SumAncestry;
  double Distance;
  //const double* theta;
  double h;
  double beta;
  double beta_shape;
  double beta_rate;

//   double sumrho; // ? necessary
//   double sumlogrho; // ? necessary
}LambdaArguments;

///Struct to hold arguments for sampling hyperparameters of sumintensities in hapmixmodel
typedef struct {
    unsigned NumIntervals;
    // const std::vector<double>* rho;
    //double priormeans[3];
    //double priorvars[3];
  double sumlambda;
//  double sumlogrho;
}LambdaPriorArguments;

///Class to hold and update population admixture and sumintensities parameters and their priors
class Latent
{
public:
  Latent( AdmixOptions*, Genome* );
  
  ~Latent();
  
  void Initialise(int Numindividuals, const Vector_s&  PopulationLabels, LogWriter& Log);  
  
  void InitializeOutputFile(const Vector_s&  PopulationLabels);
  
  void UpdateGlobalSumIntensities(const IndividualCollection* const IC, bool sumlogtheta);

  void SampleHapMixLambda(const int* SumAncestry, bool sumlogrho) ;

  void UpdatePopAdmixParams(int iteration, const IndividualCollection* const, LogWriter &Log);
  void UpdateGlobalTheta(int iteration, IndividualCollection* individuals);
  
  void OutputParams(int iteration, LogWriter &Log);
  void OutputParams(ostream* out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const std::vector<double> &getalpha0()const;
  const std::vector<std::vector<double> > &getalpha()const;
  double getrhoalpha()const;
  double getrhobeta()const;
  double getglobalrho()const;
  const vector<double>& getrho()const;
  const vector<double>& getSumLogRho()const;
  const double *getpoptheta()const;
  const double* getGlobalTheta()const;
  
  void printAcceptanceRates(LogWriter &Log);
  
  float getAlphaSamplerAcceptanceRate()const;
  float getAlphaSamplerStepsize()const;
  
  double getRhoSamplerAccRate()const;
  double getRhoSamplerStepsize()const;
  void resetStepSizeApproximator(int k); 
  
private:
  
  double sampleForRho();
  void OpenOutputFiles();
  
  int K;///< number of subpopulations / block states
  std::vector<double> rho;
  std::vector<double> rhoproposal;
  double rhoalpha;
  double rhobeta;
  double rhobeta0;
  double rhobeta1;
  std::vector<double> SumLogRho; //ergodic sum of log(rho)

  LambdaArguments HapMixLambdaArgs;
  LambdaPriorArguments LambdaPriorArgs;
  HamiltonianMonteCarlo* HapMixLambdaSampler;
  
  //RWM sampler for global rho
  StepSizeTuner TuneRhoSampler;
  int w, NumberOfUpdates;
  double step;
  double step0;

  std::vector<std::vector<double> > alpha; //population admixture Dirichlet parameters
  std::vector<double> SumAlpha; //ergodic sums of alphas
  //sampler for alpha
  DirichletParamSampler PopAdmixSampler;
  double *poptheta;    //ergodic average of population admixture, used to centre the values of individual admixture 
                       //in the regression model

  double* globaltheta;//global admixture proportions in a hapmixmodel
  //double* globalthetaproposal;//for random walk update
  //StepSizeTuner ThetaTuner;
  //double thetastep;

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 

  void ConjugateUpdateGlobalTheta(const vector<int> sumLocusAncestry);
  void UpdateGlobalThetaWithRandomWalk(IndividualCollection* IC);
  void Accept_Reject_Theta( double logpratio, int Populations);
  void SampleHapmixLambdaPriorParameters();

  static double LambdaEnergy(const double* const x, const void* const vargs);
  static void LambdaGradient( const double* const x, const void* const vargs, double* g );

  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */

