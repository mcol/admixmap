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

#include "rand.h"

#include "AdmixOptions.h"
#include "Individual.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "LogWriter.h"

#include "StepSizeTuner.h"//for sampling globalrho and globaltheta
#include "DirichletParamSampler.h"//for sampling pop admix

class InputData;
class IndividualCollection;

class Latent
{
public:
  Latent( AdmixOptions*, const Genome* const);
  
  ~Latent();
  
  void Initialise(int Numindividuals, const std::string* const PopulationLabels, LogWriter& Log);  
  
  void InitializeOutputFile(const std::string* const);
  
  void UpdateGlobalSumIntensities(const IndividualCollection* const IC, Chromosome **C);
  void SampleSumIntensities(const std::vector<unsigned> &SumNumArrivals, unsigned n, Chromosome **C); 
  void UpdatePopAdmixParams(int iteration, const IndividualCollection* const, LogWriter &Log,  bool anneal);
  void UpdateGlobalTheta(int iteration, IndividualCollection* individuals, Chromosome** C);
  
  void OutputParams(int iteration, LogWriter &Log);
  void OutputParams(ostream* out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const std::vector<double> &getalpha0()const;
  const std::vector<std::vector<double> > &getalpha()const;
  double getrhoalpha()const;
  double getrhobeta()const;
  double getglobalrho()const;
  const vector<double>& getrho()const;
  double getSumLogRho()const;
  const double *getpoptheta()const;

  void printAcceptanceRates(LogWriter &Log);
  
  float getAlphaSamplerAcceptanceRate()const;
  float getAlphaSamplerStepsize()const;
  
  double getRhoSamplerAccRate()const;
  double getRhoSamplerStepsize()const;
  
private:
  
  double sampleForRho();
  
  void OpenOutputFiles();
  
  /*
   * rho is the sumofintensities parameter. It has a beta prior with shape and scale parameters rhoalpha and rhobeta.
   * rhobeta has a beta hyperprior with parameters rhobeta0 and rhobeta1
   */
  int K;
  std::vector<double> rho; 
  double rhoalpha;
  double rhobeta; 
  double rhobeta0;
  double rhobeta1;
  double SumLogRho; //ergodic sum of log(rho)
  
  //RWM sampler for global rho
  StepSizeTuner TuneRhoSampler;
  int w, NumberOfUpdates;
  double step, step0;

  std::vector<std::vector<double> > alpha; //population admixture Dirichlet parameters
  std::vector<double> SumAlpha; //ergodic sums of alphas
  //sampler for alpha
  
  DirichletParamSampler PopAdmixSampler;

  double *poptheta;    //ergodic average of population admixture, used to centre the values of individual admixture 
                       //in the regression model

  double* globaltheta;//global admixture proportions in a hapmixmodel
  double* globalthetaproposal;//for random walk update
  StepSizeTuner ThetaTuner;
  double thetastep;

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  const Genome* Loci; 

  void ConjugateUpdateGlobalTheta(const vector<int> sumLocusAncestry);
  void UpdateGlobalThetaWithRandomWalk(IndividualCollection* IC, Chromosome** C);
  void Accept_Reject_Theta( double logpratio, int Populations);

  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */

