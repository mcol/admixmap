// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Latent.h 
 *   header file for Latent class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef LATENT_H
#define LATENT_H 1

// ** define which sampler to use for pop admixture Dirichlet parameters
#define POPADMIXSAMPLER 2 //1 = original Adaptive Rejection sampler, 
                          //2 = DirichletParamSampler, 
                          //3 = HamiltonianMonteCarlo

#include "common.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_cdf.h>

#include "rand.h"

#include "AdmixOptions.h"
#include "Individual.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "LogWriter.h"

#include "StepSizeTuner.h"

#if POPADMIXSAMPLER==1
#include "AdaptiveRejection.h"
#endif

#if POPADMIXSAMPLER == 2
#include "DirichletParamSampler.h"
#elif POPADMIXSAMPLER == 3
#include "HamiltonianMonteCarlo.h"

typedef struct{
  int dim;
  int n;
  double eps0;
  double eps1;
  const double* sumlogtheta;
}AlphaSamplerArgs;

#endif

class InputData;
class IndividualCollection;

class Latent
{
public:
  Latent( AdmixOptions*, const Genome* const, LogWriter *);
  
  ~Latent();
  
  void Initialise(int Numindividuals, const std::string* const PopulationLabels);  
  
  void InitializeOutputFile(const std::string* const);
  
  //Updating every iteration
  void UpdateRhoWithRW(const IndividualCollection* const IC, Chromosome **C, double LogL);
  void Update(int iteration, const IndividualCollection* const, bool anneal);
  
  void OutputParams(int iteration);
  void OutputParams(int iteration, ostream* out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const std::vector<double> &getalpha0()const;
  const std::vector<std::vector<double> > &getalpha()const;
  double getrhoalpha()const;
  double getrhobeta()const;
  double getrho()const;
  double getSumLogRho()const;
  const double *getpoptheta()const;
  
#if POPADMIXSAMPLER == 2 
  float getEtaSamplerAcceptanceRate()const;
  float getEtaSamplerStepsize()const;
  float getMuSamplerAcceptanceRate()const;
  float getMuSamplerStepsize()const;
#endif
  
  //#if POPADMIXSAMPLER == 3 
  float getAlphaSamplerAcceptanceRate()const;
  float getAlphaSamplerStepsize()const;
  //#endif
  
  double getRhoSamplerAccRate()const;
  double getRhoSamplerStepsize()const;
  
private:
  
  double sampleForRho();
  
  void OpenOutputFiles();
  
  /*
   * rho is the sumofintensities parameter. It has a beta prior with shape and scale parameters rhoalpha and rhobeta.
   * rhobeta has a beta hyperprior with parameters rhobeta0 and rhobeta1
   */
  double rho; 
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
  
#if POPADMIXSAMPLER == 1 //DARS sampler
  double AlphaParameters[5];
  AdaptiveRejection** DirParamArray;
#elif POPADMIXSAMPLER == 2//DirichletParamSampler
  double eta;
  double *mu;
  unsigned int obs;
  DirichletParamSampler PopAdmixSampler;
  int *SumLocusAncestry;

#elif POPADMIXSAMPLER == 3 //Hamiltonian sampler
  AlphaSamplerArgs AlphaArgs;
  double *logalpha;
  HamiltonianMonteCarlo AlphaSampler;
  double initialAlphaStepsize;
  float targetAlphaAcceptRate;
#endif


  double *poptheta;    //ergodic average of population admixture, used to centre the values of individual admixture 
                       //in the regression model

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  const Genome* Loci; 
  LogWriter *Log; 

// these methods are 'private static'
  // i.e. private helper methods that cannot
  // access any object variables, nor be used
  // outside of Latent.cc
  //
#if POPADMIXSAMPLER == 1
  static double logf( double, const void* const );
  
  static double  dlogf( double, const void* const );
  
  static double  ddlogf( double, const void* const );
#elif POPADMIXSAMPLER == 3
  static double findE(const double* const theta, const void* const args);
  static void gradE(const double* const theta, const void* const args, double *g);
#endif  
  
  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */

