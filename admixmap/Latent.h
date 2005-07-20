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
#define POPADMIXSAMPLER 2 //1 = original DARS sampler, 
                          //2 = DirichletParamSampler, 
                          //3 = HamiltonianMonteCarlo

#define GLOBALRHOSAMPLER 1 //1 = DARS
                           //2 = Random-Walk Metropolis

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

#if GLOBALRHOSAMPLER == 2
#include "TuneRW.h"
#endif
#if POPADMIXSAMPLER==1 || GLOBALRHOSAMPLER == 1
#include "DARS.h"
#endif

#if POPADMIXSAMPLER == 2
#include "DirichletParamSampler.h"
#elif POPADMIXSAMPLER == 3
#include "HamiltonianMonteCarlo.h"
#endif


class InputData;
class IndividualCollection;

class Latent
{
public:
  Latent( AdmixOptions* ,Genome*,LogWriter *);
  
  ~Latent();

  void Initialise(int Numindividuals,
		  std::string *PopulationLabels);  

  void  InitializeOutputFile(std::string *);

  //Updating every iteration
  void UpdateRhoWithRW(IndividualCollection *IC, Chromosome **C);
  void Update(int iteration,IndividualCollection *);

  void OutputParams(int iteration); 

  void OutputErgodicAvg( int, std::ofstream *avgstream);

  std::vector<double> &getalpha0();
  std::vector<std::vector<double> > &getalpha();
  double getrhoalpha();
  double getrhobeta();
  double getrho();
  const double *getpoptheta();
  float getAlphaSamplerAcceptanceRate();
  float getAlphaSamplerStepsize();

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
  double SumRho; //ergodic sum of rho

#if GLOBALRHOSAMPLER == 1
  //DARS sampler for global rho
  double RhoParameters[4];
  Matrix_i rhodata_i;
  Matrix_d rhodata_d;
  DARS* RhoDraw;
#elif GLOBALRHOSAMPLER == 2
  //RWM sampler for global rho
  TuneRW TuneRhoSampler;
  int w, NumberOfUpdates;
  double step, step0;
  int NumberAccepted;
  double SumAcceptanceProb;
#endif  

  std::vector<std::vector<double> > alpha; //population admixture Dirichlet parameters
  std::vector<double> SumAlpha; //ergodic sums of alphas
  //sampler for alpha

#if POPADMIXSAMPLER == 1 //DARS sampler
  double AlphaParameters[5];
  DARS** DirParamArray;
#elif POPADMIXSAMPLER == 2//DirichletParamSampler
  double eta;
  double *mu;
  unsigned int obs;
  DirichletParamSampler PopAdmixSampler;
#elif POPADMIXSAMPLER == 3 //Hamiltonian sampler
  double **AlphaArgs;
  double *logalpha;
  HamiltonianMonteCarlo AlphaSampler;
  double initialAlphaStepsize;
  float targetAlphaAcceptRate;
#endif


  double *poptheta;    //ergodic average of population admixture, used to centre the values of individual admixture 
                       //in the regression model

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 
  LogWriter *Log; 

// these methods are 'private static'
  // i.e. private helper methods that cannot
  // access any object variables, nor be used
  // outside of Latent.cc
  //
#if POPADMIXSAMPLER == 1
  static double logf( const double* , Matrix_i&, Matrix_d& , double );
  
  static double  dlogf( const double* , Matrix_i&, Matrix_d& , double );
  
  static double  ddlogf( const double* , Matrix_i&, Matrix_d& , double );
#elif POPADMIXSAMPLER == 3
  static double findE(unsigned dim, const double* const theta, const double* const*args);
  static void gradE(unsigned dim,const double* const theta, const double* const* args, double *g);
#endif  
  static double frho( const double* , Matrix_i&, Matrix_d& , double );
  
  static double dfrho( const double* , Matrix_i&, Matrix_d& , double );
  
  static double ddfrho( const double* , Matrix_i&, Matrix_d& , double );
  
  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */

