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
                          //3 = HMCMC


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <algorithm>
#include <gsl/gsl_cdf.h>

#include "rand.h"
#include "vector.h"
#include "vector_i.h"
#include "vector_d.h"
#include "matrix.h"
#include "matrix_i.h"
#include "matrix_d.h"

#include "AdmixOptions.h"
#include "DARS.h"
#include "Individual.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "LogWriter.h"

#if POPADMIXSAMPLER == 2
#include "DirichletParamSampler.h"
#elif POPADMIXSAMPLER == 3
#include "HMCMC.h"
#endif


class InputData;
class IndividualCollection;

class Latent
{
public:
  Latent( AdmixOptions* ,Genome*,LogWriter *);
  
  ~Latent();

  void Initialise(int Numindividuals, std::ofstream *, std::vector<bool> *_admixed, bool *_symmetric,
		  Vector_d *poptheta, std::string *PopulationLabels);  

  void  InitializeOutputFile(std::string *);

  //Updating every iteration
  void Update(int iteration,IndividualCollection *,
	      Vector_d *poptheta,std::ofstream *LogFileStreamPtr);

  void OutputParams(int iteration, std::ofstream *LogFileStreamPtr); 

  void OutputErgodicAvg( int, std::ofstream *avgstream);

  Vector_d *getalpha0();
  std::vector<Vector_d> getalpha();
  double getrhoalpha();
  double getrhobeta();
  double getrho();
  //Vector_d *getSumLogTheta();
  float getAlphaSamplerAcceptanceCount();

private:

  double sampleForRho();

  void OpenOutputFiles();

  bool CheckInitAlpha( Vector_d );
  
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
  //DARS sampler for global rho
  double RhoParameters[4];
  Matrix_i rhodata_i;
  Matrix_d rhodata_d;
  DARS* RhoDraw;
  
  std::vector<Vector_d> alpha; //population admixture Dirichlet parameters
  Vector_d SumAlpha; //ergodic sums of alphas
  //sampler for alpha

#if POPADMIXSAMPLER == 1 //DARS sampler
  double AlphaParameters[5];
  DARS** DirParamArray;
#elif POPADMIXSAMPLER == 2//DirichletParamSampler
  double eta;
  double *mu;
  double *sumlogtheta;
  unsigned int obs;
  DirichletParamSampler PopAdmixSampler;
#elif POPADMIXSAMPLER == 3 //HMCMC sampler
  double *sumlogtheta;
  HMCMC SampleAlpha;
#endif


  //std::vector<bool> _admixed; //population
  //bool _symmetric;
  //Vector_d poptheta; //population

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 
  LogWriter *Log; 

// these methods are 'private static'
  // i.e. private helper methods that cannot
  // access any object variables, nor be used
  // outside of Latent.cc
  //

  static double
  logf( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  static double
  dlogf( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  static double
  ddlogf( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  static double
  frho( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  static double
  dfrho( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  static double
  ddfrho( Vector_d & , Matrix_i&, Matrix_d& , double );
  
  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */

