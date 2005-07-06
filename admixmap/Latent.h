// *-*-C++-*-*
#ifndef LATENT_H
#define LATENT_H 1

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
#include "DirichletParamSampler.h"

class InputData;
class IndividualCollection;

class Latent
{
public:
  Latent( AdmixOptions* ,Genome*,LogWriter *);
  
  ~Latent();

  void Initialise(IndividualCollection *,std::ofstream *, std::vector<bool> *_admixed, bool *_symmetric,
		  Vector_d *poptheta, std::string *PopulationLabels);  

  void  InitializeOutputFile(std::string *);//possibly move this into PreUpate

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

private:

  double sampleForRho();

  void OpenOutputFiles();

  bool CheckInitAlpha( Vector_d );
  
  /*
   * rho is the sumofintensities parameter. It has a beta prior with shape and scale parameters rhoalpha and rhobeta.
   * rho beta has a beta hyperprior with parameters rhobeta0 and rhobeta1
   */
  double rho; 
  double rhoalpha;
  double rhobeta; 
  double rhobeta0;
  double rhobeta1;
  double SumRho;

   double eta;
   double *mu;
   double *sumlogtheta;
   unsigned int obs;
   DirichletParamSampler PopAdmixSampler;
   std::vector<Vector_d> alpha; //population admixture
  Vector_d SumAlpha; //

  //std::vector<bool> _admixed; //population
  //bool _symmetric;
  //Vector_d poptheta; //population

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 

  LogWriter *Log; 

  //These next members were previously declared in Latent::Sampler
  double RhoParameters[4];

  Matrix_i rhodata_i;
  Matrix_d rhodata_d;
 
  DARS* RhoDraw;

// these methods are 'private static'
  // i.e. private helper methods that cannot
  // access any object variables, nor be used
  // outside of Latent.cc
  //
  
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
