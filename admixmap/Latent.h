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
#include "gaussian_d.h"
#include "MatrixArray.h"
#include "MatrixArray_d.h"
#include "MatrixArray_i.h"
#include "VectorLoop.h"

#include "TuneRW.h"
#include "MetropolisHastings.h"
#include "HMM.h"
#include "AdmixOptions.h"

#include "DARS.h"

#include "Individual.h"
#include "IndividualCollection.h"

#include "Genome.h"
#include "Chromosome.h"
#include "CompositeLocus.h"
#include "LocusVisitor.h"
#include "AlleleFreqOutputter.h"

#include "chib.h"

#include "LogWriter.h"


class InputData;

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
  //Vector_d *getSumLogTheta();

private:
  //Initialise loop variables   
  void PreUpdate(IndividualCollection *);
  
  static double
  sampleForRho(Vector_d&,DARS*,MatrixArray_i&,MatrixArray_d&);

  void
  OpenOutputFiles();

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
  double SumRho; //

  std::vector<Vector_d> alpha; //population admixture
  Vector_d SumAlpha; //

  //std::vector<bool> _admixed; //population
  //bool _symmetric;
  //Vector_d f; //population, summary of correlation in ancestry between loci
  //Vector_d poptheta; //population

  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 

  LogWriter *Log; 

  //These next members were previously declared in Latent::Sampler
  Vector_d AlphaParameters;
  Vector_d RhoParameters;

  MatrixArray_i rhodata_i;
  MatrixArray_d rhodata_d;
 
  DARS** DirParamArray;
  DARS* RhoDraw;


// these methods are 'private static'
  // i.e. private helper methods that cannot
  // access any object variables, nor be used
  // outside of Latent.cc
  //
  // these could be little calculations,
  // but might also be methods that are going to
  // be moved out of Latent - i.e. making the
  // method static allows us to see what arguments
  // it needs outside this class
  //
  // mostly, the methods at the top of this list
  // are small calculations, whereas the methods at the
  // bottom of the list are on their way out of this class
  // NB some have already been moved
  static double
  logf( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  dlogf( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  ddlogf( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  frho( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  dfrho( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  ddfrho( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  strangExp( double );
  
//   static void
//   s2c( char*, std::string );


  // UNIMPLEMENTED
  // to avoid use
  Latent();
  Latent(const Latent&);
  Latent& operator=(const Latent&);

};

#endif /* !defined LATENT_H */
