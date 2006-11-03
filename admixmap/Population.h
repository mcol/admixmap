// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Population.h (formerly Latent.h)
 *   header file for Population base class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef POPULATION_H
#define POPULATION_H 1

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


class InputData;
class IndividualCollection;

///Class to hold and update population admixture and sumintensities parameters and their priors
class Population
{
public:
  Population( AdmixOptions* op, Genome* loci){
    options = op;
    Loci = loci;
  }

  virtual ~Population(){};
  
  virtual void Initialise(int Numindividuals, const Vector_s&  PopulationLabels, LogWriter& Log) = 0;  
  
  virtual void InitializeOutputFile(const Vector_s&  PopulationLabels) = 0;
  
  virtual void OutputParams(int iteration, LogWriter &Log) = 0;
  virtual void OutputParams(ostream* out) = 0; 
  virtual void OutputErgodicAvg( int, std::ofstream *avgstream) = 0 ;
  
  virtual const std::vector<double> &getalpha0()const = 0;
  virtual const std::vector<std::vector<double> > &getalpha()const = 0;
  virtual double getrhoalpha()const = 0;
  virtual double getrhobeta()const = 0;
  virtual double getglobalrho()const = 0;
  virtual const vector<double>& getrho()const = 0;
  virtual const vector<double>& getSumLogRho()const = 0;
   
  virtual void printAcceptanceRates(LogWriter &Log) = 0;

  //functions specif to PopulationHapMix
  virtual const double* getGlobalTheta()const = 0;
  virtual void SampleHapMixLambda(const int* , bool )  = 0;
    virtual void OutputLambda(const char*) const = 0;

  //functions specific to PopulationAdmixmap
  virtual void UpdateGlobalSumIntensities(const IndividualCollection* const IC, bool sumlogtheta) = 0;
  virtual void UpdatePopAdmixParams(int iteration, const IndividualCollection* const, LogWriter &Log) = 0;
  virtual void resetStepSizeApproximator(int k) = 0; 
  virtual const double *getpoptheta()const = 0;


protected:
  void OpenOutputFiles();
  int K;///< number of subpopulations / block states
  std::ofstream outputstream;//output to paramfile

  AdmixOptions *options;
  Genome* Loci; 
  
private:
  // UNIMPLEMENTED
  // to avoid use
  Population();
  Population(const Population&);
  Population& operator=(const Population&);

};

#endif /* !defined POPULATION_H */

