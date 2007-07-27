// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   PopAdmix.h 
 *   header file for PopAdmix class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef POPADMIX_H
#define POPADMIX_H 1

#include "AdmixOptions.h"
#include "Genome.h"
#include "bclib/StepSizeTuner.h"//for sampling globalrho
#include "bclib/DirichletParamSampler.h"//for sampling pop admix
#include "bclib/DelimitedFIleWriter.h"

class InputData;
class AdmixIndividualCollection;

///Class to hold and update population admixture and sumintensities parameters and their priors
class PopAdmix
{
public:
  PopAdmix(const AdmixOptions& op, Genome& loci);
  
  ~PopAdmix();
  
  void Initialise(int Numindividuals, const Vector_s&  PopulationLabels, bclib::LogWriter& Log);  
  
  void InitializeOutputFile(const Vector_s&  PopulationLabels);
  
  void UpdateGlobalSumIntensities(const AdmixIndividualCollection* const IC, bool sumlogtheta);

  void UpdatePopAdmixParams(int iteration, const AdmixIndividualCollection* const, bclib::LogWriter &Log);
  
  void OutputParams();
  void OutputParams(bclib::Delimitedostream& out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const std::vector<double> &getalpha0()const;
  const std::vector<std::vector<double> > &getalpha()const;
  double getrhoalpha()const;
  double getrhobeta()const;
  double getglobalrho()const;
  const vector<double> &getrho()const;

  const vector<double>& getSumLogRho()const;
  const double *getpoptheta()const;
  
  void printAcceptanceRates(bclib::LogWriter &Log);
  
  void resetStepSizeApproximator(int k); 
  
  //functions required by base class, not used
//   const double* getGlobalTheta()const{return 0;};
//   void SampleHapMixLambda(const int* , bool ) {};
//     void OutputLambda(const char*) const{return;};

private:
  const AdmixOptions& options;
  Genome& Loci; 
  const int K;///< number of subpopulations / block states

  std::vector<double> rho;
  std::vector<double> rhoproposal;
  double rhoalpha;
  double rhobeta;
  double rhobeta0;
  double rhobeta1;
  std::vector<double> SumLogRho; //ergodic sum of log(rho)

  //RWM sampler for global rho
  bclib::StepSizeTuner TuneRhoSampler;
  int w, NumberOfUpdates;
  double step;
  double step0;

  std::vector<std::vector<double> > alpha; //population admixture Dirichlet parameters
  std::vector<double> SumAlpha; //ergodic sums of alphas
  //sampler for alpha
  bclib::DirichletParamSampler PopAdmixSampler;
  double *poptheta;    //ergodic average of population admixture, used to centre the values of individual admixture 
                       //in the regression model

  bclib::DelimitedFileWriter outputstream;//output to paramfile


  // UNIMPLEMENTED
  // to avoid use
  PopAdmix();
  PopAdmix(const PopAdmix&);
  PopAdmix& operator=(const PopAdmix&);

};

#endif /* !defined POPADMIX_H */

