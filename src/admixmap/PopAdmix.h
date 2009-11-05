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
#include "bclib/DelimitedFileWriter.h"
#include "bclib/cvector.h"

class InputData;
class AdmixIndividualCollection;


/** \addtogroup admixmap
 * @{ */


///Class to hold and update population admixture and sumintensities parameters and their priors
class PopAdmix
{
public:
  typedef genepi::cvector<double> PopThetaType; // Must match PedBase::PopThetaType -- break this out from both classes.

  PopAdmix(const AdmixOptions& op, Genome& loci);
  
  ~PopAdmix();
  
  void Initialise(int Numindividuals, const Vector_s&  PopulationLabels, bclib::LogWriter& Log);  
  
  void InitializeOutputFile(const Vector_s&  PopulationLabels);
  
  void UpdateGlobalSumIntensities( const AdmixIndividualCollection & IC, bool sumlogtheta );

  void UpdatePopAdmixParams(int iteration, const AdmixIndividualCollection* const, bclib::LogWriter &Log);
  
  void OutputParams();
  void OutputParams(bclib::Delimitedostream& out); 
  
  void OutputErgodicAvg( int, std::ofstream *avgstream);
  
  const genepi::cvector<double> & getalpha0() const { return alpha[0]; }
  const genepi::cvector<genepi::cvector<double> > & getalpha() const { return alpha; }
  double getrhoalpha() const { return rhoalpha; }
  double getrhobeta() const { return rhobeta; }
  double getglobalrho() const { return rho[0]; }
  const genepi::RhoType & getrho() const { return rho; }

  const genepi::RhoType & getSumLogRho() const { return SumLogRho; }
  const PopThetaType & getpoptheta() const { return poptheta; }
  
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

  genepi::RhoType rho;
  genepi::RhoType rhoproposal;
  double rhoalpha;
  double rhobeta;
  double rhobeta0;
  double rhobeta1;
  genepi::RhoType SumLogRho; //ergodic sum of log(rho)

  //RWM sampler for global rho
  bclib::StepSizeTuner TuneRhoSampler;
  int w, NumberOfUpdates;
  double step;
  double step0;

  genepi::cvector<genepi::cvector<double> > alpha; //population admixture Dirichlet parameters
  genepi::cvector<double> SumAlpha; //ergodic sums of alphas
  //sampler for alpha
  bclib::DirichletParamSampler PopAdmixSampler;
  PopThetaType poptheta; // ergodic average of population admixture, used to
			 // centre the values of individual admixture in the
			 // regression model

  bclib::DelimitedFileWriter outputstream;//output to paramfile


  // UNIMPLEMENTED
  // to avoid use
  PopAdmix();
  PopAdmix(const PopAdmix&);
  PopAdmix& operator=(const PopAdmix&);

};


/** @} */


#endif /* !defined POPADMIX_H */
