/* 
 *   Annealer.cc
 *   Classs to implement simulated annealing and thermodynamic integration
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Annealer.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

Annealer::Annealer(){
  IntervalRatio = 1.03; // size of increments of coolness increases geometrically
  IntervalWidths = 0;
  Coolnesses = 0;
  NumAnnealedRuns = 0;
  Thermo = false;
  _samples = 0;
  _burnin = 0;
  SumEnergy = 0.0;
  SumEnergySq = 0.0;
  LogEvidence = 0.0;
  MeanEnergy = 0.0;
  VarEnergy = 0.0;
  LastMeanEnergy = 0.0;
}

void Annealer::Initialise(bool thermo, unsigned numAnnealedRuns, unsigned samples, unsigned burnin, const char* filename){

  NumAnnealedRuns = numAnnealedRuns; // number of annealed runs excluding last run at coolness of 1
  Thermo = thermo;
  IntervalWidths = 0;
  Coolnesses = 0;
  _samples = samples;
  _burnin = burnin;

  //open output file
  if(Thermo){		
    annealstream.open(filename);
    annealstream << "Coolness\tMeanEnergy\tVarEnergy\tlogEvidence" << std::endl
		 << std::setiosflags(std::ios::fixed) << std::setprecision(3);
  }
}

Annealer::~Annealer(){
  if(annealstream.is_open())annealstream.close();
  delete[] IntervalWidths;
  delete[] Coolnesses;
}

void Annealer::PrintRunLengths(LogWriter& Log, bool testoneindiv){
  int annealedrunlength = _samples;
  if(!Thermo && NumAnnealedRuns > 0) annealedrunlength = 1;
  int finalrunlength = _samples;
  if(Thermo)finalrunlength = _samples*2;//last run is twice as long with thermo option

  if(!testoneindiv && NumAnnealedRuns > 0) {
    Log << On << NumAnnealedRuns << " annealing runs of " << annealedrunlength 
	<< " iteration(s) followed by final run of "; 
  }
  Log << finalrunlength << " iterations at ";
  
  if( testoneindiv) {
    Log << NumAnnealedRuns+1 
	<<" coolnesses for test individual. Other individuals at ";
  }
  Log << "coolness of 1\n";

}

void Annealer::SetAnnealingSchedule(){
  // set annealing schedule
  IntervalWidths = new double[NumAnnealedRuns+1];
  Coolnesses = new double[NumAnnealedRuns + 1];
  Coolnesses[0] = 1.0;//default for unnannealed run
  
  if(NumAnnealedRuns > 0) {
    Coolnesses[0] = 0.0; // change this if you want annealing to start somewhere other than 0;
    // set initial increment of coolness so that geometric series of NumAnnealedRuns increments 
    // will sum to 1 - Coolnesses[0] after NumAnnealedRuns + 1 terms
    IntervalWidths[0] = (1.0 - Coolnesses[0]) * (1.0 - IntervalRatio) /(1.0 - pow(IntervalRatio, NumAnnealedRuns)); 
    //Coolnesses[1] = Coolnesses[0] + IntervalWidths[0];
    if(NumAnnealedRuns > 1) {
      for(int run=1; run <= NumAnnealedRuns; ++run) {
	IntervalWidths[run] = IntervalWidths[run - 1] * IntervalRatio; // geometric increments in interval width
	Coolnesses[run] = Coolnesses[run - 1] + IntervalWidths[run - 1];  
      }
    }
    Coolnesses[NumAnnealedRuns] = 1.0;
  }
}

bool Annealer::SetRunLengths(int run, unsigned* samples, unsigned* burnin, double* coolness){
  bool AnnealedRun = false;

  if(run == NumAnnealedRuns) {
    AnnealedRun = false;
    *coolness = 1.0;
    if(Thermo) {
      *samples = 2*_samples ; // last run is longer
    }
    else{
      *samples = _samples;
    }
    *burnin = _burnin;
  } 
  else {
    if(!Thermo && NumAnnealedRuns > 0) {
    *samples = 1;
    *burnin = 1;
  }

    AnnealedRun = true; 
    *coolness = Coolnesses[run];
  }
  return AnnealedRun;
}

void Annealer::CalculateLogEvidence(int run, double coolness, double SumEnergy, double SumEnergySq, unsigned samples){
  //calculate mean and variance of energy at this coolness
  MeanEnergy = SumEnergy / ((double)samples);
  VarEnergy  = SumEnergySq / ((double)samples) - MeanEnergy * MeanEnergy;
  if(Thermo){// calculate thermodynamic integral  
    //TDI.CalculateLogEvidence(MeanEnergy, VarEnergy );
    annealstream << coolness << "\t" << MeanEnergy << "\t" << VarEnergy;
    if(run > 0) { // use trapezium rule to approximate integral
      LogEvidence -= 0.5*(LastMeanEnergy + MeanEnergy) * IntervalWidths[run];
    } 
    annealstream <<"\t"<< LogEvidence << std::endl; 
    LastMeanEnergy = MeanEnergy;
  }
  std::cout << "\tMeanEnergy = " << MeanEnergy << "        " << std::flush;
}

void Annealer::CalculateLogEvidence(double *SumEnergy, double*SumEnergySq, unsigned size){
  LastMeanEnergy = 0.0;
  double *MeanEner = new double[size];
  double *VarEner = new double[size];
  LogEvidence = 0.0;
  for(unsigned ii = 0; ii < size; ++ii) { // loop over coolnesses to evaluate integral
    //calculate mean and variance of energy at each coolness
    MeanEner[ii] =  SumEnergy[ii] / (double)(_samples - _burnin);
    VarEner[ii] = SumEnergySq[ii] /  (double)(_samples - _burnin) - MeanEner[ii]*MeanEner[ii];
    annealstream << Coolnesses[ii] << "\t" << MeanEner[ii] << "\t" << VarEner[ii];
    // use trapezium rule to approximate integral
    LogEvidence -= 0.5*(LastMeanEnergy + MeanEner[ii]) * IntervalWidths[ii]; 
    annealstream <<"\t"<< LogEvidence << std::endl; 
    LastMeanEnergy = MeanEner[ii];
  }
  MeanEnergy = MeanEner[size-1];//mean at coolness of 1;
  VarEnergy = VarEner[size-1];//var at   ""
  delete[] MeanEner;
  delete[] VarEner;
}

void Annealer::PrintResults(LogWriter& Log, double D_hat){
  double Information = -LogEvidence - MeanEnergy;
  double MeanDeviance = 2.0 * MeanEnergy; 
  double VarDeviance = 4.0 * VarEnergy;
  Log << Quiet << "\nMeanDeviance(D_bar)\t" << MeanDeviance << "\n"
      << "VarDeviance(V)\t" << VarDeviance << "\n"
      << "PritchardStat(D_bar+0.25V)\t" << MeanDeviance + 0.25*VarDeviance << "\n";
  double pD = MeanDeviance - D_hat;
  double DIC = MeanDeviance + pD;
  Log << Quiet << "DevianceAtPosteriorMean(D_hat)\t" << D_hat << "\n"
      << "EffectiveNumParameters(pD)\t" << pD << "\n"
      << "DevianceInformationCriterion\t" << DIC << "\n\n";

  if(Thermo){
    Log << "thermodynamic integration for marginal likelihood yields:\n";
    Log << "LogEvidence " <<  LogEvidence << "\n"; 
    Log << "Information (negative entropy, measured in nats) " << Information << "\n";
  } 
}

const double* Annealer::GetCoolnesses()const{
  return Coolnesses;
}
