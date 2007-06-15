/* 
 *   Model.cc
 *   Abstract base class for top-level model objects
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "Model.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

Model::Model(){
 AISsumlogz = 0.0;
}

Model::~Model(){
  for(unsigned r = 0; r < R.size(); ++r)delete R[r];
  R.clear();
  delete IC;
  if(avgstream.is_open())avgstream.close();
}

void Model::InitialiseGenome(Genome& G, const Options& options, InputData& data, LogWriter& Log){
  G.Initialise(&data, options.getPopulations(), Log);//reads locusfile and creates CompositeLocus objects
  //print table of loci for R script to read
  string locustable = options.getResultsDir();
  locustable.append("/LocusTable.txt");
  G.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1), data.getUnitOfDistanceAsString());
  locustable.clear();
}

void Model::InitialiseRegressionObjects(Options& options, InputData& data,  LogWriter& Log){
  
  Regression::OpenOutputFile(options.getNumberOfOutcomes(), options.getRegressionOutputFilename(), Log);  
  Regression::OpenExpectedYFile(options.getEYFilename(), Log);

  for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
    //determine regression type and allocate regression objects
    if( data.getOutcomeType(r)== Binary ) R.push_back( new LogisticRegression() );
    else if( data.getOutcomeType(r)== Continuous ) R.push_back( new LinearRegression());
    else if( data.getOutcomeType(r)== CoxData ) R.push_back(new CoxRegression());
    
    if(R[r]->getRegressionType()==Cox)
      R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(),data.getCoxOutcomeVarMatrix(), Log);
    else
      R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(), IC->getOutcomeMatrix(), Log);
    
    R[r]->InitializeOutputFile(data.getCovariateLabels(), options.getNumberOfOutcomes());
  }
}

void Model::Start(Options& options, InputData& data, LogWriter& Log, int NumAnnealedRuns){
  // ******************* Initialize test objects and ergodicaveragefile *******************************
  InitialiseTests(options, data, Log);
  
  //open file to output loglikelihood
  string s = options.getResultsDir()+"/loglikelihoodfile.txt";
  loglikelihoodfile.open(s.c_str());
  //set 3 decimal places of precision
  loglikelihoodfile.setf(ios::fixed); 
  loglikelihoodfile.precision(3);
  
  // ******************* Set annealing schedule ************************************************
  
  s = options.getResultsDir()+"/annealmon.txt";
  A.Initialise(options.getThermoIndicator(), NumAnnealedRuns, options.getTotalSamples(), options.getBurnIn(), s.c_str());
  AISsumlogz = 0.0; //for computing marginal likelihood by Annealed Importance Sampling
  
  A.PrintRunLengths(Log, options.getTestOneIndivIndicator());
  A.SetAnnealingSchedule();
  
  //Write initial values
  //     if(options.getIndAdmixHierIndicator()  ){
  //       if(options.getDisplayLevel()>2)Log.setDisplayMode(On);
  //       else Log.setDisplayMode(Quiet);
  //       //Log << "InitialParameterValues:\n"
  //       //OutputParameters(-1, &L, &A, R, &options, Log);
  //       //Log << "\n";
  //     }
}

void Model::Run(Options& options, InputData& data, LogWriter& Log,  
		int NumAnnealedRuns){
  int samples = options.getTotalSamples();
  int burnin = options.getBurnIn();
  double coolness = 1.0; // default

  Start(options, data, Log, NumAnnealedRuns);

  for(int run=0; run <= NumAnnealedRuns; ++run) { //loop over coolnesses from 0 to 1
    // should call a posterior mode-finding algorithm before last run at coolness of 1
    //resets for start of each run
    double SumEnergy = 0.0;//cumulative sum of modified loglikelihood
    double SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood
    bool AnnealedRun = A.SetRunLengths(run, &samples, &burnin, &coolness);
    
    if(NumAnnealedRuns > 0) {
      cout <<"\rSampling at coolness of " << coolness << "       "<< flush;
      // reset approximation series in step size tuners
      int resetk = NumAnnealedRuns; //   
      if(samples < NumAnnealedRuns) {// samples=1 if annealing without thermo integration
	resetk = 1 + run;
      }
      ResetStepSizeApproximators(resetk);
      
    }
    // accumulate scalars SumEnergy and SumEnergySq at this coolness
    // array Coolnesses is not used unless TestOneIndivIndicator is true
    Iterate(samples, burnin, A.GetCoolnesses(), run, options, data, Log, SumEnergy, SumEnergySq,
	    AnnealedRun);
    
    //calculate mean and variance of energy at this coolness
    A.CalculateLogEvidence(run, coolness, SumEnergy, SumEnergySq, samples - burnin);
    
  } // end loop over coolnesses

  Finish(options, data, Log);    
}

void Model::TestIndivRun(Options& options, InputData& data, LogWriter& Log, 
			 int NumAnnealedRuns){
  double SumEnergy = 0.0;//cumulative sum of modified loglikelihood
  double SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood
  int samples = options.getTotalSamples();
  int burnin = options.getBurnIn();

  Start(options, data, Log, NumAnnealedRuns);

  // call with argument AnnealedRun false - copies of test individual will be annealed anyway  
  Iterate(samples, burnin, A.GetCoolnesses(), NumAnnealedRuns, options, data, Log, 
	  SumEnergy, SumEnergySq, false);
  
  // arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
  A.CalculateLogEvidence(getSumEnergy(), getSumEnergySq(), options.getNumAnnealedRuns());

  Finish(options, data, Log);    
}

void Model::GetEnergy(const double* Coolnesses, unsigned coolness, const Options & options, double & SumEnergy, double & SumEnergySq, 
		      double& AISz, bool AnnealedRun, int iteration){
  //accumulate energy as minus loglikelihood, calculated using unannealed genotype probs

  double Energy = IC->getEnergy(&options, R, AnnealedRun); // should store loglikelihood if not AnnealedRun
  
  SumEnergy += Energy;
  SumEnergySq += Energy*Energy;
  if(coolness>0)AISz += exp((Coolnesses[coolness]-Coolnesses[coolness-1])*(-Energy));
  // write to file if not AnnealedRun
  if(!AnnealedRun){
    loglikelihoodfile << iteration<< "\t" << Energy <<endl;
    if(options.getDisplayLevel()>2 && !(iteration%options.getSampleEvery()))cout << Energy << "\t";
  }
}

void Model::ResetStepSizeApproximators(int resetk){
  IC->resetStepSizeApproximators(resetk); 
  pA->resetStepSizeApproximator(resetk);
}

// double Model::getDevianceAtPosteriorMean(const Options* const options, Genome* Loci, LogWriter& Log){
//   return IC->getDevianceAtPosteriorMean(options, R, Loci, Log, L->getSumLogRho(), Loci->GetNumberOfChromosomes(), &A);
// }


void Model::OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq) {
  double EAvDeviance, EVarDeviance;
  EAvDeviance = 2.0*SumEnergy / (double) samples;//ergodic average of deviance
  EVarDeviance = 4.0 * SumEnergySq / (double)samples - EAvDeviance*EAvDeviance;//ergodic variance of deviance 
  avgstream << EAvDeviance << " "<< EVarDeviance <<" ";
}

/// Output at end
void Model::Finish(Options& options, InputData& data, LogWriter& Log){
  cout<< "\nIterations completed                       \n" << flush;

  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
  else Log.setDisplayMode(Quiet);
  
  Finalize(options, Log, data);
  
  A.PrintResults(Log, getDevianceAtPosteriorMean(&options, Log));
  
  if(options.getThermoIndicator()){
    Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
  }
  
  //Expected Outcome
  if(options.getNumberOfOutcomes() > 0){
    Regression::FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), 
                                         data.getOutcomeLabels());
  }
  
  cout << "Output to files completed\n" << flush;
  
  // ******************* acceptance rates - output to screen and log ***************************
  PrintAcceptanceRates(options, Log);
}
