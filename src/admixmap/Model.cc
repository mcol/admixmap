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
  if(Comms::isFreqSampler()){
    //print table of loci for R script to read
    string locustable = options.getResultsDir();
    locustable.append("/LocusTable.txt");
    G.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1));
    locustable.clear();
  }
}

void Model::InitialiseRegressionObjects(Options& options, InputData& data,  LogWriter& Log){
  if( Comms::isMaster() ){
    Regression::OpenOutputFile(options.getNumberOfOutcomes(), options.getRegressionOutputFilename(), Log);  
    Regression::OpenExpectedYFile(options.getEYFilename(), Log);
  }
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
    //determine regression type and allocate regression objects
    if( data.getOutcomeType(r)== Binary ) R.push_back( new LogisticRegression() );
    else if( data.getOutcomeType(r)== Continuous ) R.push_back( new LinearRegression());
    else if( data.getOutcomeType(r)== CoxData ) R.push_back(new CoxRegression());
    
    if(Comms::isMaster()) {
      if(R[r]->getRegressionType()==Cox)
	R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(),data.getCoxOutcomeVarMatrix(), Log);
      else
	R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(), IC->getOutcomeMatrix(), Log);
    }
    else R[r]->Initialise(r, IC->GetNumCovariates());
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
    
    if(Comms::isMaster()){	
      //calculate mean and variance of energy at this coolness
      A.CalculateLogEvidence(run, coolness, SumEnergy, SumEnergySq, samples - burnin);
    } 
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
  if(Comms::isMaster()){
    // arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
    A.CalculateLogEvidence(getSumEnergy(), getSumEnergySq(), options.getNumAnnealedRuns());
  }

  Finish(options, data, Log);    
}
void Model::Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
		    Options & options, InputData & data, LogWriter& Log, 
		    double & SumEnergy, double & SumEnergySq, 
		    bool AnnealedRun) {
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  //Accumulate Energy
  double Energy = 0.0;
  double AISz = 0.0;
  if(isMaster && !AnnealedRun) cout << endl;

  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if( (isMaster || isWorker) && iteration > burnin) {
      //accumulate energy as minus loglikelihood, calculated using unnanealed genotype probs
      if( !options.getTestOneIndivIndicator() ) {
	Energy = IC->getEnergy(&options, R, AnnealedRun); // should store loglikelihood if not AnnealedRun

 	if(isMaster){
	  SumEnergy += Energy;
	  SumEnergySq += Energy*Energy;
	  if(coolness>0)AISz += exp((Coolnesses[coolness]-Coolnesses[coolness-1])*(-Energy));
	  // write to file if not AnnealedRun
	  if(!AnnealedRun){
	    loglikelihoodfile << iteration<< "\t" << Energy <<endl;
	      if(options.getDisplayLevel()>2 && !(iteration%options.getSampleEvery()))cout << Energy << "\t";
	  }
	}
      } else {  
	IC->accumulateEnergyArrays(&options);
      }
    }

     //Write Iteration Number to screen
    if( isMaster && !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }

    //Sample Parameters
    UpdateParameters(iteration, &options, Log, data.GetPopLabels(), Coolnesses, Coolnesses[coolness], AnnealedRun);
    SubIterate(iteration, burnin, options, data, Log, SumEnergy, SumEnergySq,
	       AnnealedRun);
	
  }// end loop over iterations
  //use Annealed Importance Sampling to calculate marginal likelihood
  if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

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
  if(Comms::isMaster())
    cout<< "\nIterations completed                       \n" << flush;

  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
  else Log.setDisplayMode(Quiet);
  
  Finalize(options, Log, data);
  
  A.PrintResults(Log, getDevianceAtPosteriorMean(&options, Log));
  
  if(options.getThermoIndicator()){
    Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
  }
  
  //Expected Outcome
  if(Comms::isMaster() && options.getNumberOfOutcomes() > 0){
    Regression::FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), 
                                         data.getOutcomeLabels());
  }
  
  if(Comms::isMaster())cout << "Output to files completed\n" << flush;
  
  // ******************* acceptance rates - output to screen and log ***************************
  PrintAcceptanceRates(options, Log);
}
