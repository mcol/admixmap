/** 
 *   ADMIXMAP
 *   admixmap.cc 
 *   Top-level source file
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "admixmap.h"
#include <fstream>

using namespace std;

int ReadArgsFromFile(char* filename, int* xargc, char **xargv);
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const string* const PopulationLabels);
void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      const Genome *Loci, Chromosome **Chrm, LogWriter& Log, const std::string* const PopulationLabels, 
		      double coolness, bool anneal);
void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      LogWriter& Log);
void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(); 

void doIterations(const int & samples, const int & burnin, IndividualCollection *IC, Latent & L, AlleleFreqs  & A, Regression *R, 
		  AdmixOptions & options, 
		  const Genome  & Loci, Chromosome **chrm, LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  const double coolness, bool AnnealedRun, ofstream & loglikelihoodfile, 
		  ScoreTests & Scoretest, DispersionTest & DispTest, StratificationTest & StratTest, 
		  MisSpecAlleleFreqTest & AlleleFreqTest, HWTest & HWtest, ofstream & avgstream, InputData & data, const double* Coolnesses);

void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq, std::ofstream *avgstream);

void PrintOptionsMessage();

int main( int argc , char** argv ){
  int    xargc = argc;
  char **xargv = argv;    
  
  if (argc < 2) {
    PrintOptionsMessage();
    exit(1); 
  } else if (argc == 2) {     // using options text file        
    xargc = 1;//NB initialise to 1 to mimic argc (arg 0 is prog name), otherwise first option is ignored later
    xargv = new char*[50];  // change 50 to max number of options
    ReadArgsFromFile(argv[1], &xargc, xargv);        
  }
  
  // ******************* PRIMARY INITIALIZATION ********************************************************************************
  //read user options
  AdmixOptions options(xargc, xargv);
  if(options.getDisplayLevel()>0)
    PrintCopyrightNotice();
  
  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  Log.StartMessage();
  
  smyrand( options.getSeed() );  // Initialise random number seed
  
  InputData data; //read data files and check (except allelefreq files)
  data.readData(&options, Log);//also sets 'numberofoutcomes' and 'populations' options
  
  //check user options
  options.checkOptions(Log, data.getNumberOfIndividuals());
  
  //print user options to args.txt; must be done after all options are set
  options.PrintOptions();
  
  Genome Loci;
  Loci.loadAlleleStatesAndDistances(&options, &data);//reads locusfile and creates CompositeLocus objects
  
  AlleleFreqs A(&Loci);
  A.Initialise(&options, &data, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  
  Chromosome **chrm = 0; //Note: array of pointers to Chromosomes
  chrm = Loci.GetChromosomes(options.getPopulations());  //create Chromosome objects
  Loci.SetSizes(Log);//prints length of genome, num loci, num chromosomes
  
  IndividualCollection *IC = new IndividualCollection(&options, &data, Loci, chrm);//NB call after A Initialise
  IC->LoadData(&options, &data);                             //and before L and R Initialise
  
  Latent L( &options, &Loci);    
  L.Initialise(IC->getSize(), data.GetPopLabels(), Log);
  
  Regression R[2];
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].Initialise(r, IC, Log);
  Regression::OpenOutputFile(&options, IC, data.GetPopLabels(), Log);  
  
  if( options.isGlobalRho() ) {
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ) {
      if( !chrm[j]->isXChromosome() ) chrm[j]->InitialiseLociCorr(L.getrho());
      else chrm[j]->InitialiseLociCorr(0.5*L.getrho());
    }
  }
  cout << flush; 
  IC->Initialise(&options, &Loci, data.GetPopLabels(), L.getalpha(), L.getrhoalpha(), L.getrhobeta(), Log);
  
  //set expected Outcome
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
    R[r].SetExpectedY(IC);
  
  data.Delete();

  //  ******** single individual, one population, fixed allele frequencies  ***************************
  if( IC->getSize() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
    // nothing to do except calculate likelihood
    IC->getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
  else {
    // ******************* INITIALIZE TEST OBJECTS and ergodicaveragefile *******************************
    DispersionTest DispTest;
    StratificationTest StratTest;
    ScoreTests Scoretest;
    MisSpecAlleleFreqTest AlleleFreqTest;
    HWTest HWtest;
    std::ofstream avgstream; //output to ErgodicAverageFile
    
    if( options.getTestForDispersion() ){
      DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
    }
    if( options.getStratificationTest() )
      StratTest.Initialize( &options, Loci, chrm, IC, Log);
    if( options.getScoreTestIndicator() )
      Scoretest.Initialise(&options, IC, &Loci, chrm,data.GetPopLabels(), Log);
    if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Initialise(&options, &Loci, Log );  
    if( options.getHWTestIndicator() )
      HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);
    
    InitializeErgodicAvgFile(&options, IC, Log, &avgstream,data.GetPopLabels());
    
    string s = options.getResultsDir()+"/loglikelihoodfile.txt";
    ofstream loglikelihoodfile(s.c_str());
    
    // ******************* initialise stuff for annealing ************************************************
    double IntervalRatio = 1.03; // size of increments of coolness increases geometrically
    int NumAnnealedRuns = options.getNumAnnealedRuns(); // number of annealed runs excluding last run at coolness of 1
    int samples = 1; // default for annealing runs - overridden for final unannealed run or if AnnealIndicator = true 
    int burnin = 1; // for annealing runs we only want burnin
    
    double SumEnergy = 0.0, SumEnergySq = 0.0, LogEvidence = 0.0;
    double MeanEnergy = 0.0, VarEnergy = 0.0;
    double LastMeanEnergy = 0.0;
    
    double coolness = 1.0; // default
    bool AnnealedRun = false;
    std::ofstream annealstream;//for monitoring energy when annealing
    
    // set annealing schedule
    double *IntervalWidths = 0;
    double *Coolnesses = 0; // 
    IntervalWidths = new double[NumAnnealedRuns + 1];
    Coolnesses = new double[NumAnnealedRuns + 1];

    if(NumAnnealedRuns > 0) {
      // initial increment of coolness from 0 is set so that geometric series of increments will sum to 1 
      // after NumAnnealedRuns additional terms
      IntervalWidths[1] = (IntervalRatio - 1.0) /(pow(IntervalRatio, NumAnnealedRuns) - 1.0 ); 
      Coolnesses[1] = IntervalWidths[1];
      if(NumAnnealedRuns > 1) {
	for(int run=2; run < NumAnnealedRuns+1; ++run) {
	  IntervalWidths[run] = IntervalWidths[run - 1] * IntervalRatio; // geometric increments in interval width
	  Coolnesses[run] = Coolnesses[run - 1] + IntervalWidths[run];  
	}
      }
    }
    Coolnesses[NumAnnealedRuns] = 1.0;
    
    if(options.getThermoIndicator()) { // set up output for thermodynamic integration
      string s = options.getResultsDir()+"/annealmon.txt";
      annealstream.open(s.c_str());
      annealstream << "Coolness\tMeanEnergy\tVarEnergy\tlogEvidence" << endl;
      samples = options.getTotalSamples();
      burnin = options.getBurnIn();
    }
    Log.setDisplayMode(On);
    if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
    if(NumAnnealedRuns > 0) {
      Log << NumAnnealedRuns << " annealing runs of " << samples 
	  << " iteration(s) followed by final run of "; 
    }
    Log << options.getTotalSamples() << " iterations at ";
    if( options.getTestOneIndivIndicator() )Log << options.getNumAnnealedRuns()+1 
						<<" coolnesses for test individual. Other individuals at ";
    Log << "coolness of 1\n";
    
    //Write initial values
    if(options.getIndAdmixHierIndicator()  ){
      if(options.getDisplayLevel()>2)Log.setDisplayMode(On);
      else Log.setDisplayMode(Quiet);
      //Log << "InitialParameterValues:\n"
      //OutputParameters(-1, IC, &L, &A, R, &options, Log);
      //Log << "\n";
    }
    try {
      if(!options.getTestOneIndivIndicator()) {  
	for(int run=0; run < NumAnnealedRuns + 1; ++run) { //loop over coolnesses from 0 to 1
	  // should call a posterior mode-finding algorithm before last run at coolness of 1
	  //resets for start of each run
	  SumEnergy = 0.0;//cumulative sum of modified loglikelihood
	  SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood

	  if(run == NumAnnealedRuns) {
	    AnnealedRun = false;
	    samples = options.getTotalSamples();//redundant?
	    burnin = options.getBurnIn();
	  } else AnnealedRun = true; 
	  coolness = Coolnesses[run];
	  if(NumAnnealedRuns > 0) {
	    cout <<"\rSampling at coolness of " << coolness << "         " << flush;
	    IC->resetStepSizeApproximators(NumAnnealedRuns); // reset k <= NumAnnealedRuns in step size tuners
	  }
	  // accumulate scalars SumEnergy and SumEnergySq at this coolness
	  // array Coolnesses is not used unless TestOneIndivIndicator is true
	  doIterations(samples, burnin, IC, L, A, R, options, Loci, chrm, Log, SumEnergy, SumEnergySq, coolness, AnnealedRun, 
		       loglikelihoodfile, Scoretest, DispTest, StratTest, AlleleFreqTest, HWtest, avgstream, data, Coolnesses);
	  
	  //calculate mean and variance of energy at this coolness
	  MeanEnergy = SumEnergy / ((double)options.getTotalSamples() - options.getBurnIn());
	  VarEnergy  = SumEnergySq / ((double)options.getTotalSamples() - options.getBurnIn()) - MeanEnergy * MeanEnergy;
	  if(options.getThermoIndicator()){// calculate thermodynamic integral
	    annealstream << coolness << "\t" << MeanEnergy << "\t" << VarEnergy;
	    if(run > 0) { // use trapezium rule to approximate integral
	      LogEvidence -= 0.5*(LastMeanEnergy + MeanEnergy) * IntervalWidths[run];
	    } 
	    annealstream <<"\t"<< LogEvidence << endl; 
	    LastMeanEnergy = MeanEnergy;
	  } 
	} // end loop over coolnesses
      } else { // evaluate energy for test individual only at all coolnesses simultaneously
	// call with argument AnnealedRun false - copies of test individual will be annealed anyway  
	doIterations(samples, burnin, IC, L, A, R, options, Loci, chrm, Log, SumEnergy, SumEnergySq, 1.0, false, 
		     loglikelihoodfile, Scoretest, DispTest, StratTest, AlleleFreqTest, HWtest, avgstream, data, Coolnesses);
	// arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
	double *MeanEner = IC->getSumEnergy(); 
	double *VarEner = IC->getSumEnergySq();
	
	double LastMeanEnergy = 0.0; 
	for(int ii = 0; ii < options.getNumAnnealedRuns()+1; ++ii) { // loop over coolnesses to evaluate integral
	  //calculate mean and variance of energy at each coolness
	  MeanEner[ii] /=  ((double)options.getTotalSamples() - options.getBurnIn());
	  VarEner[ii] = VarEner[ii] /  ((double)options.getTotalSamples() - options.getBurnIn()) - MeanEner[ii]*MeanEner[ii];
	  annealstream << Coolnesses[ii] << "\t" << MeanEner[ii] << "\t" << VarEner[ii];
	  // use trapezium rule to approximate integral
	  LogEvidence -= 0.5*(LastMeanEnergy + MeanEner[ii]) * IntervalWidths[ii]; 
	  annealstream <<"\t"<< LogEvidence << endl; 
	  LastMeanEnergy = MeanEner[ii];
	}
	MeanEnergy = MeanEner[options.getNumAnnealedRuns()];//mean at coolness of 1;
	VarEnergy = VarEner[options.getNumAnnealedRuns()];//var at   ""
      } // end evaluation of test individual
    } catch (string msg) {
      Log << msg;
    }
  
    delete[] IntervalWidths;
    delete[] Coolnesses;
    cout << "\nIterations completed                       \n" << flush;
    
    // *************************** OUTPUT AT END ***********************************************************
    if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
    else Log.setDisplayMode(On);
    if( options.getChibIndicator()) {
      IC->OutputChibEstimates(options.isRandomMatingModel(), Log, options.getPopulations());
      //MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
      if(IC->getSize()==1) IC->OutputChibResults(Log);
    }

    double Information = -LogEvidence - MeanEnergy;
    double MeanDeviance = 2.0 * MeanEnergy; 
    double VarDeviance = 4.0 * VarEnergy;
    Log << "\nMeanDeviance(D_bar)\t" << MeanDeviance << "\n"
      << "VarDeviance(V)\t" << VarDeviance << "\n"
      << "PritchardStat(D_bar+0.25V)\t" << MeanDeviance + 0.25*VarDeviance << "\n";
    double D_hat = IC->getDevianceAtPosteriorMean(&options, chrm, R, Log, L.getSumLogRho(), Loci.GetNumberOfChromosomes());
    double pD = MeanDeviance - D_hat;
    double DIC = MeanDeviance + pD;
    Log << "DevianceAtPosteriorMean(D_hat)\t" << D_hat << "\n"
      << "EffectiveNumParameters(pD)\t" << pD << "\n"
      << "DevianceInformationCriterion\t" << DIC << "\n\n"; 

    if(options.getThermoIndicator()){
      Log << "thermodynamic integration for marginal likelihood yields:\n";
      Log << "LogEvidence: " <<  LogEvidence << "\n"; 
      Log << "Information (negative entropy, measured in nats): " << Information << "\n";
    }

    //Residuals
    if(options.getNumberOfOutcomes() > 0)
      IC->OutputResiduals(options.getResidualFilename(), data.getOutcomeLabels(), options.getTotalSamples()-options.getBurnIn());
    //FST
    if( strlen( options.getHistoricalAlleleFreqFilename() ) ){
      A.OutputFST();
    }
    //stratification test
    if( options.getStratificationTest() ) StratTest.Output(Log);
    //dispersion test
    if( options.getTestForDispersion() )  DispTest.Output(options.getTotalSamples() - options.getBurnIn(), Loci, data.GetPopLabels());
    //tests for mis-specified allele frequencies
    if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
      AlleleFreqTest.Output(options.getTotalSamples() - options.getBurnIn(), &Loci, data.GetPopLabels()); 
    //test for H-W eq
    if( options.getHWTestIndicator() )
      HWtest.Output(data.getLocusLabels()); 
    //finish writing score test output as R objects
    if( options.getScoreTestIndicator() ) Scoretest.ROutput();
    
    //output to likelihood ratio file
    if(options.getTestForAffectedsOnly())
      Individual::OutputLikRatios(options.getLikRatioFilename(), options.getTotalSamples()-options.getBurnIn(), data.GetPopLabels());
    
    if(annealstream.is_open())annealstream.close();
    if(avgstream.is_open())avgstream.close();
  }//end else
  cout << "Output to files completed\n" << flush;

  // *************************** CLEAN UP ******************************************************  
  for(unsigned i = 0; i < Loci.GetNumberOfChromosomes(); i++){
    delete chrm[i];
  }
  A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), data.GetPopLabels());
  delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
  delete []chrm;
  
  // ******************* acceptance rates - output to screen and log ***************************
  if( options.getIndAdmixHierIndicator() ){
    if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
    else Log.setDisplayMode(On);
    if(options.getPopulations() > 1){
      L.printAcceptanceRates(Log);
      
      if( options.isGlobalRho() ){
	Log << "Expected acceptance rate in global sumintensities sampler: "
	    << L.getRhoSamplerAccRate()
	    << "\nwith final step size of "
	    << L.getRhoSamplerStepsize()
	    << "\n";
      }
    }
    if(options.getCorrelatedAlleleFreqs()){
      Log<< "Expected acceptance rates in sampler for allele frequency proportion parameters: \n";
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log << A.getAlphaSamplerAcceptanceRate(i) << " ";
      }
      Log<< "Expected acceptance rate in sampler for allele frequency dispersion parameter: \n";
      Log << A.getEtaRWSamplerAcceptanceRate(0)
	  << "\nwith final step sizes of \n";
      for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	if(Loci(i)->GetNumberOfStates()>2)
	  Log << A.getAlphaSamplerStepsize(i) << " " ;
      }
      Log <<  A.getEtaRWSamplerStepsize(0) << "\n" ;
    }
    
#if ETASAMPLER ==1
    if( strlen( options.getHistoricalAlleleFreqFilename() )){
      Log << "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ";
      for(int k = 0; k < options.getPopulations(); ++k){Log << A.getEtaRWSamplerAcceptanceRate(k)<< " " ;}
      Log << "\nwith final step sizes of ";
      for(int k = 0; k < options.getPopulations(); ++k){Log <<  A.getEtaRWSamplerStepsize(k) << " ";}
      Log << "\n";
    }
#endif
  }
  
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  Log.ProcessingTime();
  if(options.getDisplayLevel()>0){
    cout << "Finished\n";
    for(unsigned i = 0; i < 80; ++i)cout<<"*";
    cout<<endl;
  }
  return 0;
} //end of main

void doIterations(const int & samples, const int & burnin, IndividualCollection *IC, Latent & L, AlleleFreqs & A, 
		  Regression *R, AdmixOptions & options, 
		  const Genome & Loci, Chromosome **chrm, LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  double coolness, bool AnnealedRun, ofstream & loglikelihoodfile, 
		  ScoreTests & Scoretest, DispersionTest & DispTest, StratificationTest & StratTest, 
		  MisSpecAlleleFreqTest & AlleleFreqTest, HWTest & HWtest, ofstream & avgstream, InputData & data, const double* Coolnesses) {
  double Energy = 0.0;
  if(!AnnealedRun) cout << endl;
  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if(iteration > burnin) {
      //accumulate energy as minus loglikelihood, calculated using unnanealed genotype probs
      if( !options.getTestOneIndivIndicator() ) {
	Energy = IC->getEnergy(&options, chrm, R, AnnealedRun); // should store loglikelihood if not AnnealedRun
	SumEnergy += Energy;
	SumEnergySq += Energy*Energy;
	// write to file if not AnnealedRun
	if(!AnnealedRun)loglikelihoodfile << iteration<< "\t" << Energy <<endl;
      } else {  
	IC->accumulateEnergyArrays(&options, chrm);
      }
    }
    if( !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }
    
    // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
    if(AnnealedRun || options.getTestOneIndivIndicator()) IC->annealGenotypeProbs(chrm, Loci.GetNumberOfChromosomes(), coolness, Coolnesses); 
    
    UpdateParameters(iteration, IC, &L, &A, R, &options, &Loci, chrm, Log, data.GetPopLabels(), coolness, AnnealedRun);
    Log.setDisplayMode(Quiet);
    if(!AnnealedRun){
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()) )
	OutputParameters(iteration, IC, &L, &A, R, &options, Log);
      
      // ** set merged haplotypes for allelic association score test 
      if( iteration == burnin && options.getTestForAllelicAssociation() ){
	Scoretest.SetAllelicAssociationTest(L.getalpha0());
      }
      
      //Updates and Output after BurnIn     
      if( iteration > burnin ){
	//dispersion test
	if( options.getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A);
	//stratification test
	if( options.getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus(), 
								 options.getPopulations());
	//score tests
	if( options.getScoreTestIndicator() )
	  Scoretest.Update(R[0].getDispersion());//score tests evaluated for first outcome var only
	//tests for mis-specified allelefreqs
	if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	  AlleleFreqTest.Update(IC, &A, &Loci);
	//test for Hardy-Weinberg eq
	if( options.getHWTestIndicator() )
	  HWtest.Update(IC, chrm, &Loci);

	// output every 'getSampleEvery() * 10' iterations (still after BurnIn)
	if (!(iteration % (options.getSampleEvery() * 10))){    
	  //Ergodic averages
	  Log.setDisplayMode(On);
	  if ( strlen( options.getErgodicAverageFilename() ) ){
	    int samples = iteration - burnin;
	    if( options.getIndAdmixHierIndicator() ){
	      L.OutputErgodicAvg(samples,&avgstream);
	      for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
		R[r].OutputErgodicAvg(samples, &avgstream);
	      A.OutputErgodicAvg(samples, &avgstream);
	    }
	    OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq, &avgstream);
	    if(options.getChibIndicator()) IC->OutputErgodicChib(&avgstream);
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options.getScoreTestIndicator() )  Scoretest.Output(iteration, data.GetPopLabels());
	}//end "if every'*10" block
      }//end "if after BurnIn" block
    } // end "if not AnnealedRun" block
  }// end loop over iterations
}

int ReadArgsFromFile(char* filename, int* xargc, char **xargv){
  int  _maxLineLength=1024;
  ifstream fin(filename);
  std::string str;
  //read in line from file
  while (getline(fin,str,'\n')){

    if(str.length()){    //skip blank lines. **should skip also lines with only whitespace
      str.erase(0,str.find_first_not_of(" \t\n\r"));//trim leading whitespace
      if(str[0]!= '#'){ //ignore lines commented out with #
	if(str.find_first_of("#")<str.length())str.erase(str.find_first_of("#"));//ignore #comments
	//trim remaining whitespace
	str.erase(str.find_last_not_of(" \t\n\r")+1);//trailing whitespace
	if(str.find_first_of(" \t\n\r") <= str.length()){//check for any whitespace left
	  if(str.find_first_of(" \t\n\r") < str.find("="))//check for space before '='
	    str.erase(str.find_first_of(" \t\n\r"),str.find("=")-str.find_first_of(" \t\n\r"));//before '='
	  str.erase(str.find_first_of(" \t\n\r"),str.find_last_of(" \t\n")-str.find_first_of(" \t\n\r")+1);//after '='
	}
	//add line to xargv
	xargv[*xargc]=new char[_maxLineLength];
	strcpy(xargv[*xargc],"--");
	strcat(xargv[*xargc],str.c_str());
	++(*xargc);
      }}}
  fin.close();
  return 1;
}

//this function is here because three different objects have to write to avgstream
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const std::string* const PopulationLabels){
  Log.setDisplayMode(Quiet);
  //Open ErgodicAverageFile  
  if( strlen( options->getErgodicAverageFilename() ) ) {
    avgstream->open( options->getErgodicAverageFilename(), ios::out );
    if( !*avgstream ){
      Log.setDisplayMode(On);
      Log << "ERROR: Couldn't open Ergodic Average file\n";
      //exit( 1 );
    } else {
      Log << "Writing ergodic averages of parameters to "
	  << options->getErgodicAverageFilename() << "\n\n";
    }
    
    // Header line of ergodicaveragefile
    if( options->getIndAdmixHierIndicator() ){
      for( int i = 0; i < options->getPopulations(); i++ ){
	*avgstream << PopulationLabels[i] << "\t";
      }
      if( options->isGlobalRho() ) *avgstream << "sumIntensities\t";
      else *avgstream << "sumIntensities.mean\t";
      
      // Regression parameters
      if( options->getNumberOfOutcomes() > 0 ){
	for(int r = 0; r < individuals->getNumberOfOutcomeVars(); ++r){
	  *avgstream << "intercept\t";
	  if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
	    for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	      *avgstream << individuals->getCovariateLabels(i) << "\t";
	    }
	  }
	  if( !options->getTestForAdmixtureAssociation() ){
	    for( int k = 1; k < options->getPopulations(); k++ ){
	      *avgstream << PopulationLabels[k] << "\t";
	    }
	  }
	  if( individuals->getOutcomeType(r)==0 )//linear regression
	    *avgstream << "precision\t";
	}
      }
      
      // dispersion parameters
      if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	for( int k = 0; k < options->getPopulations(); k++ ){
	  *avgstream << "eta" << k << "\t";
	}
      }
    }
    *avgstream << "MeanDeviance\tVarDeviance\t";
    if(options->getChibIndicator()){// chib calculation
      *avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
    }
    *avgstream << "\n";
  } else {
    Log << "No ergodicaveragefile given\n";
  }
}

void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      const Genome *Loci, Chromosome **Chrm, LogWriter& Log, const std::string* const PopulationLabels,
		      double coolness, bool anneal){
  A->ResetAlleleCounts();
  // ** update global sumintensities conditional on genotype probs and individual admixture proportions
  if((options->getPopulations() > 1) && (IC->getSize() > 1) && 
     options->getIndAdmixHierIndicator() && (Loci->GetLengthOfGenome() + Loci->GetLengthOfXchrm() > 0.0))
    L->UpdateSumIntensities(IC, Chrm); // should leave individuals with HMM probs bad, stored likelihood ok
  // this function also sets ancestry correlations
  
  // ** Update individual-level parameters, sampling locus ancestry states, jump indicators, number of arrivals, 
  // individual admixture and sum-intensities 
  // no need to pass global sumintensities parameter as global ancestry correlations have already been set
  IC->Update(iteration, options, Chrm, A, R, L->getpoptheta(), PopulationLabels, L->getalpha(), //L->getrho(), 
	     L->getrhoalpha(), L->getrhobeta(), //Log, 
	     anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level  
  
  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( A->IsRandom() ) {
    bool thermoSampler = (anneal && options->getThermoIndicator() && !options->getTestOneIndivIndicator());
    A->Update(IC, (iteration > options->getBurnIn() && !anneal), coolness, thermoSampler);
  }
  
  if(A->IsRandom() || anneal) { // even for fixed allele freqs, must reset annealed genotype probs as unnannealed  
    IC->setGenotypeProbs(Chrm, Loci->GetNumberOfChromosomes()); // sets unannealed probs ready for getEnergy
    IC->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
    // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
    // or from update of individual-level parameters otherwise
  
  //update population admixture Dirichlet parameters conditional on individual admixture
  L->UpdatePopAdmixParams(iteration, IC, Log, anneal);
  
  // ** update regression parameters (if regression model) conditional on individual admixture
  for(int r = 0; r < options->getNumberOfOutcomes(); ++r)
    R[r].Update((!anneal && iteration > options->getBurnIn()), IC, coolness);
}

void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, Regression *R, const AdmixOptions *options, 
		      LogWriter& Log){
  // fix so that params can be output to console  
  Log.setDisplayMode(Quiet);
  if(options->getIndAdmixHierIndicator()  ){
    //output population-level parameters only when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1) L->OutputParams(iteration, Log);
    //** dispersion parameter (if dispersion model)
	A->OutputEta(iteration, options, Log);
    // ** regression parameters
    for(int r = 0; r < options->getNumberOfOutcomes(); ++r)
      R[r].Output(iteration, options, Log);
    
    // ** new line in log file but not on screen 
    if( iteration == 0 ) {
      Log.setDisplayMode(Off);
      Log << "\n";
      Log.setDisplayMode(Quiet);
    }
  }
  //if( options->getDisplayLevel()>2 ) cout << endl;
  if( iteration > options->getBurnIn() ){
    // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
    if ( strlen( options->getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
    if(options->getOutputAlleleFreq())A->OutputAlleleFreqs();
  }
  // cout << endl;
}

void WriteIterationNumber(const int iteration, const int width, int displayLevel) {
  if( displayLevel > 2 ) { // display parameters on same line
    cout << setiosflags( ios::fixed );
    cout.width(width );
    cout << "\r"<< iteration << " ";  
  }
  else if( displayLevel > 1 && iteration >0) { // display iteration counter only
    cout << "\rIterations so far: " << iteration;
  }
  cout.flush(); 
}


void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq, std::ofstream *avgstream) {
  double EAvDeviance, EVarDeviance;
  EAvDeviance = 2.0*SumEnergy / (double) samples;//ergodic average of deviance
  EVarDeviance = 4.0 * SumEnergySq / (double)samples - EAvDeviance*EAvDeviance;//ergodic variance of deviance 
  *avgstream << EAvDeviance << " "<< EVarDeviance <<" ";
}
void PrintCopyrightNotice(){
  cout << "\n-------------------------------------------------------" << endl;
  cout << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << "Copyright(c) 2002-2006 " <<endl;
  cout << "David O'Donnell, Clive Hoggart and Paul McKeigue"<<endl;
  cout << "Send any comments or queries to david . odonnell@ucd.ie"<<endl;
  cout << "-------------------------------------------------------"<<endl;
  cout << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl;
  cout << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl;
  cout << "-------------------------------------------------------" << endl;
}

void PrintOptionsMessage() {
  cout << "You must specify an options file or list of arguments on command line\n"
       << "Usage:\n"
       << "1. (not recommended) admixmap --[optionname]=[value] ...\n"
       << "2. admixmap [optionfile], where optionfile is a text file containg a list of user options\n"
       << "3. use a Perl script to call the program with command-line arguments. \nSee sample script supplied with this program.\n"
       << "Consult the manual for a list of user options."
       << endl;
}
